# library -----------------------------------------------------------------
set.seed(123)
library(sf)
library(osmdata)
library(ggplot2)
library(mapview)
# Reorder OpenStreetMap as default
all_basemaps <- c("OpenStreetMap", 
                  "CartoDB.Positron",  
                  "CartoDB.DarkMatter", 
                  "Esri.WorldImagery", 
                  "Esri.WorldStreetMap") 
# Set default basemap to OSM
mapviewOptions(basemaps = all_basemaps)
mapviewOptions(fgb = FALSE)
library(rnaturalearth)
library(tidyverse)
library(terra)
# -------------------------
# a. Madrid map (the city) -----------------------------------------------------------------------
# -------------------------
map <- geodata::gadm(country = "ESP", level = 4, path = tempdir())
map <- map |> 
  st_as_sf() |> 
  filter(NAME_4 == "Madrid")
mapview(map)
view(map)
# -------------------------
# b. Covariate 1: Elevation -----------------------------------------------------------------------
# -------------------------

# 2. Get covariate 1: elevation for Spain (SpatRaster)
r <- geodata::elevation_30s(country = "ESP", path = tempdir())

# 3. Convert sf polygon to SpatVector (terra format)
madrid_vect <- vect(map)

# 4. Crop and mask the raster to Madrid
madrid_crop <- crop(r, madrid_vect)
elev <- mask(madrid_crop, madrid_vect)

# 5. Plot
plot(elev)

# -------------------------
# c. Covariate 2: Roads data (OSM) -----------------------------------------------------------------------
# -------------------------
# We'll request OSM "highway" features (this includes major roads, minor roads, residential, etc.)
bb <- st_bbox(map)

# 2) Define the highway groups (strings for OSM)
major_vals     <- c("motorway","motorway_link","trunk","trunk_link","primary","primary_link")
secondary_vals <- c("secondary","secondary_link","tertiary","tertiary_link")
local_vals     <- c("residential","living_street","service","unclassified")

# Combined values you want to download
all_vals <- c(major_vals, secondary_vals, local_vals)

# 3) Query OSM (only these highway types)
q1 <- opq(bbox = bb) %>%
  add_osm_feature(key = "highway", value = major_vals)
q2 <- opq(bbox = bb) %>%
  add_osm_feature(key = "highway", value = secondary_vals)
q3 <- opq(bbox = bb) %>%
  add_osm_feature(key = "highway", value = local_vals)

# 4) Get data
osm_roads1 <- osmdata_sf(q1)
osm_roads2 <- osmdata_sf(q2)
osm_roads3 <- osmdata_sf(q3)

roads_sf1 <-
  osm_roads1$osm_lines %>%
  st_transform(25830) %>%
  select(name)
roads_sf2 <-
  osm_roads2$osm_lines %>%
  st_transform(25830) %>%
  select(name)
roads_sf3 <-
  osm_roads3$osm_lines %>%
  st_transform(25830) %>%
  select(name)

roads_major_sf <- st_union(roads_sf1$geometry)
roads_secondary_sf <- st_union(roads_sf2$geometry)
roads_local_sf <- st_union(roads_sf3$geometry)

# -------------------------
# c. Covariate 3: NDVII  -----------------------------------------------------------------------
# -------------------------

ndvi <- terra::rast("data/NDVI/NDVI_20200601T000000Z.tif")
ndvi <- terra::crop(ndvi, map) |> mask(map)

# ----------------
---------
# d. Response Variable - NOISE -----------------------------------------------------------------------
# -------------------------

# station location
library(readr)
station <- read_delim(
  "data/station.csv",
  delim = ";",
  locale = locale(encoding = "ISO-8859-1"),
  trim_ws = TRUE
)
# ETRS89 geographic coords
station <- station |> 
  # as_tibble() |> 
  st_as_sf(coords = c("X", "Y"), crs=25830)
  # |>  st_transform(crs=4326)

# noise data
monthly <- read_delim(
  "data/monthly.csv",
  delim = ";",
  locale = locale(encoding = "ISO-8859-1", decimal_mark = ","),
  na = c("N/D"),
  trim_ws = TRUE
)

yearly <-
  monthly |> 
  group_by(Estación, Año) |> 
  summarise(LAeq = mean(LAeq, na.rm=T))

noise24 <- filter(yearly, Año == 2024)

# -------------------------
# e. Join all data -----------------------------------------------------------------------
# -------------------------
y <- left_join(noise24, station, join_by(Estación == ESTACIÓN)) |> 
  select(Estación, NOMBRE, geometry, LAeq) |> 
  st_as_sf(crs=25830)

y$elev <- terra::extract(elev, vect(y))$ESP_elv_msk
y$dist_major <- as.numeric(st_distance(y, roads_major_sf))
y$dist_secondary <- as.numeric(st_distance(y, roads_secondary_sf))
y$dist_local <- as.numeric(st_distance(y, roads_local_sf))
y$ndvi <- terra::extract(ndvi, vect(y))$NDVI_20200601T000000Z

# -------------------------
# e. Prediction points -----------------------------------------------------------------------
# -------------------------
yp <- terra::crds(elev)
yp <- as.data.frame(yp)
dim(yp)
yp_sf <- yp |> 
  st_as_sf(coords = c("x","y"), crs=4326) |> 
  st_transform(25830)
yp$elev <- terra::extract(elev, yp)$ESP_elv_msk  
yp$dist_major <- as.numeric(st_distance(yp_sf, roads_major_sf))
yp$dist_secondary <- as.numeric(st_distance(yp_sf, roads_secondary_sf))
yp$dist_local <- as.numeric(st_distance(yp_sf, roads_local_sf))
yp$ndvi <- terra::extract(ndvi, yp[,1:2])$NDVI_20200601T000000Z

# save(y, yp,
#      roads_major_sf, roads_secondary_sf, roads_local_sf,
#      file = "data/data.rda")
# save(map, file = "data/map.rda")
# writeRaster(elev, "data/elev.tif", overwrite = TRUE)

# -------------------------
# 0. Load data -----------------------------------------------------------------------
# -------------------------
load("data/map.rda")
load("data/data.rda")
elev <- terra::rast("data/elev.tif")
ndvi <- terra::rast("data/NDVI/NDVI_20200601T000000Z.tif")
ndvi <- terra::crop(ndvi, map) |> mask(map)


# -------------------------
# Mesh -----------------------------------------------------------------------
# -------------------------
library(INLA)
# cropped boundary map
t <- inla.nonconvex.hull(y, convex= -0.13)
plot(t)
points(y$geometry)

# crop prediction points
library(fmesher)
t_sf <- fm_as_sfc(t) |>
  st_set_crs(25830) |> 
  st_transform(crs=4326)
st_crs(t_sf)
# elev_c <- terra::crop(elev, t_sf) |> mask(vect(t_sf))
# coop <- terra::crds(elev) need to crop
yp_sf <- st_as_sf(yp, coords = c("x","y"), crs=4326)
yp2 <- st_filter(yp_sf, t_sf)
nrow(yp2)

# mesh
coo <- st_coordinates(y$geometry)
mesh <- inla.mesh.2d(
  boundary = t, # remove this if want whole of madrid
  loc = coo, max.edge = c(500, 2000),
  cutoff = 100
)
mesh$n
plot(mesh)
points(coo, col="red")

spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

# Projection matrix
A <- inla.spde.make.A(mesh = mesh, loc = coo)
dim(A)

# Prediction matrix
coop <- st_coordinates(yp2$geometry)
# coop <- yp[,1:2] |> as.matrix() # if want whole of madrid
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
dim(Ap)

# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est",
  data = list(y = y$LAeq), 
  A = list(1, A),
  effects = list(data.frame(b0 = 1, 
                            elev = y$elev,
                            ndvi = y$ndvi,
                            dist_major = y$dist_major,
                            dist_secondary = y$dist_secondary,
                            dist_local = y$dist_local), 
                 s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = 1, 
                            elev = yp2$elev,
                            ndvi = yp2$ndvi,
                            dist_major = yp2$dist_major,
                            dist_secondary = yp2$dist_secondary,
                            dist_local = yp2$dist_local),
                 s = indexs
  )
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

# -------------------------
# 1. EDA -----------------------------------------------------------------------
# -------------------------
library(RColorBrewer)
pal <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
mapview(y, zcol = "LAeq", col.regions = pal, direction=-1, layer.name="LAeq (dB)") +
  mapview(map, alpha.regions=0.4)

library(tmap)
# ensure projected CRS (important for distance-based polygons)
y_proj <- st_transform(y, 4326)  # change EPSG if needed
map_proj <- st_transform(t_sf, 32632)

# convert to terra vectors
y_vect <- terra::vect(y_proj)
map_vect <- terra::vect(t_sf)

# make Voronoi polygons clipped to study area
v <- terra::voronoi(y_vect, bnd = map_vect)

v <- st_as_sf(v)
tm_shape(v) +
  tm_fill(col = "LAeq", palette = "magma", alpha = 0.6)

# -------------------------
# 2. Covariates Plot -----------------------------------------------------------------------
# -------------------------
elev_c <- terra::crop(elev, vect(t_sf)) |> mask(vect(t_sf))
mapview(elev_c, alpha = 0.6, layer.name="elevation")
ndvi <- terra::crop(ndvi, vect(t_sf)) |> mask(vect(t_sf))
mapview(ndvi, alpha = 0.5, layer.name="NDVI")
mapview(yp2, alpha = 0.1, layer.name="Prediction Points") + 
  mapview(roads_major_sf, color = "red") +
  mapview(roads_secondary_sf, color = "blue") +
  mapview(roads_local_sf, color = "lightblue")

# -------------------------
# 2. Model -----------------------------------------------------------------------
# -------------------------
formula <- y ~ 0 + b0 + ndvi + elev + dist_major + dist_secondary + 
  dist_local + f(s, model = spde)
# formula <- y ~ 0 + b0 + ndvi + dist_major + dist_secondary + f(s, model = spde)
res <- inla(formula, family = "gaussian",
            data = inla.stack.data(stk.full),
            control.predictor = list(compute = TRUE,
                                     A = inla.stack.A(stk.full)),
            control.compute = list(return.marginals.predictor = TRUE))
res$summary.fixed
summary(res)

index <- inla.stack.index(stack = stk.full, tag = "pred")$data

noise_mean <- res$summary.fitted.values[index, "mean"]
noise_ll <- res$summary.fitted.values[index, "0.025quant"]
noise_ul <- res$summary.fitted.values[index, "0.975quant"]
noise_sd <- res$summary.fitted.values[index, "sd"]

# check residual 
index_est <- inla.stack.index(stack = stk.full, tag = "est")$data
noise_mean_est <- res$summary.fitted.values[index_est, "mean"]
residual <- y$LAeq - noise_mean_est
residual_sf <- bind_cols(y, residual=residual)
mapview(residual_sf, zcol="residual",
        col.regions = rev(viridis::inferno(33)),
        at = seq(-8,8,0.5))
rmse <- sqrt(mean(residual^2, na.rm = TRUE))
rmse

# Create SpatVector object
sv <- terra::vect(coop, atts = data.frame(noise_mean = noise_mean,
                                          noise_ll = noise_ll, 
                                          noise_ul = noise_ul,
                                          noise_sd = noise_sd),
                  crs = "+proj=longlat +datum=WGS84")

# rasterize
r_noise_mean <- terra::rasterize(
  x = sv, y = elev, field = "noise_mean",
  fun = mean
)
mapview(r_noise_mean, alpha=0.6, layer.name="predicted mean noise")

r_noise_sd <- terra::rasterize(
  x = sv, y = elev, field = "noise_sd",
  fun = mean
)
r_noise_sd2 <- ifel(r_noise_sd >= 0 & r_noise_sd <= 3, r_noise_sd, NA)
mapview(r_noise_sd, alpha=0.6, at = c(0, 1, 2, 3, 40), layer.name="sd")
# -------------------------
# Exceedence Probabiliity -----------------------------------------------------------------------
# -------------------------
excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 62, marginal = marg)})
sv$excprob <- excprob

r_excprob <- terra::rasterize(
  x = sv, y = elev,
  field = "excprob",
  fun = mean
)

pal <- colorNumeric("magma", c(0, 1), na.color = "transparent")
mapview(r_excprob, alpha=0.6, layer.name="P(n>62)") + 
  mapview(y, zcol="LAeq", hide=TRUE, layer.name="Station")

# -------------------------
# Non-Exceedence Probabiliity -----------------------------------------------------------------------
# -------------------------

nonexcprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){inla.pmarginal(q = 56, marginal = marg)})

sv$nonexcprob <- nonexcprob
r_nonexcprob <- terra::rasterize(
  x = sv, y = elev, field = "nonexcprob",
  fun = mean
)

pal <- colorNumeric("magma", c(0, 1), na.color = "transparent")
mapview(r_nonexcprob, alpha=0.6, layer.name="P(n<56)") + 
  mapview(y, zcol="LAeq", hide=TRUE, layer.name="Station") 

