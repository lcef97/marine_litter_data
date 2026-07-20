##' Preliminary settings -------------------------------------------------------
##' Packages must be already loaded in the working environment.

library(INLA)
library(inlabru)
#library(tidyverse)
#library(ggplot2)
library(magrittr)
#library(lubridate)
library(sf)
library(raster)
library(rlang)
library(scico)
library(patchwork)
library(sp)
#library(geodata)
library(spatstat.geom)


##' Define coordinate reference system
kmproj <- CRS("+proj=utm +zone=33  +ellps=WGS84 +units=km +no_defs")

##' Reference grid
r0 <- raster(xmn= 450, ymn= 4000, xmx = 900,ymx = 4550, resolution = 1,
             crs = kmproj)

##' Load all needed data here
input.dir <- "https://github.com/lcef97/marine_litter_data/tree/main/Input"
library(rvest)

filenames <- xml2::read_html(input.dir) %>%
  rvest::html_elements("a[title]") %>%
  rvest::html_attr("title")

filenames <- filenames[grepl("\\.RData$", filenames)]

files <- paste0("https://raw.githubusercontent.com/lcef97/marine_litter_data/main/Input/", filenames)

for (f in files) {
  cat("\nLOADING:", f, "\n")
  tryCatch({
    load(url(f))
    cat("OK:", f, "\n")
  }, error = function(e) {
    cat("FAILED:", f, "\n")
    print(e)
    stop(e)
  })
}

##' Options for INLA and inlabru ----------------------------------------------#
##' First, run without internal.opt=FALSE. This will
##' impact on model replicability, but it's still a long way ahead
##' before results eligible to be published.

options(INLA.expert = FALSE)
bru_options_set(bru_verbose = 1, debug = TRUE)
options(error = function(e) {
  traceback(4)
  quit(status=1)
})
c.c <- list(dic=TRUE, waic=TRUE, config=TRUE, cpo = TRUE, internal.opt = FALSE)


##' GUESS: this should be how `lonlatproj` was defined
##'
lonlatproj <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")


border.ll <- SpatialPolygons(list(Polygons(list(Polygon(poly3)),"0")),proj4string=lonlatproj)
border <- spTransform(border.ll, kmproj)
##' extent - useful for plots
r0.ext <- raster::extent(r0)
##' Standardise covariates -----------------------------------------------------

sd.u <- sd(df$u)
sd.v <- sd(df$v)
sd.logfe <- sd(df$logfe)
sd.fe <- sd(df$fe)
sd.popRadius <- sd(df$popRadius)

sd.depth <- sd(depth_SPDF@data$depth, na.rm = T)
sd.slope <- sd(slope_SPDF@data$slope, na.rm = T)
sd.dcoast <- sd(dist_coast_SPDF@data$dist_coast, na.rm = T)
sd.driver <- sd(dist_river_SPDF@data$dist_river, na.rm = T)
sd.dharbour <- sd(dist_harbour_SPDF@data$dist_harbour, na.rm = T)

df_scaled <- df %>%
  dplyr::mutate(u = .data$u/sd.u) %>%
  dplyr::mutate(v = .data$v/sd.v) %>%
  dplyr::mutate(logfe = .data$logfe/sd.logfe) %>%
  dplyr::mutate(fe = .data$fe/sd.fe) %>%
  dplyr::mutate(popRadius = .data$popRadius/sd.popRadius) %>%
  dplyr::mutate(depth = .data$depth/sd.depth) %>%
  dplyr::mutate(slope = .data$slope/sd.slope) %>%
  dplyr::mutate(dist_coast= .data$dist_coast/sd.dcoast) %>%
  dplyr::mutate(dist_river= .data$dist_river/sd.driver) %>%
  dplyr::mutate(dist_harbour = .data$dist_harbour/sd.dharbour)


cov.spdf <- list(depth_SPDF, dist_coast_SPDF, dist_harbour_SPDF,
                 dist_river_SPDF,  slope_SPDF)

for(i in seq_along(cov.spdf)){
  tmp_SPDF <- cov.spdf[[i]]
  tmp_SPDF@data[,1] <- tmp_SPDF@data[,1]/sd(tmp_SPDF@data[,1], na.rm=T)
  assign(paste0(names(tmp_SPDF@data)[1], "_SPDF_scaled"), tmp_SPDF)
}

pxl_all_scaled <- pxl_all  %>%
  dplyr::mutate(u = .data$u/sd.u) %>%
  dplyr::mutate(v = .data$v/sd.v) %>%
  dplyr::mutate(logfe = .data$logfe/sd.logfe) %>%
  dplyr::mutate(fe = .data$fe/sd.fe) %>%
  dplyr::mutate(popRadius = .data$popRadius/sd.popRadius) %>%
  dplyr::mutate(depth = .data$depth/sd.depth) %>%
  dplyr::mutate(slope = .data$slope/sd.slope) %>%
  dplyr::mutate(dist_coast= .data$dist_coast/sd.dcoast) %>%
  dplyr::mutate(dist_river= .data$dist_river/sd.driver) %>%
  dplyr::mutate(dist_harbour = .data$dist_harbour/sd.dharbour)




##' Mesh building -----------------.--------------------------------------------


##' The sp object will also be useful later
df.sp <- as(df_scaled, "Spatial")
crs(df.sp) <- kmproj
coord.df <- coordinates(df.sp)
colnames(coord.df) <- c("x","y")

get_xy <- function(df) cbind(df$x, df$y)

mesh <- fmesher::fm_mesh_2d_inla(
  boundary = st_as_sf(border),
  loc= coord.df,
  max.edge = c(13,30),
  min.angle = 25,
  cutoff = 5,
  offset = c(10, 40),
  crs = kmproj)

prior.range <- c(100, .6)
spde_gamma <- inla.spde2.pcmatern(mesh, prior.range = prior.range, #  Pr(practic.range<150 km)=0.5
                                  prior.sigma = c(1, .1))  #  P(sigma>1)=0.5


spde_bin <- inla.spde2.pcmatern(mesh, prior.range = prior.range,  # Pr(practic.range<150 km)=0.5
                                prior.sigma = c(1, .1))  #  P(sigma>1)=0.5



values  <- sort(unique(inla.group(df$depth)))
values2 <- sort(unique(inla.group(df$dist_river)))
values3 <- sort(unique(inla.group(df$dist_coast)))

## As adviced from inlabru vignette ==> this should be
## seriously an unambiguous routine, if the system
## still fails it's worrying.

get_xy <- function(df) cbind(df$x, df$y)

## define the model components --------------------------------------------------

mod.year <- list(theta=list(prior="loggamma",fixed = T, initial = log(0.001)))

##' Joint component. Includes literally everything.

cmp_joint <-    ~   -1 +
  ##' year - specific intercepts
  year_gamma_plast(year, model = "iid", hyper=mod.year) +
  year_gamma_Nplast(year, model = "iid", hyper=mod.year) +
  year_bin_plast(year, model = "iid", hyper=mod.year) +
  year_bin_Nplast(year, model = "iid", hyper = mod.year) +
  ##' matern fields
  field_z1(cbind(df_scaled$x, df_scaled$y), model = spde_gamma, group = year,
           control.group = list(model = "iid")) + # yp
  field_z2(cbind(df_scaled$x, df_scaled$y), model = spde_gamma, group = year,
           control.group = list(model = "iid")) + # yo
  field_z3(cbind(df_scaled$x, df_scaled$y), model = spde_gamma, group = year,
           control.group = list(model = "iid")) + # zp
  field_z4(cbind(df_scaled$x, df_scaled$y), model = spde_gamma, group = year,
           control.group = list(model = "iid")) + # zo
  #field_common1(cbind(df_scaled$x, df_scaled$y), copy="field_z1", fixed=F, group=year,
  #              hyper = list(beta = list(prior = "gaussian", param = c(0,10)))) + # yp <> yo
  #field_common2(cbind(df_scaled$x, df_scaled$y),  copy="field_z1", fixed=F, group = year,
  #              hyper = list(beta = list(prior = "gaussian", param = c(0,10)))) + # yp <> zp
  #field_common3(cbind(df_scaled$x, df_scaled$y),  copy="field_z2", fixed=F, group = year,
  #              hyper = list(beta = list(prior = "gaussian", param = c(0,10)))) + # yo <> zp
  #field_common4(cbind(df_scaled$x, df_scaled$y), copy="field_z1", fixed=F, group = year,
  #              hyper = list(beta = list(prior = "gaussian", param = c(0,10)))) + # yp <> zo
  #field_common5(cbind(df_scaled$x, df_scaled$y), copy="field_z2", fixed=F, group = year,
  #              hyper = list(beta = list(prior = "gaussian", param = c(0,10)))) + # yo <> zo
  #field_common6(cbind(df_scaled$x, df_scaled$y), copy="field_z3", fixed=F, group = year,
  #              hyper = list(beta = list(prior = "gaussian", param = c(0,10)))) + # zp <> z0
  ##' nonlinear covariates for density
  depth_gamma_plast(depth_SPDF,  model = "rw2", main_layer = "depth",
                    values = values, scale.model = TRUE) +
  depth_gamma_Nplast(depth_SPDF,  model = "rw2", main_layer = "depth",
                     values = values, scale.model = TRUE) +
  ##' linear covariates
  driver_gamma_plast(dist_river_SPDF_scaled,main_layer =   "dist_river") +
  dcoast_gamma_plast(dist_coast_SPDF_scaled,main_layer = "dist_coast") +
  dharbour_gamma_plast(dist_harbour_SPDF_scaled,main_layer = "dist_harbour") +
  #slope_gamma_plast(slope_SPDF_scaled,main_layer = "slope") +
  u_gamma_plast(u,main_layer = "u") +
  v_gamma_plast(v,main_layer = "v") +
  logfe_gamma_plast(logfe,main_layer = "logfe") +
  pop_radius_gamma_plast(popRadius, main_layer = "popRadius")+
  driver_gamma_Nplast(dist_river_SPDF_scaled,main_layer =   "dist_river") +
  dcoast_gamma_Nplast(dist_coast_SPDF_scaled,main_layer = "dist_coast") +
  dharbour_gamma_Nplast(dist_harbour_SPDF_scaled,main_layer = "dist_harbour") +
  #slope_gamma_Nplast(slope_SPDF_scaled,main_layer = "slope") +
  u_gamma_Nplast(u,main_layer = "u") +
  v_gamma_Nplast(v,main_layer = "v") +
  logfe_gamma_Nplast(logfe,main_layer = "logfe")+
  pop_radius_gamma_Nplast(popRadius, main_layer = "popRadius")+
  #
  depth_bin_plast(depth_SPDF_scaled, main_layer = "depth") +
  driver_bin_plast(dist_river_SPDF_scaled, main_layer = "dist_river") +
  dcoast_bin_plast(dist_coast_SPDF_scaled, main_layer = "dist_coast") +
  dharbour_bin_plast(dist_harbour_SPDF_scaled,main_layer = "dist_harbour") +
  #slope_bin_plast(slope_SPDF_scaled,  main_layer = "slope") +
  u_bin_plast(u,main_layer = "u") +
  v_bin_plast(v,main_layer = "v") +
  logfe_bin_plast(logfe,main_layer = "logfe")+
  pop_radius_bin_plast(popRadius, main_layer = "popRadius")+
  depth_bin_Nplast(depth_SPDF_scaled, main_layer = "depth") +
  driver_bin_Nplast(dist_river_SPDF_scaled, main_layer = "dist_river") +
  dcoast_bin_Nplast(dist_coast_SPDF_scaled,  main_layer = "dist_coast") +
  dharbour_bin_Nplast(dist_harbour_SPDF_scaled,main_layer = "dist_harbour") +
  #slope_bin_Nplast(slope_SPDF_scaled, main_layer = "slope") +
  u_bin_Nplast(u,main_layer = "u") +
  v_bin_Nplast(v,main_layer = "v") +
  logfe_bin_Nplast(logfe,main_layer = "logfe") +
  pop_radius_bin_Nplast(popRadius, main_layer = "popRadius")



##' Formulas - adapted from above sections

formula_gamma_plast  <- y_plast ~
  year_gamma_plast +
  field_z1 +
  depth_gamma_plast +
  driver_gamma_plast +
  dcoast_gamma_plast +
  dharbour_gamma_plast +
  #slope_gamma_plast +
  pop_radius_gamma_plast +
  u_gamma_plast + v_gamma_plast + logfe_gamma_plast



#formula for y_Nplastic - as above
formula_gamma_Nplast  <- y_Nplast ~
  year_gamma_Nplast +
  field_z2 +
  depth_gamma_Nplast +
  driver_gamma_Nplast +
  dcoast_gamma_Nplast +
  dharbour_gamma_Nplast +
  #slope_gamma_Nplast +
  pop_radius_gamma_Nplast +
  u_gamma_Nplast +
  v_gamma_Nplast +
  logfe_gamma_Nplast


#formula for z
formula_bin_plast <- z_plast ~
  year_bin_plast +
  field_z3 +
  depth_bin_plast +
  #slope_bin_plast +
  driver_bin_plast +
  dcoast_bin_plast +
  dharbour_bin_plast +
  u_bin_plast +
  v_bin_plast +
  logfe_bin_plast +
  pop_radius_bin_plast +
  offset(ssa)


formula_bin_Nplast <- z_Nplast ~
  year_bin_Nplast +
  field_z4 +
  depth_bin_Nplast +
  #slope_bin_Nplast +
  driver_bin_Nplast +
  u_bin_Nplast +
  v_bin_Nplast +
  dcoast_bin_Nplast +
  dharbour_bin_Nplast +
  logfe_bin_Nplast +
  pop_radius_bin_Nplast +
  offset(ssa)



##' Likelihoods. Likewise, defined as in previous parts

lik_gamma_plast <- bru_obs("gamma",
                           formula = formula_gamma_plast,
                           samplers = border,
                           domain = list(coordinates = mesh),
                           data = df_scaled)
lik_gamma_Nplast <- bru_obs("gamma",
                            formula = formula_gamma_Nplast,
                            samplers = border,
                            domain = list(coordinates = mesh),
                            data = df_scaled)
lik_bin_plast <- bru_obs("binomial",
                         formula = formula_bin_plast,
                         samplers = border,
                         domain = list(coordinates = mesh),
                         data = df_scaled)
lik_bin_Nplast <- bru_obs("binomial",
                          formula = formula_bin_Nplast,
                          samplers = border,
                          domain = list(coordinates = mesh),
                          data = df_scaled)


##' NOW comes the final boss ---------------------------------------------------
##'



withCallingHandlers({
  fit_marg0_ccd <-bru(
    cmp_joint,  lik_gamma_plast, lik_gamma_Nplast,
    lik_bin_plast, lik_bin_Nplast,
    options = list(
      control.predictor=list(link = 1),
      control.compute = c.c,
      bru_max_iter=1, verbose = T, debug = T,
      num.threads = 1))
}, error = function(e) {
  message("!!!!!!!!!!!!!!! \n!!! WARNING !!! \n!!! Author-added message: either inla or inlabru literally failed, traceback follows")
  if (file.exists("bru_crash_dump_marg0_ccd.rda")) {
    message("removing extant crash_marg0_ccd")
    file.remove("bru_crash_dump_marg0_ccd.rda")
  }
  dump.frames(dumpto = "bru_crash_dump_marg0_ccd", to.file = TRUE)
})

cat("configs length:", length(fit_marg0_ccd$misc$configs), "\n")


#' Be cautious here ==> rename according to covariates used
filename.marg0 <- paste0("fit_marg0_ccd_alldists", lubridate::today(), ".RData")

filepath.rm <- file.path(getwd(), filename.marg0)

if (file.exists(filepath.rm)) {
  file.remove(filepath.rm)
}


lgocv_marg0 <- tryCatch({
  list(lgocv_marg0_g1 = inla.group.cv(fit_marg0_ccd, num.level.sets=1),
       lgocv_marg0_g3 = inla.group.cv(fit_marg0_ccd, num.level.sets=3),
       lgocv_marg0_g5 = inla.group.cv(fit_marg0_ccd, num.level.sets=5))
}, error = function(e) {
  message("!!! ERROR !!! \n impossible to get lgocv !!!")
  message(conditionMessage(e))
  dump.frames(
    to.file = TRUE,
    dumpto = "groupcv_marg0_dump"
  )
  traceback()
  print(e)
  return(NULL)
})


#' Be cautious here ==> rename according to covariates used
fit_marg0_ccd_alldists <- list(fit_marg0_ccd, lgocv_marg0)
save(fit_marg0_ccd_alldists, file = filename.marg0)