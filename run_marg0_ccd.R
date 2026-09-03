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
input.dir <- "Input"
files <- list.files(input.dir, full.names = TRUE)


cat("\n=== FILESYSTEM DIAGNOSTIC ===\n")

cat("WD:", getwd(), "\n")

cat("Input is dir:", dir.exists(input.dir), "\n")

cat("\nListing Input BEFORE load:\n")
print( files )



for (f in files) {
  cat("\nLOADING:", f, "\n")
  tryCatch({
    load(f)
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
c.c <- list(dic=TRUE, waic=TRUE, config=TRUE, cpo = TRUE, internal.opt = FALSE,
            return.marginals.predictor = TRUE)


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
  field_z1(geometry, model = spde_gamma, group = year,
           control.group = list(model = "iid")) + # yp
  field_z2(geometry, model = spde_gamma, group = year,
           control.group = list(model = "iid")) + # yo
  field_z3(geometry, model = spde_gamma, group = year,
           control.group = list(model = "iid")) + # zp
  field_z4(geometry, model = spde_gamma, group = year,
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
                           domain = list(geometry = mesh),
                           data = df_scaled)
lik_gamma_Nplast <- bru_obs("gamma",
                            formula = formula_gamma_Nplast,
                            samplers = border,
                            domain = list(geometry = mesh),
                            data = df_scaled)
lik_bin_plast <- bru_obs("binomial",
                         formula = formula_bin_plast,
                         samplers = border,
                         domain = list(geometry = mesh),
                         data = df_scaled)
lik_bin_Nplast <- bru_obs("binomial",
                          formula = formula_bin_Nplast,
                          samplers = border,
                          domain = list(geometry = mesh),
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
fit_marg0_ccd_alldists <- list(res=fit_marg0_ccd, cv=lgocv_marg0)
save(fit_marg0_ccd_alldists, file = filename.marg0)




##' Further added things - not in HTC ------------------------------------------

 

library(excursions)
exc <- excursions.inla(
  fit_marg0_ccd,
  ind    = seq_len(1540),
  u      = log(30),
  type   = ">",
  method = "NI", verbose = T)
 

exc_joint <- excursions.inla(
  fit_marg0_ccd,
  ind    = seq_len(1540),
  u      = log(30),
  type   = ">",
  method = "NI", verbose = T)

#mesh.index <- inla.spde.make.index(name = "field", n.spde = spde_gamma$n.spde)
#nxy <- c(100,100)
#projgrid <- inla.mesh.projector(mesh, dims = nxy)

#xy.in <- splancs::inout(projgrid$lattice$loc,poly=border), cbind(PRborder[, 1], PRborder[, 2]))
#submesh = submesh.grid(matrix(xy.in, nxy[1], nxy[2]),
#                       list(loc = projgrid$lattice$loc, dims = nxy))


##' Excursion set - excruciating -----------------------------------------------
library(excursions)
library(lemon)
library(gridExtra)

#### Excursions set

lambda_plast  <- generate(
  fit_marg0_ccd, pxl_all_scaled,
  ~ exp(year_gamma_plast + field_z1 +depth_gamma_plast +
          driver_gamma_plast + dcoast_gamma_plast + dharbour_gamma_plast +
          #slope_gamma_plast +
          pop_radius_gamma_plast +u_gamma_plast + v_gamma_plast + logfe_gamma_plast),
  verbose = T)

lambda_Nplast <- generate(
  fit_marg0_ccd, pxl_all_scaled,
  ~ exp(year_gamma_Nplast + field_z2 + +depth_gamma_Nplast +
          driver_gamma_Nplast + dcoast_gamma_Nplast + dharbour_gamma_Nplast +
          #slope_gamma_Nplast +
          pop_radius_gamma_Nplast +u_gamma_Nplast + v_gamma_Nplast + logfe_gamma_Nplast), 
  verbose = T)

lambda_plast_pred  <- predict(
  fit_marg0_ccd, df_scaled,
  ~ exp(year_gamma_plast + field_z1 +depth_gamma_plast +
          driver_gamma_plast + dcoast_gamma_plast + dharbour_gamma_plast +
          #slope_gamma_plast +
          pop_radius_gamma_plast +u_gamma_plast + v_gamma_plast + logfe_gamma_plast),
  verbose = T)

lambda_Nplast_pred <- generate(
  fit_marg0_ccd, df_scaled,
  ~ exp(year_gamma_Nplast + field_z2 + +depth_gamma_Nplast +
          driver_gamma_Nplast + dcoast_gamma_Nplast + dharbour_gamma_Nplast +
          #slope_gamma_Nplast +
          pop_radius_gamma_Nplast +u_gamma_Nplast + v_gamma_Nplast + logfe_gamma_Nplast), 
  verbose = T)




prob_plast <- generate(fit_marg0_ccd, pxl_all_scaled, ~ 
                         probs( year_bin_plast + field_z3 +
                                  depth_bin_plast + driver_bin_plast + 
                                  dcoast_bin_plast + dharbour_bin_plast + u_bin_plast + v_bin_plast + 
                                  logfe_bin_plast + pop_radius_bin_plast  ))


pred <- list(lambda_plast, lambda_Nplast, lambda_plast_pred, lambda_Nplast_pred)




probs <- function(x) exp(x)/(1+exp(x))


p1 = data.frame(x = sf::st_coordinates(pxl_all)[,1],
                y = sf::st_coordinates(pxl_all)[,2],
                year = pxl_all$year,
                lambda_plast * prob_plast,
                z = apply(lambda_plast * prob_plast,1,mean))








aa <- data.frame(x = sf::st_coordinates(pxl_all)[,1],
                y = sf::st_coordinates(pxl_all)[,2],
                year = pxl_all$year,
                lambda_plast)

aa_Nplast <- data.frame(x = sf::st_coordinates(pxl_all)[,1],
                       y = sf::st_coordinates(pxl_all)[,2],
                       year = pxl_all$year,
                       lambda_Nplast)

# aa %>% dplyr::filter(year == 1) %>% ggplot() +  geom_tile(aes(x,y,fill = X1)) 

exc1 <- excursions.mc(aa[aa$year=="1", -c(1:3)],
                     alpha = 0.05,
                     u = 30, # (in our case it's natural scale)                     in log scale!!!
                     type = ">")  

exc1_Np <- excursions.mc(aa_Nplast[aa_Nplast$year=="1", - c(1:3)],
                        alpha = 0.05,
                        u = 30, # (in our case it's natural scale)                     in log scale!!!
                        type = ">")  

exc1P_joint = excursions.mc(p1[p1$year=="1",-c(1:3)],
                            alpha = 0.05,
                            u = 30, # (in our case it's natural scale)                     in log scale!!!
                            type = ">")  

exc1O_joint = excursions.mc(p2[p2$year=="1",-c(1:3)],
                            alpha = 0.05,
                            u = 10, # (in our case it's natural scale)                     in log scale!!!
                            type = ">")                    


exc2 = excursions.mc(aa[aa$year=="2",-c(1:3)], 
                     alpha = 0.05, 
                     u = 30, # in log scale!!! 
                     type = ">") 

exc2_Np = excursions.mc(aa_Nplast[aa_Nplast$year=="2",-c(1:3)],
                        alpha = 0.05,
                        u = 30, # (in our case it's natural scale)                     in log scale!!!
                        type = ">")  


exc2P_joint = excursions.mc(p1[p1$year=="2",-c(1:3)],
                            alpha = 0.05,
                            u = 30, # in log scale!!!
                            type = ">") 
exc2O_joint = excursions.mc(p2[p2$year=="2",-c(1:3)],
                            alpha = 0.05,
                            u = 10, # (in our case it's natural scale)                     in log scale!!!
                            type = ">")
exc3 =excursions.mc(aa[aa$year=="3",-c(1:3)], 
                    alpha = 0.05, 
                    u = 30, # in log scale!!! 
                    type = ">") 

exc3_Np = excursions.mc(aa_Nplast[aa_Nplast$year=="3",-c(1:3)],
                        alpha = 0.05,
                        u = 30, # (in our case it's natural scale)                     in log scale!!!
                        type = ">")  

exc3P_joint =excursions.mc(p1[p1$year=="3",-c(1:3)],
                           alpha = 0.05,
                           u = 30, # in log scale!!!
                           type = ">") 


exc3O_joint =excursions.mc(p2[p2$year=="3",-c(1:3)],
                           alpha = 0.05,
                           u = 10, # in log scale!!!
                           type = ">")                  

exc4 = excursions.mc(aa[aa$year=="4",-c(1:3)], 
                     alpha = 0.05, 
                     u = 30, # in log scale!!!
                     type = ">")  

exc4_Np = excursions.mc(aa_Nplast[aa_Nplast$year=="4",-c(1:3)],
                        alpha = 0.05,
                        u = 30, # (in our case it's natural scale)                     in log scale!!!
                        type = ">")  



exc4P_joint =excursions.mc(p1[p1$year=="4",-c(1:3)],
                           alpha = 0.05,
                           u = 30, # in log scale!!!
                           type = ">") 

exc4O_joint =excursions.mc(p2[p2$year=="4",-c(1:3)],
                           alpha = 0.05,
                           u = 10, # in log scale!!!
                           type = ">")  


exc5 = excursions.mc(aa[aa$year=="5",-c(1:3)], 
                     alpha = 0.05, 
                     u = 30, # in log scale!!! 
                     type = ">")  

exc5_Np = excursions.mc(aa_Nplast[aa_Nplast$year=="5",-c(1:3)],
                        alpha = 0.05,
                        u = 30, # (in our case it's natural scale)                     in log scale!!!
                        type = ">")  

exc5P_joint =excursions.mc(p1[p1$year=="5",-c(1:3)],
                           alpha = 0.05,
                           u = 30, # in log scale!!!
                           type = ">") 

exc5O_joint =excursions.mc(p2[p2$year=="5",-c(1:3)],
                           alpha = 0.05,
                           u = 10, # in log scale!!!
                           type = ">") 

exc6 = excursions.mc(aa[aa$year=="6",-c(1:3)], 
                     alpha = 0.05, 
                     u = 30, # in log scale!!! 
                     type = ">")  

exc6_Np = excursions.mc(aa_Nplast[aa_Nplast$year=="6",-c(1:3)],
                        alpha = 0.05,
                        u = 30, # (in our case it's natural scale)                     in log scale!!!
                        type = ">")  

exc6P_joint =excursions.mc(p1[p1$year=="6",-c(1:3)],
                           alpha = 0.05,
                           u = 30, # in log scale!!!
                           type = ">")

exc6O_joint =excursions.mc(p2[p2$year=="6",-c(1:3)],
                           alpha = 0.05,
                           u = 10, # in log scale!!!
                           type = ">")                  

exc7 = excursions.mc(aa[aa$year=="7",-c(1:3)], 
                     alpha = 0.05, 
                     u = 30, # in log scale!!! 
                     type = ">") 

exc7_Np = excursions.mc(aa_Nplast[aa_Nplast$year=="7",-c(1:3)],
                        alpha = 0.05,
                        u = 30, # (in our case it's natural scale)                     in log scale!!!
                        type = ">")  

exc7P_joint =excursions.mc(p1[p1$year=="7",-c(1:3)],
                           alpha = 0.05,
                           u = 30, # in log scale!!!
                           type = ">") 
exc7O_joint =excursions.mc(p2[p2$year=="7",-c(1:3)],
                           alpha = 0.05,
                           u = 10, # in log scale!!!
                           type = ">")                     

exc8 = excursions.mc(aa[aa$year=="8",-c(1:3)], 
                     alpha = 0.05, 
                     u = 30, # in log scale!!! 
                     type = ">") 

exc8_Np = excursions.mc(aa_Nplast[aa_Nplast$year=="8",-c(1:3)],
                        alpha = 0.05,
                        u = 30, # (in our case it's natural scale)                     in log scale!!!
                        type = ">")  

exc8P_joint =excursions.mc(p1[p1$year=="8",-c(1:3)],
                           alpha = 0.05,
                           u = 30, # in log scale!!!
                           type = ">")  
exc8O_joint =excursions.mc(p2[p2$year=="8",-c(1:3)],
                           alpha = 0.05,
                           u = 10, # in log scale!!!
                           type = ">")                   

myplot1PD <- data.frame(
  sf::st_coordinates(pxl_all[pxl_all$year=="1",]),z = exc1$F) %>% 
  ggplot2::ggplot() + ggplot2::geom_tile(ggplot2::aes(X, Y, fill = z))

myplot1OD <- data.frame(
  sf::st_coordinates(pxl_all[pxl_all$year=="1",]),z = exc1_Np$F) %>%
  ggplot2::ggplot() + ggplot2::geom_tile(aes(X, Y,fill = z))

 
myplot2PD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="2",]),z = exc2$F) %>%  ggplot2::ggplot() + geom_tile(aes(X,Y,fill = z))
myplot2OD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="2",]),z = exc2_Np$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))

myplot2O=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="2",]),z = exc2O_joint$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))

myplot3PD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="3",]),z = exc3$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))
myplot3OD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="3",]),z = exc3_Np$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))


myplot3O=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="3",]),z = exc3O_joint$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))


myplot4PD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="4",]),z = exc4$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))
myplot4OD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="4",]),z = exc4_Np$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))

myplot4O=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="4",]),z = exc4O_joint$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))

myplot5PD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="5",]),z = exc5$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))
myplot5OD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="5",]),z = exc5_Np$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))

myplot5O=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="5",]),z = exc5O_joint $F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))

myplot6PD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="6",]),z = exc6$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))
myplot6OD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="6",]),z = exc6_Np$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))


myplot6O=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="6",]),z = exc6O_joint$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))

myplot7PD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="7",]),z = exc7$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))
myplot7OD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="7",]),z = exc7_Np$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))


myplot7O=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="7",]),z = exc7O_joint$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))

myplot8PD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="8",]),z = exc8$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))
myplot8OD=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="8",]),z = exc8_Np$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))

myplot8O=data.frame(sf::st_coordinates(pxl_all[pxl_all$year=="8",]),z = exc8O_joint$F) %>%  ggplot2::ggplot() + ggplot2::geom_tile(aes(X,Y,fill = z))



year_plast1=myplot1PD + ggplot2::coord_equal() + 
  scico::scale_fill_scico(palette = "lajolla",direction=1)+   #theme_map+
  ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor =
          ggplot2::element_blank(),legend.position = c(0.91, 0.24),
        plot.title = ggplot2::element_text(size = 15, face = "bold"))+ 
  ggplot2::geom_sf(data = italy_sf, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggplot2::ggtitle("2013") 
#ggsave("exc2013.png", width = 13, height = 11, units = "cm")  
#year_other1=myplot1O + coord_equal() + 
#scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
#theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
#plot.title = element_text(size = 15, face = "bold"))+ 
#geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
#coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
#ggtitle("2013") 

year_Nplast1=myplot1OD + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2013") 
#ggsave("exc2013.png", width = 13, height = 11, units = "cm")  
year_other1=myplot1O + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2013") 

year_plast2=myplot2PD + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2014") 
# ggsave("exc2014.png", width = 13, height = 11, units = "cm") 

year_other2=myplot2O + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2014") 

year_plast2=myplot2PD + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2014") 
# ggsave("exc2014.png", width = 13, height = 11, units = "cm") 

year_other2=myplot2OD + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2014") 



year_plast3=myplot3P + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2015") 
# ggsave("exc2015.png", width = 13, height = 11, units = "cm") 

year_other3=myplot3O + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2015")


year_plast4=myplot4P + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2016") 
# ggsave("exc2016.png", width = 13, height = 11, units = "cm") 

year_other4=myplot4O + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2016") 

year_plast5=myplot5P + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2017") 
#ggsave("exc2017.png", width = 13, height = 11, units = "cm")

year_other5=myplot5O + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2017") 

year_plast6=myplot6P + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2018") 
# ggsave("exc2018.png", width = 13, height = 11, units = "cm") 

year_other6=myplot6O + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2018") 

year_plast7=myplot7P + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2019") 
# ggsave("exc2019.png", width = 13, height = 11, units = "cm") 

year_other7=myplot7O + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2019") 

year_plast8=myplot8P + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1)+   theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2020") 
# ggsave("exc2020.png", width = 13, height = 11, units = "cm") 

theme_map2 =  theme_light() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 

year_other8=myplot8O + coord_equal() + 
  scale_fill_scico(palette = "lajolla",direction=1) +theme_map+
  theme(panel.grid.major = element_blank(), panel.grid.minor =                 element_blank(),legend.position = c(0.91, 0.24),
        plot.title = element_text(size = 15, face = "bold"))+ 
  geom_sf(data = italy, alpha = 0.3, fill = 'white')+ 
  coord_sf(xlim=c(280, 815), ylim=c(4050, 4480))+
  ggtitle("2020")


nt <- theme(legend.position='none')
exc_plast=grid_arrange_shared_legend(year_plast1+nt,year_plast2+nt,year_plast3+nt,year_plast4+nt, year_plast5+nt,year_plast6+nt,year_plast7+nt,
                                     year_plast8+nt, ncol = 4, nrow = 2, position='bottom')
ggsave("exc_plast.png", exc_plast, width = 1280/72, height = 800/72, dpi = 72)

ggsave("exc_plast2.png", exc_plast, width = 9, height = 5, dpi=300)


exc_other=grid_arrange_shared_legend(year_other1+nt,year_other2+nt,year_other3+nt,year_other4+nt, year_other5+nt,year_other6+nt,year_other7+nt,
                                     year_other8+nt, ncol = 4, nrow = 2, position='bottom')
ggsave("exc_other.png", exc_other, width = 1280/72, height = 800/72, dpi = 72)







##' Other - to be relocated  ---------------------------------------------------
 

fit_marg0_ccd$summary.random$year_gamma_plast %>%
  dplyr::mutate(ID = as.character(c(2013:2021, 2023, 2024))) %>%
  ggplot2::ggplot() + 
  ggplot2::geom_errorbar(ggplot2::aes(ID, ymin = `0.025quant`, ymax = `0.975quant`),
                         alpha = 0.8, col="black",size=0.6) +
  ggplot2::theme_bw()+
  ggplot2::labs(y= "mean", x = "year") +
  ggplot2::geom_point(aes(ID,mean), col="red", size=2.2)+
  ggplot2::ylim(2.5,5.5)+
  ggplot2::ggtitle("Plastic") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))


fit_marg0_ccd$summary.random$year_gamma_Nplast %>%
  dplyr::mutate(ID = as.character(c(2013:2021, 2023, 2024))) %>%
  ggplot2::ggplot() + 
  ggplot2::geom_errorbar(ggplot2::aes(ID, ymin = `0.025quant`, ymax = `0.975quant`),
                         alpha = 0.8, col="black",size=0.6) +
  ggplot2::theme_bw()+
  ggplot2::labs(y= "mean", x = "year") +
  ggplot2::geom_point(aes(ID,mean), col="red", size=2.2)+
  ggplot2::ylim(2.5,5)+
  ggplot2::ggtitle("Other") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))











