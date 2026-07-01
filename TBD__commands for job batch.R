


##'``
##'`R` 

##' Commands to run remotely ---------------------------------------------------
##' 
##' First line ----------------------------------------------------------------#
install.packages(
  c(
    "inlabru",
    "terra",
    "tidyverse",
    "lubridate",
    "sf",
    "raster",
    "scico",
    "patchwork",
    "sp",
    "geodata",
    "spatstat.geom",
    "rvest",
    "xml2",
    "knitr",
    "magrittr",
    "rlang"
  )
)
##' Second line ---------------------------------------------------------------#
install.packages(
  "INLA",
  repos=c(
    getOption("repos"),
    INLA="https://inla.r-inla-download.org/R/testing"
  ),
  dep=TRUE
)