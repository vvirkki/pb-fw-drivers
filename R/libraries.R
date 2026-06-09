# ------------------------------------------------------------------------------
#
# libraries
#
# Call libraries required for executing other scripts and functions therein.
# Tested configurations of libraries are listed in environment.txt.
#
# Analysis script for publication:
#
# "Regionally divergent drivers behind transgressions of the freshwater change
# planetary boundary"
#
# Vili Virkki, Lauren Seaby Andersen, Sofie te Wierik, Dieter Gerten,
# Miina Porkka
#
# Published in Nature Communications:
# https://doi.org/10.1038/s41467-026-73051-x
#
# Output data available in a Zenodo repository:
# https://doi.org/10.5281/zenodo.19663530
#
# Code availability:
# https://github.com/vvirkki/pb-fw-drivers
#
# Corresponding author & script author
# Vili Virkki (vili.virkki@uef.fi)
#
# Date:
# April 24th, 2026
#
# ------------------------------------------------------------------------------

library(qmap)
library(terra)
library(sf)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(ggplot2)
library(tmap)
library(viridisLite)
library(changepoint)
library(zoo)
library(doParallel)
library(gridExtra)
library(zyp)
library(exactextractr)
library(lubridate)
library(biscale)

options(dplyr.summarise.inform = FALSE)
