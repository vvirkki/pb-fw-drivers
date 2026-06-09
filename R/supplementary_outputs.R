# ------------------------------------------------------------------------------
#
# supplementary_outputs
#
# Supplementary intermediate output script that creates outputs that are further
# post-processed into figures used in the manuscript.
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

# Load scripts -----------------------------------------------------------------

source("R/libraries.R")

# data export & plotting scripts
source("R/00_configure_run.R")
source("R/output01_land_area_with_local_deviations.R")

# save sessionInfo()
sink(paste0("sessionInfo ", Sys.time(), ".txt"))
sessionInfo()
sf::sf_extSoftVersion()
sink()

.find_ensmems <- function(cfg) {

  dirs <- list.dirs(paste0("data/", cfg$variable, "/", cfg$period_label[2],
                           "/", cfg$variable, "_monthly"))

  imdirs <- dirs %>%
    lapply(str_split, "/") %>%
    lapply(unlist) %>%
    lapply(tail, 1) %>%
    unlist()

  imdirs <- dirs[which(imdirs %in% (str_replace_all(cfg$impactmodel, "-", "_")))]

  emdirs <- c()
  for (i in 1:length(imdirs)) {
    emdirs <- c(emdirs, list.dirs(imdirs[i])[-1])
  }

  ensmems_lst <- emdirs %>%
    lapply(str_split, "/") %>%
    lapply(unlist) %>%
    lapply(tail, 2) %>%
    lapply(str_replace_all, "_", "-")

  return (ensmems_lst)

}

# Create outputs ---------------------------------------------------------------

if (!dir.exists("output")) { dir.create("output") }

# ISIMIP 3a --------------------------------------------------------------------

run_variable <- commandArgs(trailingOnly = TRUE)

# template config
cfg <- configure_run(paste0("configs/", run_variable, " main isimip3a"))

scenarios <- c("obsclim_histsoc_default",
               "obsclim_1901soc_default",
               "counterclim_histsoc_default")

if (run_variable == "dis") {

  cfg$impactmodel <- c("h08", "hydropy", "jules-w2", "lpjml5-7-10-fire",
                       "miroc-integ-land", "watergap2-2e")
  cfg$forcing <- c("gswp3-w5e5", "20crv3-era5", "20crv3-w5e5")

} else if (run_variable == "rootmoist") {

  cfg$impactmodel <- c("hydropy", "lpjml5-7-10-fire",
                       "miroc-integ-land", "web-dhm-sg")
  cfg$forcing <- c("gswp3-w5e5", "20crv3-era5", "20crv3-w5e5")

}

# outputs global ---------------------------------------------------------------

global_output_out <- paste0("output/", cfg$variable, "_supplementary_ensemble")
ensmem_output_out <- paste0("output/", cfg$variable, "_supplementary_single_ensemble_members")

if (!dir.exists(global_output_out)) { dir.create(global_output_out) }
if (!dir.exists(ensmem_output_out)) { dir.create(ensmem_output_out) }

storecfg <- cfg # for individual ensemble member looping
for (s in 1:length(scenarios)) {

  scenario_output_global <- paste0(global_output_out, "/", scenarios[s], "/")
  scenario_output_ensmem <- paste0(ensmem_output_out, "/", scenarios[s], "/")

  if (!dir.exists(scenario_output_global)) { dir.create(scenario_output_global) }
  if (!dir.exists(scenario_output_ensmem)) { dir.create(scenario_output_ensmem) }

  cfg$period_label <- c("counterclim_1901soc_default", scenarios[s])
  cfg$ensemble_selection <- "ensemble_all.csv"

  # area aggregations globally for Fig. S12
  cfg$areal_division <- ""
  cfg$draw_ts_var <- "landShareInClass"
  cfg$normalise_ts <- FALSE
  fldr_created <- output_land_area_with_local_deviations(cfg)

  file.rename(fldr_created,
              str_replace(fldr_created, "output/", scenario_output_global))

  # outputs ensemble members ---------------------------------------------------
  message(scenarios[s])
  ensmems <- .find_ensmems(cfg)
  for (em in 1:length(ensmems)) {

    em_out <- paste0(ensmem_output_out, "/", scenarios[s], "/",
                     paste(ensmems[[em]], collapse = "_"), "/")

    # run regional outputs for certain ensemble members for Fig. S3-S4
    if ((run_variable == "dis" & ensmems[[em]][1] == "watergap2-2e") |
        (run_variable == "rootmoist" & ensmems[[em]][1] == "miroc-integ-land")) {

      if (!dir.exists(em_out)) { dir.create(em_out) }

      message(paste(ensmems[[em]][1], ensmems[[em]][2], collapse = " "))
      cfg$impactmodel <- ensmems[[em]][1]
      cfg$forcing <- ensmems[[em]][2]
      cfg$areal_division <- "hybas4"
      cfg$draw_ts_var <- "landShareInClass"
      cfg$normalise_ts <- FALSE
      cfg$ensemble_selection <- ""

      fldr_created <- output_land_area_with_local_deviations(cfg)
      file.rename(fldr_created,
                  str_replace(fldr_created, "output/", em_out))

      cfg <- storecfg
      cfg$period_label <- c("counterclim_1901soc_default", scenarios[s])

    } else {
      next
    }

  }

}
