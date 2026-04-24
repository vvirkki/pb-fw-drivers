# ------------------------------------------------------------------------------
#
# temp02_infilling_diagnostics
#
# Query statistics created in runtime logging for missing values that were
# infilled in data preparation.
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
# ADD LINK
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

library(tidyverse)
options(dplyr.summarise.inform = FALSE)

.scenarios <- c("counterclim_1901soc_default",
                "obsclim_histsoc_default",
                "obsclim_1901soc_default",
                "counterclim_histsoc_default",
                "picontrol_1850soc_default",
                "historical_histsoc_default",
                "ssp126_2015soc_default",
                "ssp370_2015soc_default",
                "ssp585_2015soc_default")

.parse_diagnostics <- function(diag_rds_in, variable) {

  dg <- readRDS(diag_rds_in) %>%
    as_tibble()

  ncells_affected <- nrow(dg)

  dg_mat <- dg %>%
    select(starts_with(!!variable)) %>%
    as.matrix()

  nvals_negative <- length(which(dg_mat < 0))
  nvals_na <- length(which(is.na(dg_mat)))

  ret <- dg %>%
    group_by(impactmodel, forcing, month) %>%
    summarise() %>%
    mutate(scenario = .scenarios[which(sapply(.scenarios, grepl, x = diag_rds_in, fixed = TRUE))],
           ncells_affected = ncells_affected,
           nvals_negative = nvals_negative,
           nvals_na = nvals_na) %>%
    ungroup()

  return (ret)

}

logs_files <- list.files("logs", full.names = TRUE)

diag_dis <- logs_files[grepl("dis_infilled", logs_files)] %>%
  lapply(.parse_diagnostics, variable = "dis") %>%
  bind_rows()

diag_rootmoist <- logs_files[grepl("rootmoist_infilled", logs_files)] %>%
  lapply(.parse_diagnostics, variable = "rootmoist") %>%
  bind_rows()

diag_dis %>%
  group_by(impactmodel, forcing, scenario) %>%
  summarise(min_cells = min(ncells_affected),
            max_cells = max(ncells_affected),
            min_negative = min(nvals_negative),
            max_negative = max(nvals_negative),
            min_na = min(nvals_na),
            max_na = max(nvals_na))

diag_rootmoist %>%
  group_by(impactmodel, forcing, scenario) %>%
  summarise(min_cells = min(ncells_affected),
            max_cells = max(ncells_affected),
            min_negative = min(nvals_negative),
            max_negative = max(nvals_negative),
            min_na = min(nvals_na),
            max_na = max(nvals_na)) %>%
  arrange(-max_cells)
