# ------------------------------------------------------------------------------
#
# temp01_planetary_health_check_export
#
# Pick and export specific output data items to be transferred for composing
# Planetary Health check reports (https://www.planetaryhealthcheck.org).
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

.get_latest_mean <- function (file_in, tbl_in, bound) {

  ret <- tbl_in %>%
    filter(class == bound) %>%
    select(area, year, class, ensemble_median, boundary_upper_end) %>%
    group_by(area, class, boundary_upper_end) %>%
    summarise(mean_latest_10yr = mean(ensemble_median)) %>%
    mutate(area = as.numeric(str_replace(area, "hybas3_id", "")),
           variable = str_replace_all(variables[str_detect(file_in, variables)], "/", "")) %>%
    rename(hybas3_id = area)

  return (ret)

}

outdir <- "zenodo_export/planetary_health_check"
if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }

# ------------------------------------------------------------------------------
# Global outputs for PHC2025 Figure 31 & 33 / PHC2026 Figure TBD
# dis = blue water; rootmoist = green water
# ------------------------------------------------------------------------------

output_files <- list.files("output", recursive = TRUE, full.names = TRUE)
variables <- c("/dis", "/rootmoist")
bounds <- c("summed_deviations")

global_deviations <- output_files[grepl("global_land_area_with_local_deviations_annual_mean_ensemble_median_IQR.csv", output_files) &
                                  grepl("obsclim_histsoc_default", output_files) &
                                  grepl("main_ensemble", output_files)]

for (i in 1:length(global_deviations)) {

  read_delim(global_deviations[i], delim = ";", show_col_types = FALSE) %>%
    filter(class == "summed_deviations" & scenario == "obsclim_histsoc_default") %>%
    select(-ts_normalised) %>%
    mutate(variable = str_replace_all(variables[str_detect(global_deviations[i], variables)], "/", "")) %>%
    select(variable, everything()) %>%
    write.table(paste0(outdir, "/", tail(unlist(str_split(global_deviations[i], "/")), 1)),
                sep = ";", row.names = FALSE)

}

# ------------------------------------------------------------------------------
# Regional outputs for PHC2025 Figure 30 & 32 / PHC2026 Figure TBD
# dis = blue water; rootmoist = green water
# ------------------------------------------------------------------------------

hybas3_deviations <- output_files[grepl("land_area_with_local_deviations_annual_mean_ensemble_median_IQR.csv", output_files) &
                                  grepl("obsclim_histsoc_default", output_files) &
                                  grepl("main_ensemble", output_files) &
                                  grepl("hybas3", output_files)]

tbl_all <- tibble()

for (i in 1:length(hybas3_deviations)) {

  tbl_in <- read_delim(hybas3_deviations[i], delim = ";", show_col_types = FALSE) %>%
    filter(year > max(.$year - 10) & scenario == "obsclim_histsoc_default")

  for (j in 1:length(bounds)) {
    tbl_all <- bind_rows(tbl_all, .get_latest_mean(hybas3_deviations[i], tbl_in, bounds[j]))
  }

}

variables_out <- c("dis", "rootmoist")

for (j in 1:length(variables_out)) {
  for (k in 1:length(bounds)) {

    tbl_out <- tbl_all %>%
      filter(variable == variables_out[j] &
             class == bounds[k]) %>%
      select(hybas3_id, variable, class, boundary_upper_end, mean_latest_10yr)

    file_out <- paste0(outdir, "/", variables_out[j], "_", bounds[k], "_hybas3_latest_10yr_means.csv")
    write.table(tbl_out, file_out, sep = ";", row.names = FALSE)

  }
}
