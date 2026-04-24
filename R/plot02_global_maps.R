# ------------------------------------------------------------------------------
#
# plot02_global_maps
#
# Create plots and outputs related to Fig. 3, Fig. 4, Fig. 5, Fig. S3, Fig. S4,
# Fig. S5, Fig. S6, and Fig. S7.
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


source("R/libraries.R")
source("R/plot00_functions.R")

# Setup ------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
variable <- args[1]
ardiv <- args[2]
indir <- args[3] # where data is read from
outdir <- args[4] # where outputs are saved
if (!dir.exists(outdir)) { dir.create(outdir) }

wyrs <- 30 # number of years to compare in the end of the time series
supp <- ifelse(grepl("supplementary", indir), TRUE, FALSE) # only for Fig. S3-S4

exportdir <- "zenodo_export"
if (!dir.exists(exportdir)) { dir.create(exportdir) }
exportvar <- case_when(
  variable == "dis" ~ "streamflow",
  variable == "rootmoist" ~ "soilmoisture"
)

# ice-dominant and not-ice-dominant catchments to correct anthromes ice mask
.hybas4_isIce <- c(2630, 3320)
.hybas4_notIce <- c(6540, 5540, 5520)

# Regional transgression status prep -------------------------------------------
# Use for
# Fig. 3: example graphic on decomposing CRF and DHF contributions, within loop

global_files <- list.files(paste0("output/", variable, "_", indir),
                           recursive = TRUE, full.names = TRUE)

files_in <- global_files[grepl("land_area_with_local_deviations_annual_mean_ensemble_median_IQR.csv", global_files) &
                         grepl(ardiv, global_files)]

geoms <- read_sf(paste0("data/ardiv/", ardiv, ".gpkg")) %>%
  mutate(area = paste0(ardiv, "_id", id))
ids <- unique(geoms$id)

coll_area_perc <- tibble()
coll_area_data_norm <- tibble()
coll_scenario_differences <- tibble()

for (i in 1:length(ids)) {

  area_files <- files_in[grepl(paste0("id", ids[i]), files_in)]

  if (i %% 10 == 0) { message(i) }
  if (length(area_files) == 0) { next }

  area_data <- area_files %>%
    lapply(read_delim, delim = ";", show_col_types = FALSE) %>%
    bind_rows() %>%
    distinct() %>%
    filter(class != "summed_deviations")

  check_scenario <- area_data %>%
    group_by(year, scenario) %>%
    summarise(nn = n()) %>%
    pull(nn) %>%
    unique()

  # each year-scenario-region should have two entries per normalised/non-normalised run
  if (check_scenario != 2 * length(unique(area_data$ts_normalised))) { stop(i) }

  maxyrs <- area_data %>%
    group_by(scenario) %>%
    summarise(maxyr = max(year))

  area_status <- area_data %>%
    left_join(maxyrs, by = "scenario")

  # for regionally normalised transgression plots (Fig. 4 and supplements)
  # normalise according to Eq. 1-3
  area_status_mean_norm <- area_status %>%
    filter(!ts_normalised) %>%
    filter(year > maxyr - wyrs) %>%
    mutate(ensemble_median = case_when(
      ensemble_median > boundary_upper_end ~ 1 + (ensemble_median - boundary_upper_end) / boundary_upper_end,
      ensemble_median >= boundary_baseline ~ (ensemble_median - boundary_baseline) / (boundary_upper_end - boundary_baseline),
      ensemble_median < boundary_baseline ~ (ensemble_median - boundary_baseline) / boundary_baseline
    )) %>%
    group_by(area, scenario, class) %>%
    summarise(status_end = mean(ensemble_median)) %>%
    pivot_wider(id_cols = c(area, class), names_from = scenario, values_from = status_end)

  coll_area_data_norm <- coll_area_data_norm %>%
    bind_rows(area_status_mean_norm)

  if (supp) { next }

  # for scenario contributions (Fig. 5, ...)
  area_status_perc <- area_status %>%
    filter(!ts_normalised) %>%
    filter(year > maxyr - wyrs) %>%
    select(area, year, scenario, class, ensemble_median, boundary_baseline, boundary_upper_end) %>%
    pivot_wider(id_cols = c(area, year, class, boundary_baseline, boundary_upper_end),
                names_from = scenario, values_from = ensemble_median)

  coll_area_perc <- coll_area_perc %>%
    bind_rows(area_status_perc)

  area_status_perc <- area_status_perc %>%
    select(-c(boundary_baseline, boundary_upper_end))

  scen_diff_drydev <- area_status_perc %>%
    filter(class == "dry_deviations") %>%
    .get_scenario_differences(col_diff_to = "counterclim_1901soc_default") %>%
    mutate(area = unique(area_status_perc$area)) %>%
    select(area, class, everything())

  scen_diff_wetdev <- area_status_perc %>%
    filter(class == "wet_deviations") %>%
    .get_scenario_differences(col_diff_to = "counterclim_1901soc_default") %>%
    mutate(area = unique(area_status_perc$area)) %>%
    select(area, class, everything())

  coll_scenario_differences <- coll_scenario_differences %>%
    bind_rows(scen_diff_drydev, scen_diff_wetdev)

  # Fig. 3 illustration
  if (i == 147) {

    arm2 <- area_status_mean_norm %>%
      pivot_longer(-c(area, class)) %>%
      filter(!grepl("baseline", name))

    nile <- ggplot(data = area_status %>% filter(ts_normalised)) +
      geom_line(aes(x = year, y = ensemble_median, color = scenario)) +
      geom_segment(data = arm2, aes(x = 1989, xend = 2019, y = value, yend = value, color = name),
                   linetype = "dashed") +
      geom_hline(yintercept = c(0, 1), color = "#e51aba") +
      facet_wrap(~ class) +
      scale_colour_manual(values = c("#575756", "#006085", "#f3a44c", "#bb5b0f")) +
      theme_minimal()

    ggsave(paste0(outdir, "/fig3_", variable,
                  "_method_example_nile.pdf"), nile)

  }

}

# Regional transgression status ------------------------------------------------
# Use for
# Fig. 4: regional transgression status, panels done one by one within loop
# Fig. S3, S4, S5, S6 different scenario/model setups of the same

status_long <- coll_area_data_norm %>%
  pivot_longer(-c(area, class), names_to = "scenario", values_to = "value") %>%
  mutate(id = as.numeric(str_replace(area, paste0(ardiv, "_id"), ""))) %>%
  select(id, everything())

all_files <- list.files("data/", recursive = TRUE, full.names = TRUE)
ice_mask_file <- all_files[grepl("ice_cells", all_files) & grepl("obsclim_histsoc_default", all_files)][1]

ice <- readRDS(ice_mask_file) %>%
  as.matrix() %>%
  rast(type = "xyz", crs = "+proj=lonlat")

total_area_rast <- rast(nrows = 360, ncols = 720, crs = crs(ice)) %>%
  cellSize(unit = "km")

total_area_from_rast <- exact_extract(total_area_rast, geoms, fun = "sum",
                                      force_df = TRUE, append_cols = "id") %>%
  rename(total_area = sum) %>%
  group_by(id) %>%
  summarise(total_area = sum(total_area)) # for datetime line overlapping polygons

ice_area <- exact_extract(ice * crop(total_area_rast, ice), geoms, fun = "sum",
                          force_df = TRUE, append_cols = "id") %>%
  rename(ice_area = sum) %>%
  group_by(id) %>%
  summarise(ice_area = sum(ice_area))  # for datetime line overlapping polygons

ice_fraction <- total_area_from_rast %>%
  left_join(ice_area, by = "id") %>%
  as_tibble() %>%
  mutate(ice_perc = ice_area / total_area)

ice_catchments <- ice_fraction %>%
  filter((ice_perc > 0.5 | id %in% .hybas4_isIce) & !id %in% .hybas4_notIce) %>%
  pull(id)

if (!supp) {
  saveRDS(ice_catchments, paste0(outdir, "/fig4_maps_", ardiv, "_ice_catchments.rds"))
}

nodata_catchments <- geoms %>%
  filter(!area %in% unique(status_long$area) & !id %in% ice_catchments) %>%
  pull(id)

plt_values <- expand_grid(id = geoms$id,
                         class = unique(status_long$class),
                         scenario = unique(status_long$scenario)) %>%
  left_join(status_long %>% select(-area), by = c("id", "class", "scenario")) %>%
  filter(!grepl("baseline_sd", scenario)) %>%
  mutate(value = case_when(
    id %in% ice_catchments ~ -999, # ice catchment
    is.na(value) ~ NA,             # true missing value (no data)
    TRUE ~ value
  )) %>%
  mutate(facet_by = paste0(class, "_", scenario))

brk <- c(-Inf, -1, -0.75, -0.25, 0, 0.25, 0.75, 1, 1.25, 2, Inf)
pal <- c("#b5b5b5", "#80cdc1", "#b1e1da", "#d4ede9", "#f5ebd2", "#efdcad",
         "#dec17c", "#feaa38", "#ef7818", "#cc4c02")
lbl <- paste0(brk[-length(brk)], "...", brk[-1])
lbl[1] <- "ice"
lbl[2] <- str_replace(lbl[2], "0-", "-")
lbl[length(lbl)] <- str_replace(lbl[length(lbl)], "-Inf", "-")

facets <- unique(plt_values$facet_by)

for (i in 1:length(facets)) {

  facet_values <- plt_values %>%
    filter(facet_by == facets[i]) %>%
    select(id, value) %>%
    distinct()

  map <- tm_shape(geoms %>% left_join(facet_values, by = "id")) +
    tm_polygons(fill = "value",
                fill.scale = tm_scale_intervals(
                  breaks = brk,
                  values = pal,
                  labels = lbl,
                  midpoint = NA,
                  value.na = "grey50"),
                lwd = 0) +
    tm_crs("+proj=robin") +
    tm_layout(frame = FALSE,
              legend.show = FALSE)

  supp_spec <- ifelse(supp, "for_figS3S4_", "")

  prefix <- case_when(
    grepl("obsclim_histsoc_default", facets[i]) ~ "fig4",
    grepl("obsclim_1901soc_default", facets[i]) ~ "figS5",
    grepl("counterclim_histsoc_default", facets[i]) ~ "figS6",
    grepl("counterclim_1901soc_default", facets[i]) ~ "not_used"
  )

  if (supp & prefix %in% c("fig4", "not_used")) { next }

  tmap_save(map, paste0(outdir, "/", prefix, "_", supp_spec, variable, "_",
                        ardiv, "_", facets[i], ".pdf"),
            width = 15, height = 10, units = "cm")

  tmap_save(map, paste0(outdir, "/", prefix, "_", supp_spec, variable, "_",
                        ardiv, "_", facets[i], ".tiff"),
            width = 8, height = 4, units = "cm", dpi = 600)

}

plt_zenodo_export <- plt_values %>%
  distinct() %>%
  mutate(value = ifelse(value == -999, NA, value)) %>% # ice catchments to NA
  select(-facet_by) %>%
  pivot_wider(id_cols = c("id", "class"), names_from = "scenario", values_from = "value") %>%
  select(-counterclim_1901soc_default) %>%
  mutate(variable = exportvar, .before = 1) %>%
  mutate(areal_division = ardiv, .after = 1) %>%
  rename(scenario_dhf_only = counterclim_histsoc_default,
         scenario_crf_only = obsclim_1901soc_default,
         scenario_historical = obsclim_histsoc_default,
         area_id = id) %>%
  arrange(area_id, class)

if (supp) {

  # # # zenodo export # # #
  # needs to be fetched manually from figure folders
  write.table(plt_zenodo_export, paste0(outdir, "/figS3S4_", exportvar, "_plotted_values.csv"),
              sep = ";", row.names = FALSE)

  quit(save = "no")

}

annual_values_zenodo_export <- coll_area_perc %>%
  mutate(id = as.numeric(str_replace(area, paste0(ardiv, "_id"), ""))) %>%
  mutate(variable = exportvar, .before = 1) %>%
  mutate(areal_division = ardiv, .after = 1) %>%
  rename(scenario_baseline = counterclim_1901soc_default,
         scenario_dhf_only = counterclim_histsoc_default,
         scenario_crf_only = obsclim_1901soc_default,
         scenario_historical = obsclim_histsoc_default,
         boundary_median = boundary_baseline,
         area_id = id) %>%
  select(-area) %>%
  relocate(area_id, .after = 2) %>%
  relocate(starts_with("boundary"), .after = -1) %>%
  arrange(area_id, class, year)

# # # zenodo export # # #
write.table(plt_zenodo_export, paste0(exportdir, "/fig4_figS5S6_", exportvar, "_plotted_values.csv"),
            sep = ";", row.names = FALSE)

# # # zenodo export # # #
write.table(annual_values_zenodo_export, paste0(exportdir, "/fig4_fig5_figS5S6S7_", exportvar, "_annual_values.csv"),
            sep = ";", row.names = FALSE)

# Regional driver contributions to transgression status ------------------------
# Use for
# Fig. 5: scenario contribution classification (strength)
# Fig. S7: scenario contribution classification (direction)
# intermediate exports to plot03_ that draws Fig. 6

contr_relative <- coll_scenario_differences %>%
  mutate(rel_diff_means = case_when(
    sd_to == 0 ~ 0,
    TRUE ~ diff_means / sd_to)
  ) %>%
  select(area, class, scenario_from, rel_diff_means) %>%
  pivot_wider(id_cols = c(area, class),
              names_from = scenario_from, values_from = rel_diff_means) %>%
  rename(contr_dhf = counterclim_histsoc_default,
         contr_crf = obsclim_1901soc_default,
         contr_both = obsclim_histsoc_default)

# for Fig. 6 to be composed in plot03_
percranks <- contr_relative %>%
  group_by(class) %>%
  mutate(contr_dhf_percrank = percent_rank(contr_dhf),
         contr_crf_percrank = percent_rank(contr_crf),
         contr_dhf_crf_percrank = percent_rank(contr_both)) %>%
  ungroup() %>%
  mutate(id = as.numeric(str_replace(area, paste0(ardiv, "_id"), "")),
         isIce = id %in% ice_catchments,
         isNodata = id %in% nodata_catchments,
         across(ends_with("percrank"), ~ ifelse(isIce | isNodata, NA, .x))) %>%
  select(id, class, ends_with("percrank"), isIce, isNodata)

saveRDS(percranks, paste0(outdir, "/fig6_", variable, "_", ardiv,
                          "_percranks.rds"))

# for Fig. 5
scen_cntr_strength <- coll_scenario_differences %>%
  mutate(rel_diff_means = diff_means / sd_to) %>%
  select(area, class, scenario_from, rel_diff_means, wilcox_pval) %>%
  mutate(bi_pal_class = case_when(
    wilcox_pval > 0.05 | rel_diff_means <= 0 ~ 1,
    rel_diff_means < 1 ~ 2,
    rel_diff_means < 2 ~ 3,
    TRUE ~ 4
  )) %>%
  select(-c(rel_diff_means, wilcox_pval)) %>%
  pivot_wider(id_cols = c(area, class),
              names_from = scenario_from, values_from = bi_pal_class) %>%
  select(-obsclim_histsoc_default) %>% # CRF and DHF scenarios only
  mutate(impact_class = paste0(obsclim_1901soc_default, "-", counterclim_histsoc_default),
         id = as.numeric(str_replace(area, paste0(ardiv, "_id"), ""))) %>%
  select(id, class, impact_class)

export_scendiff <- coll_scenario_differences %>%
  mutate(rel_diff_means = diff_means / sd_to) %>%
  rowwise() %>%
  mutate(rel_diff_means = ifelse(wilcox_pval > 0.05, NA, rel_diff_means)) %>%
  ungroup() %>%
  select(area, class, scenario_from, rel_diff_means) %>%
  pivot_wider(id_cols = c("area", "class"), names_from = "scenario_from", values_from = "rel_diff_means") %>%
  mutate(area = as.numeric(str_replace(area, paste0(ardiv, "_id"), ""))) %>%
  rename(id = area)

# for Fig. S7
cntr_directions <- coll_scenario_differences %>%
  select(area, class, scenario_from, diff_means, wilcox_pval) %>%
  mutate(diff_means = ifelse(wilcox_pval > 0.05 | diff_means == 0, 0, diff_means)) %>%
  select(-wilcox_pval) %>%
  pivot_wider(id_cols = c(area, class),
              names_from = scenario_from, values_from = diff_means) %>%
  rename(crf = obsclim_1901soc_default,
         dhf = counterclim_histsoc_default,
         net = obsclim_histsoc_default) %>%
  mutate(countering = case_when(
    dhf < 0 & crf > 0 ~ "dhf_down_crf_up",
    dhf > 0 & crf < 0 ~ "dhf_up_crf_down",
    dhf < 0 & crf < 0 ~ "dhf_crf_down",
    dhf > 0 & crf > 0 ~ "dhf_crf_up",
    dhf == 0 | crf == 0 ~ "insign_cntr", # one or both insignificant
  )) %>%
  mutate(net_change = case_when(
    net > 0 ~ "net_increase",
    net < 0 ~ "net_decrease",
    net == 0 ~ "net_insign"
  )) %>%
  mutate(direction_class = case_when(
    countering == "dhf_crf_down" ~ "dhf_crf_down",
    countering == "insign_cntr" ~ "insign_cntr",
    net_change == "net_decrease" ~ "net_decrease",
    net_change == "net_insign" ~ "net_insign",
    TRUE ~ paste0(countering, "_", net_change)
  )) %>%
  mutate(id = as.numeric(str_replace(area, paste0(ardiv, "_id"), ""))) %>%
  select(id, class, direction_class)


# Fig. 5, Fig. S7 prep and output ----------------------------------------------

plt_values <- expand_grid(id = geoms$id,
                          class = unique(scen_cntr_strength$class)) %>%
  left_join(scen_cntr_strength, by = c("id", "class")) %>%
  left_join(cntr_directions, by = c("id", "class")) %>%
  mutate(across(c(impact_class, direction_class), ~ case_when(
    id %in% ice_catchments ~ "xx_ice",  # ice catchment
    is.na(.x) ~ NA,                     # no data
    TRUE ~ .x
  )))


bi_pal_nm <- "DkBlue2"
bi_pal(bi_pal_nm, 4) %>% ggsave(filename = paste0(outdir, "/fig5_legend.pdf")) # save legend

impact_classes <- names(bi_pal(bi_pal_nm, 4, preview = FALSE))
pal <- c(bi_pal(bi_pal_nm, 4, preview = FALSE), "#b5b5b5") %>%
  setNames(c(impact_classes, "xx_ice"))

classes <- unique(plt_values$class)
for (i in 1:length(classes)) {

  impact_class_values <- plt_values %>%
    filter(class == classes[i]) %>%
    select(id, impact_class) %>%
    distinct()

  map <- tm_shape(geoms %>% left_join(impact_class_values, by = "id")) +
    tm_polygons(fill = "impact_class",
                fill.scale = tm_scale_categorical(
                  values = pal[unique(impact_class_values$impact_class)],
                  value.na = "grey50"
                ),
                lwd = 0) +
    tm_crs("+proj=robin") +
    tm_layout(frame = FALSE,
              legend.show = TRUE)

  tmap_save(map, paste0(outdir, "/fig5_", variable, "_", ardiv,
                        "_", classes[i], "_scenario_contributions.pdf"),
            width = 30, height = 30, units = "cm")

  tmap_save(map + tm_layout(legend.show = FALSE),
            paste0(outdir, "/fig5_", variable, "_", ardiv,
                   "_", classes[i], "_scenario_contributions.tiff"),
            width = 8, height = 4, units = "cm", dpi = 600)

  direction_class_values <- plt_values %>%
    filter(class == classes[i]) %>%
    select(id, direction_class) %>%
    distinct()

  map <- tm_shape(geoms %>% left_join(direction_class_values, by = "id")) +
    tm_polygons(fill = "direction_class",
                fill.scale = tm_scale_categorical(
                  values = c("#daf5e1", "#b5ab6b", "#7fa074", "#90719f",
                             "#daf5e1", "#daf5e1", "#daf5e1", "#daf5e1", "#b5b5b5"),
                  value.na = "#b5b5b5"
                ),
                lwd = 0) +
    tm_crs("+proj=robin") +
    tm_layout(frame = FALSE,
              legend.show = TRUE)

  tmap_save(map, paste0(outdir, "/figS7_", variable, "_", ardiv,
                        "_", classes[i], "_scenario_directions.pdf"),
            width = 30, height = 30, units = "cm")

  tmap_save(map + tm_layout(legend.show = FALSE),
            paste0(outdir, "/figS7_", variable, "_", ardiv,
                   "_", classes[i], "_scenario_directions.tiff"),
            width = 8, height = 4, units = "cm", dpi = 600)

}

# # # zenodo export # # #
plt_values %>%
  left_join(export_scendiff, by = c("id", "class")) %>%
  select(-impact_class) %>%
  relocate(direction_class, .after = -1) %>%
  mutate(variable = exportvar, .before = 1) %>%
  mutate(areal_division = ardiv, .after = 1,
         direction_class = ifelse(direction_class == "xx_ice", NA, direction_class)) %>%
  rename(cntr_dhf_relative_to_bsl_sd = counterclim_histsoc_default,
         cntr_crf_relative_to_bsl_sd = obsclim_1901soc_default,
         cntr_historical_relative_to_bsl_sd = obsclim_histsoc_default,
         cntr_direction = direction_class,
         area_id = id) %>%
  mutate(across(starts_with("cntr_"), ~ ifelse(is.na(cntr_direction), NA, .x))) %>%
  arrange(area_id, class) %>%
  write.table(paste0(exportdir, "/fig5_figS7_", exportvar, "_plotted_values.csv"),
              sep = ";", row.names = FALSE)
