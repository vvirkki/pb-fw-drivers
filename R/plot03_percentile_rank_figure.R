# ------------------------------------------------------------------------------
#
# plot03_percentile_rank_figure
#
# Create plots and outputs related to Fig. 6 and Fig. S8.
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

# percentile rank combination map ----------------------------------------------
# Use for
# Fig. 6: synthesising analysis of scenario percentile ranks and auxiliary data
# Fig. S8: maps of auxiliary variables

args <- commandArgs(trailingOnly = TRUE)
ardiv <- args[1]
indir <- args[2] # where previous outputs are read from
outdir <- args[3] # where outputs are saved
if (!dir.exists(outdir)) { dir.create(outdir) }

exportdir <- "zenodo_export"
if (!dir.exists(exportdir)) { dir.create(exportdir) }

geoms <- read_sf(paste0("data/ardiv/", ardiv, ".gpkg")) %>%
  mutate(area = paste0(ardiv, "_id", id))
ids <- unique(geoms$id)

dis_percranks <- readRDS(paste0(indir, "/fig6_dis_",
                                ardiv, "_percranks.rds")) %>%
  pivot_wider(id_cols = c(id, isIce, isNodata), names_from = class,
              values_from = ends_with("percrank"), names_prefix = paste0("dis_")) %>%
  rename_with(~ paste0("dis_", .x), c(isIce, isNodata))

rootmoist_percranks <- readRDS(paste0(indir, "/fig6_rootmoist_",
                                      ardiv, "_percranks.rds")) %>%
  pivot_wider(id_cols = c(id, isIce, isNodata), names_from = class,
              values_from = ends_with("percrank"), names_prefix = paste0("rootmoist_")) %>%
  rename_with(~ paste0("rootmoist_", .x), c(isIce, isNodata))

ice_catchments <- readRDS(paste0(outdir, "/fig4_maps_", ardiv, "_ice_catchments.rds"))

percranks_together <- expand_grid(id = unique(geoms$id)) %>%
  left_join(dis_percranks, by = "id") %>%
  left_join(rootmoist_percranks, by = "id")

plt_values <- percranks_together %>%
  rowwise() %>%
  mutate(median_percrank_dhf = median(c_across(matches("contr_dhf_percrank"))),
         median_percrank_crf = median(c_across(matches("contr_crf_percrank"))),
         median_percrank_dhf_crf = median(c_across(matches("contr_dhf_crf_percrank"))),
         isIce = any(c_across(ends_with("isIce"))),
         isNodata = any(c_across(ends_with("isNodata")))) %>%
  ungroup() %>%
  select(id, isIce, isNodata, starts_with("median_percrank")) %>%
  mutate(across(starts_with("median_percrank"), ~ case_when(
    isIce | id %in% ice_catchments ~ -999,
    isNodata & !isIce ~ NA,
    TRUE ~ .x
  ))) %>%
  pivot_longer(-c(id, isIce, isNodata), names_to = "facet")


zenodo_export_percranks <- plt_values %>%
  mutate(value = ifelse(isIce | id %in% ice_catchments, NA, value)) %>%
  filter(facet != "median_percrank_dhf_crf") %>%
  select(id, facet, value) %>%
  pivot_wider(id_cols = "id", names_from = "facet", values_from = "value") %>%
  left_join(percranks_together %>%
              select(id, starts_with("contr_dhf"), starts_with("contr_crf")) %>%
              select(-starts_with("contr_dhf_crf")),
            by = "id") %>%
  rename_with(~ str_replace(.x, "contr", "cntr")) %>%
  rename_with(~ str_replace(.x, "dis", "streamflow")) %>%
  rename_with(~ str_replace(.x, "rootmoist", "soilmoisture")) %>%
  arrange(id) %>%
  mutate(areal_division = ardiv, .before = 1) %>%
  rename(area_id = id)

# # # zenodo export # # #
write.table(zenodo_export_percranks, paste0(exportdir, "/fig6ab_percranks.csv"),
            sep = ";", row.names = FALSE)

# create and save maps ---------------------------------------------------------

percrank_variability <- plt_values %>%
  filter(value > 0 & value < 1) %>%
  group_by(facet) %>%
  summarise(global_mean = mean(value),
            global_sd = sd(value))

write.table(percrank_variability,
            paste0(outdir, "/fig6_percrank_variability.csv"),
            sep = ";", row.names = FALSE)

brk <- c(-Inf, seq(0, 1, 0.2))
facets <- unique(plt_values$facet)
for (i in 1:length(facets)) {

  plt_geoms <- geoms %>%
    left_join(plt_values %>% filter(facet == facets[i]), by = "id")

  mm <- tm_shape(plt_geoms) +
    tm_polygons(fill = "value",
                fill.scale = tm_scale_intervals(
                  breaks = brk,
                  values = c("#e4e3de", "#fdfad5", "#fed98f", "#f8992e", "#d76127", "#993720"),
                  value.na = "grey50"),
                lwd = 0) +
    tm_crs("+proj=robin") +
    tm_layout(frame = FALSE)

  tmap_save(mm, paste0(outdir, "/fig6_", ardiv, "_", facets[i], ".pdf"))

  tmap_save(mm + tm_layout(legend.show = FALSE),
            paste0(outdir, "/fig6_", ardiv, "_", facets[i], ".tiff"),
            width = 8, height = 4, units = "cm", dpi = 600)

}

# percentile rank auxiliary data bins ------------------------------------------

# REQUIRED EXTERNAL DATA:
# HYDE 3.5 population data
# https://geo.public.data.uu.nl/vault-hyde/hyde35_c9_apr2025%5B1749214444%5D/original/gbc2025_7apr_base/NetCDF/
hyde_population <- rast("Data/auxdata/hyde35/population.nc") %>%
  setNames(paste0("pop", year(time(.))))

pop_in_hybas4 <- .ee_sum(hyde_population$pop2019, geoms) %>%
  group_by(id) %>%
  summarise(pop2019 = sum(pop2019)) %>%
  as_tibble()

# REQUIRED EXTERNAL DATA:
# Terrestrial MSA for 2015 by Schipper et al. (2020) (https://doi.org/10.1111/gcb.14848)
# https://www.globio.info/globio-data-downloads
msa_2015 <- rast("Data/auxdata/Globio4_TerrestrialMSA_10sec_2015/TerrestrialMSA_2015_World.tif")
msa_in_hybas4 <- .ee_area_weighted_mean(msa_2015, geoms) %>%
  as_tibble()

# REQUIRED EXTERNAL DATA:
# HANPP for 2010 by Kastner et al. (2021) (https://doi.org/10.1111/gcb.15932)
# https://zenodo.org/records/7313791
hanpp_layers <- rast("Data/auxdata/Kastner_et_al_2021_HANPP/2010_results_mean_144_runs.tif")
nm_hanpp <- names(hanpp_layers)

landmask <- sum(hanpp_layers[[nm_hanpp[grepl("area_km2", nm_hanpp)]]])
landmask[landmask == 0] <- NA

npp_pot <- sum(hanpp_layers[[nm_hanpp[grepl("NPPpot_", nm_hanpp)]]])
hanpp_luc <- sum(hanpp_layers[[nm_hanpp[grepl("HANPPluc_", nm_hanpp)]]])
hanpp_harv <- sum(hanpp_layers[[nm_hanpp[grepl("HANPPharv_", nm_hanpp)]]])
hanpp_def <- sum(hanpp_layers[[nm_hanpp[grepl("HANPPdef_", nm_hanpp)]]])

hanpp_perc <- (hanpp_luc + hanpp_harv + hanpp_def) / npp_pot
hanpp_perc[npp_pot == 0] <- 0
hanpp_perc[is.na(landmask)] <- NA

names(hanpp_perc) <- "hanpp_perc"

hanpp_in_hybas4 <- .ee_area_weighted_mean(hanpp_perc, geoms) %>%
  as_tibble()

# auxiliary variable boxplots --------------------------------------------------

percrank_bins <- plt_values %>%
  mutate(value_bin = cut(value, seq(0, 1, 0.2), include.lowest = TRUE, right = FALSE)) %>%
  left_join(pop_in_hybas4, by = "id") %>%
  left_join(msa_in_hybas4, by = "id") %>%
  left_join(hanpp_in_hybas4, by = "id") %>%
  mutate(log_pop2019 = ifelse(pop2019 != 0, log10(pop2019), 1e-6),  # avoid log(0)
         log_pop2019 = ifelse(log_pop2019 < 0, 1e-6, log_pop2019)) %>% # population < 1
  rename(msa_2015 = TerrestrialMSA_2015_World,
         scenario = facet)

zenodo_export_auxdata <- percrank_bins %>%
  filter(scenario != "median_percrank_dhf_crf") %>%
  select(-c(isIce, isNodata)) %>%
  mutate(across(-c(id, scenario, value_bin), ~ ifelse(is.na(value_bin), NA, .)),
         scenario = case_when(
           scenario == "median_percrank_dhf" ~ "dhf_only",
           scenario == "median_percrank_crf" ~ "crf_only"
         ),
         areal_division = ardiv) %>%
  rename(area_id = id,
         median_percrank = value,
         percrank_bin = value_bin,
         population_2019 = pop2019,
         log_population_2019 = log_pop2019,
         hanpp_2010 = hanpp_perc) %>%
  relocate(log_population_2019, .after = population_2019) %>%
  relocate(areal_division, .before = 1) %>%
  arrange(area_id, scenario)

# # # zenodo export # # #
write.table(zenodo_export_auxdata, paste0(exportdir, "/fig6cde_boxplots_figS8_auxdata.csv"),
            sep = ";", row.names = FALSE)

percrank_long <- percrank_bins %>%
  filter(!is.na(value_bin) & !is.na(log_pop2019) & !is.na(msa_2015) & !is.na(hanpp_perc) &
         scenario != "median_percrank_dhf_crf") %>%
  select(value_bin, scenario, log_pop2019, msa_2015, hanpp_perc) %>%
  pivot_longer(-c(value_bin, scenario), values_to = "scatter_value", names_to = "scatter_variable")

plt_boxes <- percrank_long %>%
  group_by(value_bin, scenario, scatter_variable) %>%
  summarise(p5 = quantile(scatter_value, 0.1),
            p25 = quantile(scatter_value, 0.25),
            p50 = quantile(scatter_value, 0.5),
            p75 = quantile(scatter_value, 0.75),
            p95 = quantile(scatter_value, 0.9)) %>%
  ungroup()

scale_fixer <- percrank_long %>%
  group_by(scatter_variable, scenario) %>%
  summarise() %>%
  ungroup() %>%
  mutate(scatter_value = c(0, 1, 0, 8, 0, 1),
         value_bin = factor("(1,dummy)")) %>%
  select(value_bin, scenario, scatter_variable, scatter_value)

xlabs <- plt_boxes %>%
  select(value_bin, scenario, scatter_variable) %>%
  bind_rows(scale_fixer %>% select(-scatter_value)) %>%
  arrange(scatter_variable, value_bin, scenario) %>%
  pull(value_bin)

# sample sizes
tableS1_sample_sizes <- percrank_long %>%
  group_by(value_bin, scenario) %>%
  summarise(nn = n() / 3) %>%
  arrange(scenario, value_bin) %>%
  ungroup()

write.table(tableS1_sample_sizes,
            paste0(outdir, "/tableS1_boxplot_sample_sizes.csv"),
            sep = ";", row.names = FALSE)

plt_percrank_boxplot <- ggplot() +
  geom_boxplot(data = plt_boxes,
               aes(x = interaction(scenario, value_bin, scatter_variable),
                   ymin = p5, lower = p25, middle = p50, upper = p75, ymax = p95,
                   fill = scenario),
               stat = "identity",
               outliers = FALSE,
               color = "grey50") +
  geom_point(data = scale_fixer, aes(x = value_bin, y = scatter_value, fill = scenario)) +
  facet_wrap(~ scatter_variable,
             scales = "free") +
  scale_x_discrete(labels = xlabs) +
  scale_fill_manual(values = c("#f3a44c", "#006086")) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggsave(paste0(outdir, "/fig6_", ardiv, "_boxplots.pdf"), plt_percrank_boxplot,
       width = 20, height = 10, units = "cm")

# auxiliary variable maps ------------------------------------------------------

auxvars <- c("log_pop2019", "msa_2015", "hanpp_perc")

pop_map <- geoms %>%
  left_join(percrank_bins %>%
              select(id, log_pop2019) %>%
              distinct() %>%
              mutate(log_pop2019 = ifelse(id %in% ice_catchments, NA, log_pop2019)),
            by = "id") %>%
  tm_shape() +
  tm_polygons(fill = "log_pop2019",
              fill.scale = tm_scale_intervals(
                breaks = c(0, 4, 5, 6, 7, Inf),
                values = "brewer.purples",
                value.na = "#e4e3de",
                label.na = "ice or missing data"),
              lwd = 0) +
  tm_crs("+proj=robin") +
  tm_layout(frame = FALSE)

msa_map <- geoms %>%
  left_join(percrank_bins %>%
              select(id, msa_2015) %>%
              distinct() %>%
              mutate(msa_2015 = ifelse(id %in% ice_catchments, NA, msa_2015)),
            by = "id") %>%
  tm_shape() +
  tm_polygons(fill = "msa_2015",
              fill.scale = tm_scale_intervals(
                breaks = c(0, 0.3, 0.4, 0.6, 0.8, Inf),
                values = "carto.blu_grn",
                value.na = "#e4e3de",
                label.na = "ice or missing data"),
              lwd = 0) +
  tm_crs("+proj=robin") +
  tm_layout(frame = FALSE)

hanpp_map <- geoms %>%
  left_join(percrank_bins %>%
              select(id, hanpp_perc) %>%
              distinct() %>%
              mutate(hanpp_perc = ifelse(id %in% ice_catchments, NA, hanpp_perc)),
            by = "id") %>%
  tm_shape() +
  tm_polygons(fill = "hanpp_perc",
              fill.scale = tm_scale_intervals(
                breaks = c(-Inf, 0.1, 0.2, 0.3, 0.4, Inf),
                values = "tableau.classic_orange",
                value.na = "#e4e3de",
                label.na = "ice or missing data"),
              lwd = 0) +
  tm_crs("+proj=robin") +
  tm_layout(frame = FALSE)

tmap_save(pop_map, paste0(outdir, "/figS8_population_map.pdf"))
tmap_save(msa_map, paste0(outdir, "/figS8_msa_map.pdf"))
tmap_save(hanpp_map, paste0(outdir, "/figS8_hanpp_map.pdf"))

tmap_save(pop_map + tm_layout(legend.show = FALSE),
          paste0(outdir, "/figS8_population_map.tiff"),
          width = 8, height = 4, units = "cm", dpi = 600)
tmap_save(msa_map + tm_layout(legend.show = FALSE),
          paste0(outdir, "/figS8_msa_map.tiff"),
          width = 8, height = 4, units = "cm", dpi = 600)
tmap_save(hanpp_map + tm_layout(legend.show = FALSE),
          paste0(outdir, "/figS8_hanpp_map.tiff"),
          width = 8, height = 4, units = "cm", dpi = 600)
