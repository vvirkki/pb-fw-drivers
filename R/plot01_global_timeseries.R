# ------------------------------------------------------------------------------
#
# plot01_global_timeseries
#
# Create plots and outputs related to Fig. 1, Fig. 2, Fig. S1, Fig. S2, Fig. S9,
# Fig. S10, Fig. S11, and Fig. S12.
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

source("R/libraries.R")
source("R/plot00_functions.R")

# Setup ------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
variable <- args[1]
indir <- args[2] # where data is read from
outdir <- args[3] # where outputs are saved
if (!dir.exists(outdir)) { dir.create(outdir) }

rmyrs <- 10 # number of years for rollmean
trendyrs <- 30 # number of years for assessing average pp/yr trends
diffyrs <- 30 # number of years for assessing differences between scenarios
supp <- ifelse(grepl("supplementary", indir), TRUE, FALSE) # only for Fig. S12

exportdir <- "zenodo_export"
if (!dir.exists(exportdir)) { dir.create(exportdir) }
exportvar <- case_when(
  variable == "dis" ~ "streamflow",
  variable == "rootmoist" ~ "soilmoisture"
)

# Global land area with deviations ---------------------------------------------

global_files <- list.files(paste0("output/", variable, "_", indir),
                           recursive = TRUE, full.names = TRUE)

files_in <- global_files[grepl("global_land_area_with_local_deviations_annual_mean_ensemble_median_IQR.csv",
                               global_files)]

data_in <- files_in %>%
  lapply(read_delim, delim = ";", show_col_types = FALSE) %>%
  bind_rows() %>%
  distinct() %>%
  group_by(scenario, class) %>%
  mutate(rollmean_median = rollapply(ensemble_median, rmyrs, mean, align = "right", partial = TRUE),
         rollmean_iqrmin = rollapply(IQR_min, rmyrs, mean, align = "right", partial = TRUE),
         rollmean_iqrmax = rollapply(IQR_max, rmyrs, mean, align = "right", partial = TRUE)) %>%
  ungroup() %>%
  filter(year >= 1911) # exclude 10 first years due to spinups

# REQUIRED EXTERNAL DATA:
# PB-FW estimates by Porkka et al. (2024) (https://doi.org/10.1038/s44221-024-00208-7)
# https://doi.org/10.5281/zenodo.10531807

nat_water_var <- ifelse(variable == "dis", "streamflow", "soilmoisture")
nat_water_data <- read_delim(paste0("Data/Porkka_et_al_NatWater_2024_data/Fig2_Fig4_Fig5_EDFig1_EDFig2ab/global/",
                                    nat_water_var, "/global_land_area_with_local_deviations_annual_mean_ensemble_median_IQR.csv"),
                             delim = ",", show_col_types = FALSE) %>%
  filter(class == "summed_deviations" & year >= 1911) %>%
  mutate(rollmean_median = rollapply(ensemble_median, rmyrs, mean, align = "right", partial = TRUE),
         rollmean_iqrmin = rollapply(IQR_min, rmyrs, mean, align = "right", partial = TRUE),
         rollmean_iqrmax = rollapply(IQR_max, rmyrs, mean, align = "right", partial = TRUE),
         scenario = "nat_water") %>%
  rename_with(~ paste0("nw_", .x), -c(area, year, class, scenario))

# historical timeseries plots --------------------------------------------------

# summed deviations
data_fig1fig2_summed_deviations <- data_in %>%
  filter(class == "summed_deviations")

# dry and wet deviations separately
data_fig2_separate_deviations <- data_in %>%
  filter(class != "summed_deviations")

fig1_summed_deviations <- .plt_land_area_with_deviations(data_fig1fig2_summed_deviations,
                                                         nat_water_data, ymax = 0.35)$fig1
fig2_summed_deviations <- .plt_land_area_with_deviations(data_fig1fig2_summed_deviations,
                                                         ymax = 0.36)$fig2
fig2_separate_deviations <- .plt_land_area_with_deviations(data_fig2_separate_deviations,
                                                           ymax = 0.2)$fig2

# pdf outputs ------------------------------------------------------------------

if (!supp) {
  ggsave(paste0(outdir, "/fig1_", variable, ".pdf"),
         fig1_summed_deviations,
         width = 20, height = 10, units = "cm")
}

supp_spec <- ifelse(supp, "for_figS12_", "")
ggsave(paste0(outdir, "/fig2_", supp_spec, variable, ".pdf"),
       arrangeGrob(fig2_summed_deviations, fig2_separate_deviations),
       width = 20, height = 10, units = "cm")

# csv prep ---------------------------------------------------------------------

fig1fig2_numbers <- data_fig1fig2_summed_deviations %>%
  bind_rows(data_fig2_separate_deviations) %>%
  group_by(scenario, class) %>%
  filter(year == max(year)) %>%
  select(area, year, scenario, class, rollmean_median, starts_with("boundary")) %>%
  mutate(across(c(rollmean_median, starts_with("boundary")), ~ round(.x * 100, 3))) %>%
  arrange(class, year)

data_fig1fig2_trends <- data_fig1fig2_summed_deviations %>%
  bind_rows(data_fig2_separate_deviations)

trend_ends <- c(2019)
scenarios <- unique(data_fig1fig2_trends$scenario)
classes <- unique(data_fig1fig2_trends$class)

coll_trends <- tibble()
coll_scenario_diffs <- tibble()

for (i in 1:length(trend_ends)) {
  for (k in 1:length(classes)) {
    for (j in 1:length(scenarios)) {

      trend_data <- data_fig1fig2_trends %>%
        filter(year > (trend_ends[i] - trendyrs) & year <= trend_ends[i] &
               scenario == scenarios[j] & class == classes[k])

      linear_trend <- .get_trend(trend_data)
      coll_trends <- coll_trends %>%
        bind_rows(linear_trend)

    }

    scenario_difference <- data_fig1fig2_trends %>%
      select(area, year, scenario, class, ensemble_median) %>%
      filter(year > (trend_ends[i] - diffyrs) & year <= trend_ends[i] &
               class == classes[k]) %>%
      pivot_wider(id_cols = c(area, year, class),
                  names_from = scenario, values_from = ensemble_median)

    scendiff <- .get_scenario_differences(scenario_difference,
                                          "counterclim_1901soc_default")

    coll_scenario_diffs <- coll_scenario_diffs %>%
      bind_rows(scendiff)
  }
}

fig1_fig2_trends <- coll_trends %>%
  arrange(scenario, class, trend_begin_year)

fig1_fig2_scendiffs <- coll_scenario_diffs %>%
  arrange(scenario_from, class, diff_begin_year)

# csv outputs ------------------------------------------------------------------

write.table(fig1fig2_numbers, paste0(outdir, "/fig1fig2_numbers_",  supp_spec, variable, ".csv"),
            sep = ";", row.names = FALSE)

if (supp) {

  # # # zenodo export # # #
  write.table(data_in %>%
                select(-starts_with("rollmean"), -ts_normalised) %>%
                mutate(variable = exportvar, .before = 1,
                       scenario = case_when(
                         scenario == "counterclim_1901soc_default" ~ "baseline",
                         scenario == "obsclim_histsoc_default" ~ "historical",
                         scenario == "counterclim_histsoc_default" ~ "dhf_only",
                         scenario == "obsclim_1901soc_default" ~ "crf_only"
                       )) %>%
                rename(boundary_median = boundary_baseline,
                       value_ensemble_median = ensemble_median,
                       ensemble_IQR_min = IQR_min,
                       ensemble_IQR_max = IQR_max),
              paste0(exportdir, "/figS12_", exportvar, ".csv"), sep = ";", row.names = FALSE)
  quit(save = "no")

}

write.table(fig1_fig2_trends, paste0(outdir, "/fig1fig2_ppyr_trends_", variable, ".csv"),
            sep = ";", row.names = FALSE)

write.table(fig1_fig2_scendiffs, paste0(outdir, "/fig1fig2_scenario_differences_", variable, ".csv"),
            sep = ";", row.names = FALSE)

# # # zenodo export # # #
write.table(data_in %>%
              select(-starts_with("rollmean"), -ts_normalised) %>%
              mutate(variable = exportvar, .before = 1,
                     scenario = case_when(
                       scenario == "counterclim_1901soc_default" ~ "baseline",
                       scenario == "obsclim_histsoc_default" ~ "historical",
                       scenario == "counterclim_histsoc_default" ~ "dhf_only",
                       scenario == "obsclim_1901soc_default" ~ "crf_only"
                     )) %>%
              rename(boundary_median = boundary_baseline,
                     value_ensemble_median = ensemble_median,
                     ensemble_IQR_min = IQR_min,
                     ensemble_IQR_max = IQR_max),
            paste0(exportdir, "/fig1_fig2_", exportvar, ".csv"), sep = ";", row.names = FALSE)

# Ensemble member wise land area with deviations -------------------------------

ensmem_files <- list.files(paste0("output/", variable, "_main_single_ensemble_members"),
                           recursive = TRUE, full.names = TRUE)

# do all ensmems that are available in outputs
files_in <- ensmem_files[grepl("global_land_area_with_local_deviations_annual_mean_ensemble_median_IQR.csv",
                               ensmem_files)]

ensmems <- files_in %>%
  lapply(str_split, "/") %>%
  lapply(unlist) %>%
  lapply(function(x) {x[4]}) %>%
  unlist() %>%
  unique()

coll_plots <- vector(mode = "list", length = length(ensmems))
coll_ensmems <- vector(mode = "list", length = length(ensmems))
for (i in 1:length(ensmems)) {

  ensmem_data <- files_in[grepl(ensmems[i], files_in)] %>%
    lapply(read_delim, delim = ";", show_col_types = FALSE) %>%
    bind_rows() %>%
    distinct() %>%
    group_by(scenario, class) %>%
    mutate(rollmean_median = rollapply(ensemble_median, rmyrs, mean, align = "right", partial = TRUE),
           rollmean_iqrmin = rollapply(IQR_min, rmyrs, mean, align = "right", partial = TRUE),
           rollmean_iqrmax = rollapply(IQR_max, rmyrs, mean, align = "right", partial = TRUE)) %>%
    ungroup() # %>%
    # filter(year >= 1911) # draw also 10 first years in these

  coll_ensmems[[i]] <- ensmem_data %>%
    mutate(ensmem = ensmems[i])

  plt_data_ensmem <- ensmem_data %>%
    filter(class == "summed_deviations" &
           scenario %in% c("counterclim_1901soc_default", "obsclim_histsoc_default"))

  coll_plots[[i]] <- .plt_land_area_with_deviations(plt_data_ensmem)$figS10S11 +
    ggtitle(ensmems[i]) +
    scale_color_manual(values = c("#575756", "#bb5b0f")) +
    scale_linetype_manual(values = c("solid", "solid")) +
    scale_y_continuous(limits = c(0, 0.412),
                       breaks = seq(0, 0.3, 0.1),
                       labels = c(0, 10, 20, 30),
                       expand = c(0,0)) +
    scale_x_continuous(breaks = c(1901, 1930, 1960, 1990, 2019),
                       expand = c(0,0)) +
    theme(strip.text = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 8),
          legend.position = "none",
          panel.grid.minor = element_blank())

}

arrangeGrob(grobs = coll_plots, nrow = 5, ncol = 3) %>%
  ggsave(paste0(outdir, "/figS10S11_", variable, ".pdf"), .,
         width = 15, height = 18, units = "cm")

# # # zenodo export # # #
coll_ensmems %>%
  bind_rows() %>%
  filter(scenario %in% c("counterclim_1901soc_default", "obsclim_histsoc_default") &
           class == "summed_deviations") %>%
  mutate(variable = exportvar,
         scenario = case_when(
           scenario == "counterclim_1901soc_default" ~ "baseline",
           scenario == "obsclim_histsoc_default" ~ "historical"
         )) %>%
  select(variable, ensmem, area, year, scenario, class, ensemble_median, boundary_baseline) %>%
  rename(boundary_median = boundary_baseline,
         ensemble_member = ensmem,
         value = ensemble_median) %>%
  write.table(paste0(exportdir, "/figS10S11_", exportvar, ".csv"),
              sep = ";", row.names = FALSE)

# ensemble member comparisons for supplementary figures ------------------------

ensmem_filter <- case_when(
  variable == "dis" ~ "watergap",
  variable == "rootmoist" ~ "miroc-integ-land"
)

data_figS1S2_summed_deviations <- coll_ensmems %>%
  bind_rows() %>%
  filter(grepl(ensmem_filter, ensmem) & class == "summed_deviations")

data_figS1S2_separate_deviations <- coll_ensmems %>%
  bind_rows() %>%
  filter(grepl(ensmem_filter, ensmem) & class != "summed_deviations")

figS1S2_summed_deviations <- .plt_land_area_with_deviations(data_figS1S2_summed_deviations)$fig2 +
  scale_linetype_manual(values = c("dashed", "dotdash", "dotted"))
figS1S2_summed_deviations$layers[c(2,3,5)] <- NULL # remove ensemble range, boundaries

figS1S2_separate_deviations <- .plt_land_area_with_deviations(data_figS1S2_separate_deviations)$fig2 +
  scale_linetype_manual(values = c("dashed", "dotdash", "dotted")) +
  theme(legend.position = "none")
figS1S2_separate_deviations$layers[c(2,3,5)] <- NULL # remove ensemble range, boundaries

if (variable == "dis") {

  nat_water_data <- read_delim(paste0("Data/Porkka_et_al_NatWater_2024_data/EDFig3_EDFig4/WaterGAP2/streamflow/",
                                      "global_land_area_with_local_deviations_annual_mean_ensemble_median_IQR.csv"),
                               delim = ",", show_col_types = FALSE) %>%
    filter(year >= 1911) %>%
    group_by(class) %>%
    mutate(rollmean_median = rollapply(ensemble_median, rmyrs, mean, align = "right", partial = TRUE),
           rollmean_iqrmin = rollapply(IQR_min, rmyrs, mean, align = "right", partial = TRUE),
           rollmean_iqrmax = rollapply(IQR_max, rmyrs, mean, align = "right", partial = TRUE),
           scenario = "nat_water") %>%
    rename_with(~ paste0("nw_", .x), -c(area, year, class, scenario))

  figS1S2_summed_deviations <- figS1S2_summed_deviations +
    geom_ribbon(data = nat_water_data %>% filter(class == "summed_deviations"),
                aes(x = year, ymin = nw_rollmean_iqrmin, ymax = nw_rollmean_iqrmax, fill = scenario),
                alpha = 0.5)

  figS1S2_separate_deviations <- figS1S2_separate_deviations +
    geom_ribbon(data = nat_water_data %>% filter(class != "summed_deviations"),
                aes(x = year, ymin = nw_rollmean_iqrmin, ymax = nw_rollmean_iqrmax, fill = scenario),
                alpha = 0.5) +
    facet_wrap(~ class)

}

ggsave(paste0(outdir, "/figS1S2_", variable, ".pdf"),
       arrangeGrob(figS1S2_summed_deviations, figS1S2_separate_deviations),
       width = 20, height = 10, units = "cm")

# # # zenodo export # # #
data_figS1S2_summed_deviations %>%
  bind_rows(data_figS1S2_separate_deviations) %>%
  filter(year >= 1911) %>%
  mutate(variable = exportvar,
         scenario = case_when(
           scenario == "counterclim_1901soc_default" ~ "baseline",
           scenario == "obsclim_histsoc_default" ~ "historical",
           scenario == "counterclim_histsoc_default" ~ "dhf_only",
           scenario == "obsclim_1901soc_default" ~ "crf_only"
         )) %>%
  select(variable, ensmem, area, year, scenario, class, ensemble_median) %>%
  rename(ensemble_member = ensmem,
         value = ensemble_median) %>%
  write.table(paste0(exportdir, "/figS1S2_", exportvar, ".csv"),
              sep = ";", row.names = FALSE)

# future projections plots -----------------------------------------------------

future_files <- list.files(paste0("output/", variable, "_future_ensemble"),
                                  recursive = TRUE, full.names = TRUE)

future_files_in <- future_files[grepl("global_land_area_with_local_deviations_annual_mean_ensemble_median_IQR.csv",
                                      future_files)]

future_data_in <- future_files_in %>%
  lapply(read_delim, delim = ";", show_col_types = FALSE) %>%
  bind_rows() %>%
  distinct() %>%
  group_by(scenario, class) %>%
  mutate(rollmean_median = rollapply(ensemble_median, rmyrs, mean, align = "right", partial = TRUE),
         rollmean_iqrmin = rollapply(IQR_min, rmyrs, mean, align = "right", partial = TRUE),
         rollmean_iqrmax = rollapply(IQR_max, rmyrs, mean, align = "right", partial = TRUE)) %>%
  ungroup() %>%
  filter(year >= 1911) # exclude 10 first years due to spinups

# future projections summed deviations
data_figS9_future_summed_deviations <- future_data_in %>%
  filter(class == "summed_deviations")

# future projections dry and wet deviations separately
data_figS9_future_separate_deviations <- future_data_in %>%
  filter(class != "summed_deviations")

figS9_future_summed_deviations <- .plt_land_area_with_deviations(data_figS9_future_summed_deviations,
                                                                 ymax = 0.5)$fig2
figS9_future_separate_deviations <- .plt_land_area_with_deviations(data_figS9_future_separate_deviations,
                                                                   ymax = 0.3)$fig2

ggsave(paste0(outdir, "/figS9_", variable, ".pdf"),
       arrangeGrob(figS9_future_summed_deviations, figS9_future_separate_deviations),
       width = 20, height = 10, units = "cm")

figS9_numbers <- data_figS9_future_summed_deviations %>%
  bind_rows(data_figS9_future_separate_deviations) %>%
  group_by(scenario, class) %>%
  filter(year == max(year)) %>%
  select(area, year, scenario, class, rollmean_median, starts_with("boundary")) %>%
  mutate(across(c(rollmean_median, starts_with("boundary")), ~ round(.x * 100, 3))) %>%
  arrange(class, year)

write.table(figS9_numbers, paste0(outdir, "/figS9_numbers_", variable, ".csv"),
            sep = ";", row.names = FALSE)

# # # zenodo export # # #
write.table(future_data_in %>%
              select(-starts_with("rollmean"), -starts_with("boundary"), -ts_normalised) %>%
              mutate(variable = exportvar, .before = 1,
                     scenario = case_when(
                       scenario == "picontrol_1850soc_default" ~ "baseline",
                       scenario == "historical_histsoc_default" ~ "historical",
                       scenario == "ssp126_2015soc_default" ~ "ssp1_rcp26",
                       scenario == "ssp370_2015soc_default" ~ "ssp3_rcp70",
                       scenario == "ssp585_2015soc_default" ~ "ssp5_rcp85"
                     )) %>%
              rename(value_ensemble_median = ensemble_median,
                     ensemble_IQR_min = IQR_min,
                     ensemble_IQR_max = IQR_max),
            paste0(exportdir, "/figS9_", exportvar, ".csv"), sep = ";", row.names = FALSE)
