# ------------------------------------------------------------------------------
#
# plot00_functions
#
# Provide helper functions common for plotting and other derivation of final
# outputs. Should only be sourced, not called with parameters.
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

.get_scenario_differences <- function(scenarios_data, col_diff_to) {

  comp_to <- scenarios_data[col_diff_to] %>% pull(1)

  y1 <- min(scenarios_data$year)
  y2 <- max(scenarios_data$year)
  cls <- unique(scenarios_data$class)

  cols_to_compare <- scenarios_data %>%
    select(-c(area, year, class, !!col_diff_to))

  nmcomp <- colnames(cols_to_compare)
  coll <- tibble()
  for (i in 1:length(nmcomp)) {

    scen_comp <- cols_to_compare[nmcomp[i]] %>% pull(1)

    mean_from <- mean(scen_comp)
    mean_to <- mean(comp_to)
    median_from <- median(scen_comp)
    median_to <- median(comp_to)
    sd_to <- sd(comp_to)

    diff_means <- mean_from - mean_to
    diff_medians <- median_from - median_to

    wilcox <- wilcox.test(scen_comp, comp_to,
                          paired = TRUE, alternative = "two.sided", exact = FALSE)
    ttest <- t.test(scen_comp, comp_to,
                    paired = TRUE, alternative = "two.sided")

    testres <- c(cls, nmcomp[i], col_diff_to, y1, y2,
                 mean_from, mean_to, median_from, median_to, sd_to,
                 diff_means, diff_medians,
                 as.numeric(wilcox["p.value"]), as.numeric(ttest["p.value"])) %>%
      setNames(c("class", "scenario_from", "scenario_to", "diff_begin_year", "diff_end_year",
                 "mean_from", "mean_to", "median_from", "median_to", "sd_to",
                 "diff_means", "diff_medians",
                 "wilcox_pval", "ttest_pval"))

    coll <- coll %>% bind_rows(testres)

  }

  coll <- coll %>%
    mutate(across(-c(class, starts_with("scenario")), ~ as.numeric(.x)))

  return (coll)

}

.plt_land_area_with_deviations <- function(plt_data, nw_data = NULL, ymax = NULL) {

  ymax <- ifelse(is.null(ymax), ceiling(max(plt_data$IQR_max) * 100) / 100, ymax)
  yby <- 0.05
  yaxis = seq(0, ymax, yby)

  xmin <- min(plt_data$year)
  xmax <- max(plt_data$year)
  xby <- 10
  xaxis <- sort(c(xmin, xmax, 2005, seq(round(xmin+xby, -1), xmax, xby)))

  xpointers <- plt_data %>%
    filter(year %in% xaxis) %>%
    select(year, scenario, rollmean_median)

  if (any(grepl("ssp", unique(plt_data$scenario)))) {
    pal <- c("#bb5b0f", "#575756", "#f794b5", "#d8527c", "#9a133d")
    xpointers <- xpointers %>% filter(scenario %in% c("historical_histsoc_default", "ssp585_2015soc_default"))
  } else {
    pal <- c("#575756", "#006085", "#f3a44c", "#bb5b0f")
    xpointers <- xpointers %>% filter(scenario == "obsclim_histsoc_default")
  }

  subplt_fig1_data <- plt_data %>%
    filter(scenario %in% c("counterclim_1901soc_default", "obsclim_histsoc_default"))

  subplt_fig1 <- ggplot() +
    geom_hline(yintercept = seq(0, ymax, 0.05)) +
    geom_hline(yintercept = unique(plt_data$boundary_upper_end), col = "#e51aba") +
    geom_hline(yintercept = unique(plt_data$boundary_baseline), col = "#e51aba", linetype = "dashed") +
    geom_segment(data = xpointers,
                 aes(x = year, xend = year, y = 0, yend = rollmean_median)) +
    geom_ribbon(data = nw_data,
                aes(x = year, ymin = nw_rollmean_iqrmin, ymax = nw_rollmean_iqrmax, fill = scenario),
                alpha = 0.5) +
    geom_ribbon(data = subplt_fig1_data,
                aes(x = year, ymin = IQR_min, ymax = IQR_max, fill = scenario),
                alpha = 0.5) +
    geom_line(data = subplt_fig1_data,
              aes(x = year, y = ensemble_median, color = scenario), lwd = 0.5) +
    geom_line(data = subplt_fig1_data,
              aes(x = year, y = rollmean_median, color = scenario), lwd = 1) +
    geom_line(data = nw_data,
              aes(x = year, y = nw_rollmean_median, color = scenario), lwd = 0.75) +
    scale_fill_manual(values = c("grey80", "salmon", "grey80")) +
    scale_colour_manual(values = c("#575756", "#0fcad4", "#bb5b0f")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 8)) +
    scale_y_continuous(limits = c(0, ymax),
                       expand = c(0,0),
                       breaks = yaxis,
                       labels = yaxis * 100) +
    scale_x_continuous(limits = c(xmin, xmax),
                       expand = c(0,0),
                       breaks = xaxis)

  if (!"ensmem" %in% colnames(plt_data)) {
    plt_data$ensmem <- "ensemble" # dummy
  }

  subplt_fig2 <- ggplot(data = plt_data) +
    geom_hline(yintercept = seq(0, ymax, 0.05)) +
    geom_hline(yintercept = unique(plt_data$boundary_upper_end), col = "#e51aba") +
    geom_hline(yintercept = unique(plt_data$boundary_baseline), col = "#e51aba", linetype = "dashed") +
    geom_segment(data = xpointers,
                 aes(x = year, xend = year, y = 0, yend = rollmean_median)) +
    geom_ribbon(aes(x = year, ymin = rollmean_iqrmin, ymax = rollmean_iqrmax, fill = scenario),
                alpha = 0.5) +
    geom_line(aes(x = year, y = rollmean_median, color = scenario, linetype = ensmem)) +
    facet_wrap(~ class) +
    scale_fill_manual(values = rep("grey80", length(unique(plt_data$scenario)))) +
    scale_colour_manual(values = pal) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 8)) +
    scale_y_continuous(limits = c(0, ymax),
                       expand = c(0,0),
                       breaks = yaxis,
                       labels = yaxis * 100) +
    scale_x_continuous(limits = c(xmin, xmax),
                       expand = c(0,0),
                       breaks = xaxis)

  subplt_figS10S11 <- ggplot(data = plt_data) +
    geom_line(aes(x = year, y = ensemble_median, color = scenario, linetype = scenario)) +
    facet_wrap(~ class) +
    theme_minimal() +
    scale_x_continuous(limits = c(xmin, xmax),
                       expand = c(0,0),
                       breaks = c(seq(1920, 2000, 20), 2019)) +
    geom_hline(yintercept = unique(plt_data$boundary_baseline), col = "#e51aba", linetype = "dashed")

  return (list(subplt_fig1, subplt_fig2, subplt_figS10S11) %>%
            setNames(c("fig1", "fig2", "figS10S11")))

}

.get_trend <- function(trend_data) {

  ts_res <- zyp.trend.vector(trend_data$ensemble_median, method = "zhang")

  ts_resid <- trend_data %>%
    select(year, ensemble_median) %>%
    mutate(estimate = ts_res["intercept"] + ts_res["trend"] * row_number(.),
           resid = ensemble_median - estimate)

  ret <- trend_data %>%
    select(area, scenario, class) %>%
    distinct() %>%
    mutate(trend_begin_year = min(trend_data$year),
           trend_end_year = max(trend_data$year),
           theil_sen_trend_ppyr = ts_res["trend"] * 100, # from 0-1 fraction to pp/year
           theil_sen_total = ts_res["trendp"] * 100,
           ols_trend_ppyr = ts_res["linear"] * 100,
           kendall_tau = ts_res["tau"],
           kendall_pval = ts_res["sig"],
           sd_vals = sd(trend_data$ensemble_median) * 100,
           sd_resid = sd(ts_resid$resid) * 100,
           shapiro_wilk_pval_resid = as.numeric(shapiro.test(ts_resid$resid)["p.value"]))

  return (ret)

}

.nmfun <- function(values, ...) {return (values)}

.ee_sum <- function(sel_rast, sel_geom) {

  ret <- exact_extract(sel_rast,
                       sel_geom,
                       fun = "sum",
                       append_cols = "id",
                       colname_fun = .nmfun,
                       force_df = TRUE,
                       max_cells_in_memory = 2^31-1)
  return (ret)

}

.ee_area_weighted_mean <- function(sel_rast, sel_geom) {

  multipart_geom_ids <- sel_geom %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(nn = n()) %>%
    filter(nn > 1) %>%
    pull(id)

  sel_geom <- sel_geom %>%
    mutate(newid = case_when(
      id %in% multipart_geom_ids ~ id + max(.$id) + row_number(.),
      TRUE ~ id
    ))

  multipart_newids <- sel_geom %>%
    filter(id %in% multipart_geom_ids) %>%
    pull(newid)

  multipart_areas <- sel_rast %>%
    cellSize(unit = "km") %>%
    .ee_sum(sel_geom %>%
              filter(id %in% multipart_geom_ids) %>%
              select(newid) %>%
              rename(id = newid))

  ret <- exact_extract(sel_rast,
                       sel_geom,
                       fun = "weighted_mean",
                       weights = "area",
                       append_cols = c("id", "newid"),
                       colname_fun = .nmfun,
                       force_df = TRUE,
                       max_cells_in_memory = 2^31-1)

  multipart_awas <- tibble()
  for (i in 1:length(multipart_geom_ids)) {

    multiparts <- ret %>%
      filter(id == multipart_geom_ids[i]) %>%
      left_join(multipart_areas, by = c("newid" = "id"))
    awa <- sum(multiparts[names(sel_rast)] * multiparts$area) / sum(multiparts$area)
    multipart_awas <- multipart_awas %>%
      bind_rows(tibble(id = multipart_geom_ids[i], awa = awa)) %>%
      rename(!!names(sel_rast) := awa)

  }

  ret <- ret %>%
    filter(!id %in% multipart_geom_ids) %>%
    select(-newid) %>%
    bind_rows(multipart_awas)

  return (ret)

}
