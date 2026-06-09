# ------------------------------------------------------------------------------
#
# output01_land_area_with_local_deviations
#
# Determine pre-industrial variability and create outputs (plot, csv) describing
# the percentage of land area with deviations, either globally or within areas
# prescribed in config.
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

.ts_boundVoc <- c("Summed dry and wet", "Dry", "Wet") %>%
  setNames(c("summed", "below", "above"))
.ts_varVoc <- c("streamflow", "root-zone soil moisture") %>%
  setNames(c("dis", "rootmoist"))
.ts_colorLst <- list(median = "black", dash = "#545454", rollmean = "red")
.ts_lineType <- list(median = "solid", dash = "dashed", rollmean = "solid")
.ts_drawVar <- c("landShareInClass") %>%
  setNames(c("landShareInClass"))

.read_visopts <- function(path) {

  # read config file; return a list of options to be used
  cfg_data <- readLines(path) %>%
    tibble() %>%
    rename(raw = 1) %>%
    separate(raw, into = c("variable", "value"), sep = ";")

  cfg <- cfg_data$value %>%
    as.list() %>%
    setNames(cfg_data$variable)

  for (i in 1:length(cfg)) {
    if (grepl(",", cfg[[i]])) {
      cfg[[i]] <- cfg[[i]] %>%
        str_split(",") %>%
        unlist() %>%
        as.vector()
    }
    suppressWarnings(
      cfg[[i]] <- ifelse(is.na(as.numeric(cfg[[i]])),
                         cfg[[i]], as.numeric(cfg[[i]]))
    )
  }

  return (cfg)

}

.compute_transgression_year <- function(boundary_values, area_shares,
                                        draw_class, npers) {

  for (i in 1:nrow(boundary_values)) {
    b <- boundary_values[i,]$boundary
    area_shares <- area_shares %>%
      mutate(!!b := boundary_values[i,]$value) %>%
      mutate(!!paste0("tg_", b) := drawVar > .[b])
  }

  cumulative_transgressions <- area_shares %>%
    mutate(across(starts_with("tg_"),
                  ~ zoo::rollsum(x = .x, k = npers, align = "right", fill = 0))) %>%
    arrange(year)

  persistent_transgression_years <- c() # year when transgression becomes persistent
  nm <- c()
  for (i in 1:nrow(boundary_values)) {
    b <- boundary_values[i,]$boundary
    tg_yr <- NA
    j <- 1
    while (j < nrow(cumulative_transgressions)) {
      yr <- cumulative_transgressions %>%
        slice(j:nrow(cumulative_transgressions)) %>%
        pull(!!paste0("tg_", b))
      if (all(yr == npers)) {
        tg_yr <- cumulative_transgressions %>%
          slice(j) %>%
          pull(year)
        persistent_transgression_years <- c(persistent_transgression_years, tg_yr)
        nm <- c(nm, paste0("tg_", b))
        break
      } else {
        j <- j + 1
      }
    }
    if (j == nrow(cumulative_transgressions)) {
      persistent_transgression_years <- c(persistent_transgression_years, tg_yr)
      nm <- c(nm, paste0("tg_", b))
    }
  }

  ret <- tibble(dptClass = draw_class,
                # first year of transgression that became persistent
                value = (persistent_transgression_years - npers + 1),
                boundary = str_replace(nm, "tg_", "persistent_transgression_year_")) %>%
    pivot_wider(id_cols = dptClass, names_from = boundary)

  return (ret)

}

.compose_output <- function(cfg, plt_data, draw_class, visopts) {

  # prepare data
  plt_data <- plt_data %>%
    filter(dptClass == draw_class)

  if (!cfg$draw_ts_var %in% .ts_drawVar) {
    stop(paste0("draw variable should be one of: ",
                paste(names(.ts_drawVar), collapse = ", ")))
  }

  if (any(is.na(plt_data[cfg$draw_ts_var])) | length(unique(plt_data$period)) == 1) {
    write("returning empty...", cfg$log_file, append = TRUE)
    return (list(ggplot(), tibble()) %>% setNames(c("plt", "currData")))
  }

  ystart <- max(min(cfg$period_begin) + cfg$cptmax, min(plt_data$year))
  yend <- max(plt_data$year)

  plt_data <- plt_data %>%
    filter(year >= ystart) %>%
    rename(drawVar = landShareInClass)

  # axis properties
  brky_max_fixed <- visopts[[paste0("ts_yrange_", draw_class)]]
  if (!as.logical(visopts$ts_fixed_yrange)) {
    brky_max <- max(brky_max_fixed, quantile(plt_data$drawVar, 0.95))
  } else {
    brky_max <- brky_max_fixed
  }
  brky_max <- 1 # dummy
  breaks_y <- seq(0, brky_max, visopts$ts_yspacing) %>%
    round(2)

  breaks_x <- seq(ystart, yend, 1)
  breaks_x <- breaks_x[which(breaks_x %% visopts$ts_xspacing == 0)]
  breaks_x[1] <- ystart
  breaks_x[length(breaks_x)] <- yend

  # boundary lines and their labels
  boundaries <- plt_data %>%
    filter(period == cfg$baseline_period &
           year >= cfg$local_bound_begin &
           year <= cfg$local_bound_end) %>%
    group_by(year) %>%
    summarise(drawVar = median(drawVar)) %>% # ensemble median
    ungroup() %>%
    summarise(baseline = quantile(drawVar, cfg$boundary_baseline),
              upper_end = quantile(drawVar, cfg$boundary_upper_end)) %>%
    mutate(dptClass = draw_class)

  if (cfg$normalise_ts %>% as.logical()) {

    plt_data <- plt_data %>%
      mutate(abv_upper_end = drawVar >= boundaries$upper_end,
             blw_baseline = drawVar <= boundaries$baseline) %>%
      mutate(drawVar = case_when(
        drawVar == 0 ~ -1, # no departing land area
        abv_upper_end ~ drawVar * (1 / boundaries$upper_end),
        blw_baseline ~ (drawVar - boundaries$baseline) / boundaries$baseline,
        TRUE ~ 1 - ((boundaries$upper_end - drawVar) / (boundaries$upper_end - boundaries$baseline))
      )) %>%
      select(-c(abv_upper_end, blw_baseline))

    boundaries <- boundaries %>%
      mutate(baseline = 0,
             upper_end = 1)

    brky_max <- max(1, min(10, max(plt_data$drawVar))) # for visuals only
    breaks_y <- c(-1, -0.5, 0, 0.5, seq(1, brky_max, 0.5))

  }

  boundary_values <- data.frame(c("baseline", "upper_end"),
                                  c(cfg$boundary_baseline,
                                    cfg$boundary_upper_end)) %>%
    setNames(c("boundary", "position"))

  plt_data_boundaries <- boundaries %>%
    pivot_longer(-dptClass, names_to = "boundary") %>%
    mutate(xbegin = ystart,
           xend = yend)

  plt_data_boundary_labs <- plt_data_boundaries %>%
    left_join(boundary_values, by = "boundary") %>%
    mutate(label = paste0(sprintf("%.1f", round(value * 100, 1)), "%"),
           label = paste0(label, " (", position, ")"),
           xend = yend) %>%
    select(xend, value, label) %>%
    setNames(c("x", "y", "label"))

  # backgrounds
  bg_below_boundaries <- data.frame(space = "aaa_below_boundaries") %>%
    mutate(ymin = 0,
           ymax = min(plt_data_boundaries$value))

  bg_between_boundaries <- data.frame(space = "aab_between_boundaries") %>%
    mutate(ymin = min(plt_data_boundaries$value),
           ymax = max(plt_data_boundaries$value))

  bg_above_boundaries <- data.frame(space = "aac_above_boundaries") %>%
    mutate(ymin = max(plt_data_boundaries$value),
           ymax = Inf)

  plt_data_backgrounds <- bg_below_boundaries %>%
    bind_rows(bg_between_boundaries) %>%
    bind_rows(bg_above_boundaries) %>%
    mutate(xmin = ystart,
           xmax = yend)

  # interquartile range shadings
  plt_data_ribbons <- plt_data %>%
    group_by(year, period) %>%
    summarise(minrib = quantile(drawVar, 0.25),
              maxrib = quantile(drawVar, 0.75)) %>% # ensemble interquartile range
    ungroup() %>%
    mutate(ribbonid = period)

  # lines and their labels
  plt_data_lines <- plt_data %>%
    group_by(year, period) %>%
    summarise(drawVar = median(drawVar)) %>% # ensemble median
    ungroup() %>%
    mutate(lineid = paste0("median_", period))

  begin_rollmean <- max(cfg$period_begin) + cfg$cptmax - visopts$ts_moving_window_yrs + 1

  plt_data_lines_rollmean <- plt_data_lines %>%
    arrange(year) %>%
    filter(!(year >= max(cfg$period_begin) & year < begin_rollmean)) %>%
    group_by(lineid) %>%
    mutate(drawVar = c(zoo::rollapply(drawVar, visopts$ts_moving_window_yrs, mean,
                                      align = "right", partial = TRUE))) %>%
    filter(!is.na(drawVar)) %>%
    mutate(lineid = paste0("rollmean_", period)) %>%
    ungroup()

  plt_data_line_labs <- plt_data_lines_rollmean %>%
    filter(year == max(.$year)) %>%
    mutate(status = sprintf("%.1f", round(drawVar * 100, 1)),
           label = paste0(status, "% (", visopts$ts_moving_window_yrs, "-yr mean)")) %>%
    select(year, drawVar, label) %>%
    setNames(c("x", "y", "label"))

  dupl_row <- plt_data_lines %>%
    filter(year == max(cfg$period_begin) + cfg$cptmax) # to continue from dashed to solid line

  plt_data_lines <- plt_data_lines %>%
    mutate(dash = year >= max(cfg$period_begin) & year <= max(cfg$period_begin) + cfg$cptmax,
           lineid = ifelse(dash, str_replace(lineid, "median", "dash"), lineid)) %>%
    select(-dash) %>%
    bind_rows(plt_data_lines_rollmean) %>%
    bind_rows(dupl_row)

  # annotations for persistent transgression times
  upper_end_boundary <- plt_data_boundaries %>%
    filter(boundary == "upper_end") %>%
    select(value, boundary)

  area_shares_rollmean <- plt_data_lines_rollmean %>%
    filter(period == cfg$period_label[cfg$period_label != cfg$baseline_period]) %>%
    select(year, drawVar)

  transgression_year <- .compute_transgression_year(upper_end_boundary,
                                                    area_shares_rollmean,
                                                    draw_class,
                                                    visopts$ts_persistent_transgression_yrs)

  plt_data_transgression_years <- transgression_year %>%
    pivot_longer(-dptClass, values_to = "year") %>%
    filter(!is.na(year)) %>%
    left_join(plt_data_lines_rollmean %>% filter(period != cfg$baseline_period),
              by = "year") %>%
    select(year, drawVar) %>%
    rename(x = year, yend = drawVar) %>%
    mutate(xend = x, y = 0)

  # other plot elements
  plt_data_labs <- plt_data_boundary_labs %>%
    bind_rows(plt_data_line_labs)

  period_switch <- plt_data %>%
    group_by(period) %>%
    summarise(ymax = max(year)) %>%
    pull(ymax) %>%
    min()

  nensmem <- plt_data %>%
    group_by(impactmodel, forcing) %>%
    n_groups()

  pal_color <- plt_data_lines$lineid %>% unique() %>% sort()
  for (i in 1:length(.ts_colorLst)) {
    pal_color[grepl(names(.ts_colorLst[i]), pal_color)] <- .ts_colorLst[[i]]
  }

  pal_linetype <- plt_data_lines$lineid %>% unique() %>% sort()
  for (i in 1:length(.ts_lineType)) {
    pal_linetype[grepl(names(.ts_lineType[i]), pal_linetype)] <- .ts_lineType[[i]]
  }

  plt_title <- paste0(.ts_boundVoc[[draw_class]], " ",
                      .ts_varVoc[[cfg$variable]],
                      " deviations (annual mean; ensemble n=",
                      nensmem, ")")

  # plot
  if (as.logical(visopts$ts_draw_backgrounds)) {
    plt <- ggplot() +
      geom_rect(data = plt_data_backgrounds,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = space),
                alpha = 0.3)
    pal_fill <- c("#3a4f8e", "#F4B326", "#cf5426", "grey50", "grey50")
  } else {
    plt <- ggplot()
    pal_fill <- c("grey50", "grey50")
  }

  plt <- plt +
    geom_ribbon(data = plt_data_ribbons,
                aes(x = year, ymin = minrib, ymax = maxrib, fill = ribbonid),
                alpha = 0.5) +
    geom_segment(data = plt_data_boundaries,
                 aes(x = xbegin, xend = xend, y = value, yend = value),
                 color = "grey25",
                 linetype = "dashed",
                 linewidth = (visopts$ts_guideline_width / ggplot2::.pt)) +
    geom_line(data = plt_data_lines,
              aes(x = year, y = drawVar, color = lineid, linetype = lineid),
              linewidth = 1.25 * (visopts$ts_guideline_width / ggplot2::.pt)) +
    geom_text(data = plt_data_transgression_years,
              aes(x = x, y = brky_max * 0.03, label = x),
              size = (visopts$ts_label_size / ggplot2::.pt),
              hjust = 0,
              nudge_x = 2) +
    geom_text(data = plt_data_labs,
              aes(x = x, y = y, label = label),
              size = (visopts$ts_label_size / ggplot2::.pt),
              hjust = 0) +
    geom_vline(xintercept = period_switch,
               color = "grey25",
               linewidth = (visopts$ts_guideline_width / ggplot2::.pt)) +
    geom_segment(aes(x = ystart, xend = yend, y = 0, yend = 0),
                 color = "grey25",
                 linewidth = (visopts$ts_guideline_width / ggplot2::.pt)) +
    geom_segment(data = plt_data_transgression_years,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "grey25",
                 linewidth = (visopts$ts_guideline_width / ggplot2::.pt) / 2) +
    scale_color_manual(values = pal_color) +
    scale_fill_manual(values = pal_fill) +
    scale_linetype_manual(values = pal_linetype) +
    scale_y_continuous(name = "Percentage of ice-free land area",
                       breaks = breaks_y,
                       # labels = abs(breaks_y * 100),
                       expand = expansion()) +
    scale_x_continuous(breaks = breaks_x,
                       limits = c(ystart, yend),
                       expand = expansion(mult = c(0, 0.15))) +
    # coord_cartesian(ylim = c(min(breaks_y), brky_max)) +
    ggtitle(plt_title) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = visopts$ts_title_size),
          axis.text = element_text(size = visopts$ts_axistext_size),
          axis.title = element_text(size = visopts$ts_annotation_size))

  # collect data on boundaries and current status for exporting
  current_status <- plt_data_lines_rollmean %>%
    filter(year == yend & period != cfg$baseline_period) %>%
    select(drawVar) %>%
    mutate(boundary = "current") %>%
    rename(value = drawVar)

  export_boundaries_current_status <- plt_data_boundaries %>%
    filter(boundary %in% c("baseline", "upper_end")) %>%
    select(value, boundary) %>%
    bind_rows(current_status) %>%
    mutate(dptClass = draw_class) %>%
    pivot_wider(id_cols = dptClass, names_from = boundary) %>%
    mutate(relStatus_baseline = current / baseline,
           relStatus_upper_end = current / upper_end,
           area = unique(plt_data$area)) %>%
    left_join(transgression_year, by = "dptClass") %>%
    rename(class = dptClass) %>%
    mutate(class = case_when(
             draw_class == "summed" ~ "summed_deviations",
             draw_class == "above" ~ "wet_deviations",
             draw_class == "below" ~ "dry_deviations")) %>%
    select(area, everything())

  export_ensemble_median_iqr <- plt_data_lines %>%
    filter(!grepl("rollmean_", lineid)) %>%
    filter(!grepl("dash_", lineid)) %>% # only when cptmin = cptmax = 0
    select(year, period, drawVar) %>%
    mutate(class = case_when(
             draw_class == "summed" ~ "summed_deviations",
             draw_class == "above" ~ "wet_deviations",
             draw_class == "below" ~ "dry_deviations"
           ),
           area = unique(plt_data$area)) %>%
    left_join(plt_data_ribbons %>% select(year, period, minrib, maxrib),
              by = c("year", "period")) %>%
    rename(scenario = period,
           ensemble_median = drawVar,
           IQR_min = minrib,
           IQR_max = maxrib) %>%
    mutate(boundary_baseline = export_boundaries_current_status$baseline,
           boundary_upper_end = export_boundaries_current_status$upper_end,
           ts_normalised = cfg$normalise_ts) %>%
    select(area, year, scenario, class, everything()) %>%
    arrange(scenario, year)

  return (list(plt,
               export_boundaries_current_status,
               export_ensemble_median_iqr) %>%
            setNames(c("plt",
                       "export_boundaries_current_status",
                       "export_ensemble_median_iqr")))

}

# public function
output_land_area_with_local_deviations <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("outputting land area with local deviations...",
        cfg$log_file, append = TRUE)

  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_") %>%
    paste(collapse = "|")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_") %>%
    paste(collapse = "|")
  ardiv_lbl <- ifelse(cfg$areal_division == "", "global", cfg$areal_division)

  visopts_cfg <- paste0("configs/", cfg$visopts_cfg)
  if (!file.exists(visopts_cfg)) {
    stop("visopts not defined for this config")
  } else {
    visopts <- .read_visopts(visopts_cfg)
  }

  if (!cfg$areal_division == "") {

    ardiv_rast_exists <- file.exists(paste0("data/ardiv/", cfg$areal_division, ".tif"))
    ardiv_poly_exists <- file.exists(paste0("data/ardiv/", cfg$areal_division, ".gpkg"))

    if (!ardiv_rast_exists & !ardiv_poly_exists) {
      stop(paste0("neither polys nor rast for areal division ",
                  cfg$areal_division, " found"))
    }
    if (ardiv_rast_exists) {
      unique_areas <- rast(paste0("data/ardiv/", cfg$areal_division, ".tif")) %>%
        values() %>%
        unique() %>%
        as.numeric()
    }
    if (ardiv_poly_exists) {
      unique_areas <- read_sf(paste0("data/ardiv/", cfg$areal_division, ".gpkg")) %>%
        pull(id) %>%
        unique()
    }
    unique_areas <- unique_areas[!is.na(unique_areas)]

  } else {
    unique_areas <- c(NA)
  }

  # create directories
  tstamp <- Sys.time()
  output_fldr <- paste0("output/", cfg$variable, "_", ardiv_lbl,
                     "_land_area_with_local_deviations_",
                     paste(cfg$impactmodel, collapse = "_"), "_", tstamp)
  if (!dir.exists(output_fldr)) {
    dir.create(output_fldr, recursive = TRUE)
  }
  if (any(!is.na(unique_areas))) {
    for (i in 1:length(unique_areas)) {
      figs_area_out <- paste0(output_fldr, "/id", unique_areas[i])
      if (!dir.exists(figs_area_out)) {
        dir.create(figs_area_out, recursive = TRUE)
      }
    }
  }

  scenario_string <- paste0(cfg$variable, " for ",
                            paste(cfg$period_label, collapse =  " and "), "; baseline ",
                            cfg$baseline_period)
  write(scenario_string, paste0(output_fldr, "/scenario.txt"))
  write(paste0(names(cfg), ";", as.character(cfg)), paste0(output_fldr, "/config.txt"))

  # query files to process
  aggr_files <- list.files(paste0("Data/", cfg$variable),
                           recursive = TRUE, full.names = TRUE)
  aggr_files <- aggr_files[grepl(paste0(ardiv_lbl, "_departure_aggregates"), aggr_files) &
                           grepl(f_out, aggr_files) &
                           grepl(i_out, aggr_files) &
                           grepl(paste(cfg$period_label, collapse = "|"), aggr_files)]

  # constrain ensemble to a narrower set according to an external csv
  if (cfg$ensemble_selection != "") {

    emsel_file <- paste0("data/ensemble_selection/", cfg$ensemble_selection)
    if (file.exists(emsel_file)) {

      selected_ensemble <- read_delim(emsel_file, show_col_types = FALSE, delim = ";") %>%
        filter(variable == cfg$variable) %>%
        mutate(across(-variable, ~ str_replace_all(., "-", "_")),
               ensmem_string = paste0("/", impactmodel, "/", forcing, "/"))

      ems <- paste(selected_ensemble$ensmem_string, collapse = "|")
      aggr_files <- aggr_files[grepl(ems, aggr_files)]

      read_delim(emsel_file, show_col_types = FALSE, delim = ";") %>%
        filter(variable == cfg$variable) %>%
        write.csv(paste0(output_fldr, "/ensemble.csv"))

    } else {
      warning("ensemble selection file missing")
    }
  }

  aggr_data <- aggr_files %>%
    lapply(readRDS) %>%
    bind_rows()

  # unique ensemble members that cover both periods to be compared
  em_cover <- aggr_data %>%
    group_by(impactmodel, forcing, period) %>%
    summarise() %>%
    mutate(dummy = 1) %>%
    summarise(em_cover = sum(dummy)) %>%
    ungroup() %>%
    filter(em_cover == 2)

  nensmem <- nrow(em_cover)
  ensmems_available <- paste0(em_cover$impactmodel, "_", em_cover$forcing)

  aggr_data <- aggr_data %>%
    filter(paste0(impactmodel, "_", forcing) %in% ensmems_available)

  # create output for each area
  for (i in 1:length(unique_areas)) {

    curr_status_data <- tibble()
    if (!is.na(unique_areas[i])) {
      write(paste0("processing area id ", unique_areas[i], "..."),
            cfg$log_file, append = TRUE)
      area_id <- paste0("id", unique_areas[i])
      filter_area <- paste0(cfg$areal_division, "_", area_id)
    } else {
      area_id <- "global"
      filter_area <- "global"
      write("processing global...", cfg$log_file, append = TRUE)
    }

    area_data <- aggr_data %>%
      filter(area == filter_area)
    area_fldr <- ifelse(!is.na(unique_areas[i]), paste0("id", unique_areas[i], "/"), "")

    # check if an area is covered by all or not all ensemble members
    ngroups <- area_data %>%
      group_by(impactmodel, forcing) %>%
      n_groups()

    if (ngroups == 0) {
      write("some small area, all ensemble members missing...",
            cfg$log_file, append = TRUE)
      next
    } else if (ngroups < nensmem) {
      write(paste0(ngroups, " ensemble members covering area..."),
            cfg$log_file, append = TRUE)
    }

    # take annual means of months
    aggregates_year <- area_data %>%
      filter(!is.na(dptClass)) %>%
      group_by(variable, impactmodel, forcing, area, period, year, dptClass) %>%
      summarise(landShareInClass = mean(landShareInClass)) %>%
      ungroup()

    aggregates_year_summed <- aggregates_year %>%
      filter(!is.na(dptClass)) %>%
      group_by(variable, impactmodel, forcing, area, period, year) %>%
      summarise(landShareInClass = sum(landShareInClass)) %>%
      ungroup() %>%
      mutate(dptClass = "summed")

    output_dry_dpts <- .compose_output(cfg, aggregates_year, "below", visopts)
    output_wet_dpts <- .compose_output(cfg, aggregates_year, "above", visopts)
    output_summed_dpts <- .compose_output(cfg, aggregates_year_summed, "summed", visopts)

    arr_annual <- gridExtra::arrangeGrob(grobs = list(output_summed_dpts$plt,
                                                      output_wet_dpts$plt,
                                                      output_dry_dpts$plt))

    plt_year_out <- paste0(output_fldr, "/", area_fldr, cfg$variable, "_",
                           ifelse(!cfg$areal_division == "", paste0(cfg$areal_division, "_"), ""),
                           area_id,
                           "_land_area_with_local_deviations_annual_mean.pdf")
    ggsave(plt_year_out, arr_annual, width = visopts$ts_out_width,
           height = 3 * visopts$ts_out_height, units = visopts$ts_out_units)
    write(paste0("saved ", plt_year_out, "..."), cfg$log_file, append = TRUE)

    export_boundaries <- output_summed_dpts$export_boundaries_current_status %>%
      bind_rows(output_wet_dpts$export_boundaries_current_status) %>%
      bind_rows(output_dry_dpts$export_boundaries_current_status)

    export_lines <- output_summed_dpts$export_ensemble_median_iqr %>%
      bind_rows(output_wet_dpts$export_ensemble_median_iqr) %>%
      bind_rows(output_dry_dpts$export_ensemble_median_iqr)

    export_boundaries_out <- paste0(output_fldr, "/", area_fldr, cfg$variable, "_",
                                    ifelse(!cfg$areal_division == "", paste0(cfg$areal_division, "_"), ""),
                                    area_id,
                                    "_boundaries_current_status_persistent_transgression_year.csv")
    write.table(export_boundaries, export_boundaries_out, sep = ";", row.names = FALSE)
    write(paste0("saved ", export_boundaries_out, "..."), cfg$log_file, append = TRUE)

    export_lines_out <- paste0(output_fldr, "/", area_fldr, cfg$variable, "_",
                               ifelse(!cfg$areal_division == "", paste0(cfg$areal_division, "_"), ""),
                               area_id,
                               "_land_area_with_local_deviations_annual_mean_ensemble_median_IQR.csv")
    write.table(export_lines, export_lines_out, sep = ";", row.names = FALSE)
    write(paste0("saved ", export_lines_out, "..."), cfg$log_file, append = TRUE)

  }

  return (output_fldr)
}
