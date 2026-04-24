# ------------------------------------------------------------------------------
#
# 04_detect_local_deviations
#
# Detect local deviations (here also termed departures) in all grid cells, save
# deviations as differences to the local baseline range (dry deviations; values
# are < -1, wet deviations; values are > 1).
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

options(dplyr.summarise.inform = FALSE)

.check_departure <- function(obs, low, up, mid) {

  if (any(obs < 0)) { stop("negative obs found (should not be)") }

  d <- tibble(obs, low, mid, up) %>%
    mutate(position = case_when(
      up == 0 ~ 0,
      up == low ~ 0,
      up == mid ~ 0,
      low == mid ~ 0,
      obs > up ~ 1 + (obs - up) / (up - mid),         # above upper bound
      obs < low ~ -1 + (low - obs) / (low - mid),     # below lower bound
      obs >= mid ~ (obs - mid) / (up - mid),          # btwn median and upper bound
      obs < mid ~ (obs - mid) / (mid - low)           # btwn lower bound and median
    ))

  if (nrow(d %>% filter(is.infinite(position) | is.na(position))) > 0) {
    stop("infinite/NA values found in normalised local departures")
  }

  return (d$position)

}

.compare <- function(cfg, file_local_bounds, compare_to) {

  local_bounds <- readRDS(file_local_bounds)
  comparison_data <- file_local_bounds %>%
    str_replace_all("local_bounds", "monthly") %>%
    str_replace_all(cfg$baseline_period, compare_to) %>%
    readRDS() %>%
    as_tibble()

  data <- local_bounds %>%
    left_join(comparison_data, by = c("impactmodel", "forcing", "cell",
                                      "cellArea_km2", "x", "y")) %>%
    as_tibble()

  departures <- data %>%
    mutate(across(starts_with(cfg$variable),
                  ~ .check_departure(.x, local_bound_low, local_bound_high,
                                     local_bound_median))) %>%
    rename_with(~ paste0("dpt_", .x), starts_with(cfg$variable))

  missing_data_cells <- departures %>%
    filter(is.na(cpt)) %>%
    pull(cell)

  if (length(missing_data_cells) > 0) {
    write(paste0(length(missing_data_cells), " cells eliminated for no comparison data"),
          cfg$log_file, append = TRUE)
    departures <- departures %>%
      filter(!cell %in% missing_data_cells)
  }

  return (departures)

}

# public function
detect_local_deviations <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("detecting local deviations...", cfg$log_file, append = TRUE)
  baseline_period <- cfg$baseline_period
  comparison_period <- cfg$period_label[which(cfg$period_label != baseline_period)]
  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_")

  # check that enough directories exist
  skip <- vector(length = length(cfg$period_label))
  for (p in 1:length(cfg$period_label)) {
    dpt_fldr <- paste0("Data/", cfg$variable, "/", cfg$period_label[p],
                       "/", cfg$variable, "_departures/")
    for (i in 1:length(i_out)) {
      for (j in 1:length(f_out)) {
        outdir <- paste0(dpt_fldr,
                         ifelse(i_out[i] == "", "", paste0(i_out[i], "/")),
                         f_out[j])
        if (!.check_outdir(outdir)) {
          dir.create(outdir, recursive = TRUE)
        } else {
          write(paste0("outputs in ", outdir, " are complete, skipping..."),
                cfg$log_file, append = TRUE)
          skip[p] <- TRUE
        }
      }
    }
  }

  # query files to process
  local_bounds <- paste0("Data/", cfg$variable, "/", baseline_period, "/",
                         cfg$variable, "_local_bounds") %>%
    list.files(recursive = TRUE, full.names = TRUE)
  i_out <- i_out %>% paste(collapse = "|")
  f_out <- f_out %>% paste(collapse = "|")
  local_bounds <- local_bounds[grepl(i_out, local_bounds) &
                               grepl(f_out, local_bounds)]

  # compute and save departures
  if (!skip[which(cfg$period_label == baseline_period)]) {

    for (i in 1:length(local_bounds)) {

      baseline_departures <- .compare(cfg, local_bounds[i], baseline_period)
      baseline_departures_out <- local_bounds[i] %>%
        str_replace_all("local_bounds", "departures")
      saveRDS(baseline_departures, baseline_departures_out)
      write(paste0("saved ", baseline_departures_out), cfg$log_file, append = TRUE)

    }
  }

  if (!skip[which(cfg$period_label == comparison_period)]) {

    for (i in 1:length(local_bounds)) {

      comparison_departures <- .compare(cfg, local_bounds[i], comparison_period)
      comparison_departures_out <- local_bounds[i] %>%
        str_replace_all("local_bounds", "departures") %>%
        str_replace_all(baseline_period, comparison_period)
      saveRDS(comparison_departures, comparison_departures_out)
      write(paste0("saved ", comparison_departures_out), cfg$log_file, append = TRUE)

    }
  }
}
