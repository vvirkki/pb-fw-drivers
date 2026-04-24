# ------------------------------------------------------------------------------
#
# 03_set_local_baseline_range
#
# Set the local baseline range (local variability bounds) using baseline period
# quantiles given in config parameters local_bound_low and local_bound_high.
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

.get_quantile <- function(x, quant, cptmin, cptmax) {

  if (cptmin == 0 & cptmax == 0) {
    # first data values indicating spinup already excluded
    return (c(x[1], x[2], quantile(x[4:length(x)], quant)))
  } else {
    stop ("cpt based exclusion should not be done at all")
    cpt <- x[3]
    cpt <- ifelse(cpt < cptmin | cpt > cptmax, cptmin, cpt)
    # x,y,cpt,quantile excluding cpt first data values
    return (c(x[1], x[2], quantile(x[(cpt+4):length(x)], quant)))
  }

}

.determine_bounds <- function(file_in, cfg) {

  out <- file_in %>%
    str_replace_all(paste0("_monthly"), "_local_bounds")

  sel_years <- paste0(cfg$variable, seq(cfg$local_bound_begin, cfg$local_bound_end, 1)) %>%
    paste(collapse = "|")

  data <- readRDS(file_in)

  lblow <- data %>%
    select(x, y, cpt, matches(sel_years)) %>%
    as.matrix() %>%
    apply(MARGIN = 1, FUN = .get_quantile, quant = cfg$local_bound_low,
          cptmin = cfg$cptmin, cptmax = cfg$cptmax, simplify = FALSE) %>%
    bind_rows() %>%
    setNames(c("x", "y", "local_bound_low"))

  lbhigh <- data %>%
    select(x, y, cpt, matches(sel_years)) %>%
    as.matrix() %>%
    apply(MARGIN = 1, FUN = .get_quantile, quant = cfg$local_bound_high,
          cptmin = cfg$cptmin, cptmax = cfg$cptmax, simplify = FALSE) %>%
    bind_rows() %>%
    setNames(c("x", "y", "local_bound_high"))

  lbmedian <- data %>%
    select(x, y, cpt, matches(sel_years)) %>%
    as.matrix() %>%
    apply(MARGIN = 1, FUN = .get_quantile, quant = 0.5,
          cptmin = cfg$cptmin, cptmax = cfg$cptmax, simplify = FALSE) %>%
    bind_rows() %>%
    setNames(c("x", "y", "local_bound_median"))

  local_bounds <- data %>%
    select(-c(starts_with(cfg$variable), cpt)) %>%
    left_join(lblow, by = c("x", "y")) %>%
    left_join(lbhigh, by = c("x", "y")) %>%
    left_join(lbmedian, by = c("x", "y"))

  saveRDS(local_bounds, out)
  write(paste0("saved ", out), cfg$log_file, append = TRUE)

}

# public function
set_local_baseline_range <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("creating dry and wet local bounds...", cfg$log_file, append = TRUE)
  vmon <- paste0(cfg$variable, "_monthly")

  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_")

  # check that enough directories exist
  local_bounds_fldr <- paste0("Data/", cfg$variable, "/", cfg$baseline_period,
                              "/", cfg$variable, "_local_bounds/")
  for (i in 1:length(i_out)) {
    for (j in 1:length(f_out)) {
      outdir <- paste0(local_bounds_fldr,
                       ifelse(i_out[i] == "", "", paste0(i_out[i], "/")),
                       f_out[j])
      if (!.check_outdir(outdir)) {
        dir.create(outdir, recursive = TRUE)
      } else {
        write(paste0("outputs in ", outdir, " are complete, skipping..."),
              cfg$log_file, append = TRUE)
        return ()
      }
    }
  }

  # query files to process
  files_in <- paste0("Data/", cfg$variable, "/",
                     cfg$baseline_period, "/", vmon) %>%
    list.files(full.names = TRUE, recursive = TRUE)
  i_out <- i_out %>% paste(collapse = "|")
  f_out <- f_out %>% paste(collapse = "|")
  files_in <- files_in[grepl(i_out, files_in) & grepl(f_out, files_in)]

  # compute and save local bounds
  for (i in 1:length(files_in)) {
    .determine_bounds(files_in[i], cfg)
  }
}
