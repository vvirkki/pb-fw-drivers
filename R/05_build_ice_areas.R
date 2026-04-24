# ------------------------------------------------------------------------------
#
# 05_build_ice_areas
#
# Query HYDE anthromes within each period described in config and parse a binary
# grid cell data frame describing whether a cell is permanent ice or not.
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

.count_ice_years <- function(x) {

  ice_cols <- x[3:length(x)]
  return (c(x[1], x[2], length(ice_cols[ice_cols > 0])) %>%
            setNames(c("x", "y", "isIce")))

}

.build_ice <- function(anthrome_in, template_raster) {

  ice_year <- readr::parse_number(anthrome_in)
  anthromes <- rast(anthrome_in)
  ice <- anthromes == 63
  ice[is.na(ice)] <- 0
  aggregate_fact <- (res(template_raster) / res(ice)) %>%
    mean() %>%
    round()

  ice_df <- ice %>%
    terra::aggregate(fact = aggregate_fact, fun = "modal") %>%
    terra::resample(template_raster, method = "near") %>%
    as.data.frame(xy = TRUE) %>%
    rename(ice := 3) %>%
    mutate(ice = as.logical(ice),
           year = readr::parse_number(anthrome_in)) %>%
    as_tibble()

  return (ice_df)

}

# public function
build_ice_areas <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("building ice areas...", cfg$log_file, append = TRUE)
  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_") %>%
    paste(collapse = "|")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_") %>%
    paste(collapse = "|")

  # query files to process
  dpt_files <- list.files(paste0("Data/", cfg$variable),
                          recursive = TRUE, full.names = TRUE)
  dpt_files <- dpt_files[grepl(paste0(cfg$variable, "_departures"), dpt_files) &
                           grepl(f_out, dpt_files) &
                           grepl(i_out, dpt_files)]

  ice_template <- readRDS(dpt_files[[1]]) %>%
    select(x, y) %>%
    mutate(value = NA) %>%
    as.matrix() %>%
    rast(type = "xyz", crs = "+proj=lonlat")

  # build ice areas
  for (p in 1:length(cfg$period_label)) {

    # check if ice for this grid & period has already been built
    out_fldr <- paste0("Data/", cfg$variable, "/", cfg$period_label[p], "/",
                       cfg$variable, "_departures/")
    out_file <- paste0(cfg$variable, "_ice_cells_", cfg$period_label[p], "_",
                       cfg$grid_label, ".rds")
    out <- paste0(out_fldr, out_file)
    if (file.exists(out)) {
      write(paste0("ice built for this period ", cfg$period_label[p],
                   " already exists..."), cfg$log_file, append = TRUE)
      next()
    }

    period_begin <- cfg$period_begin[which(cfg$period_label == cfg$period_label[p])]
    period_end <- cfg$period_end[which(cfg$period_label == cfg$period_label[p])]

    # REQUIRED: HYDE 3.5 Anthromes
    # https://doi.org/10.24416/UU01-F45D44
    anthromes <- list.files("Data/anthromes", full.names = TRUE)

    data_years <- seq(period_begin, period_end, 1) %>%
      paste(collapse = "|")
    ice_files <- anthromes[grepl(data_years, anthromes) & grepl("AD", anthromes)]
    ice_years <- readr::parse_number(ice_files)
    ice_files <- ice_files[which(ice_years >= period_begin & ice_years <= period_end)]

    ice_df <- ice_files %>%
      lapply(.build_ice, template_raster = ice_template) %>%
      bind_rows() %>%
      pivot_wider(id_cols = c(x, y), values_from = ice,
                  names_from = year, names_prefix = "ice_")

    ice_majority_vote <- ice_df %>%
      as.matrix() %>%
      apply(MARGIN = 1, FUN = .count_ice_years, simplify = FALSE) %>%
      bind_rows() %>%
      as_tibble() %>%
      filter(isIce >= (ncol(ice_df) - 2) / 2) %>%
      mutate(isIce = TRUE)

    saveRDS(ice_majority_vote, out)
    write(paste0("saved ", out), cfg$log_file, append = TRUE)

  }
}
