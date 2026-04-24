# ------------------------------------------------------------------------------
#
# 06_get_area_aggregates
#
# Aggregate local deviations globally or within areas prescribed in config and
# given as arbitrary polygons (or tif raster) in Data/ardiv to get percentage of
# land area with local deviations.
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

.classify_binary_dpt <- function(dpt) {

  # NA if no dpt
  case_when(
    dpt < -1 ~ "below",
    dpt > 1 ~ "above"
  )

}

.get_dpts_from_rast <- function(dpt_binary, ard) {

  if (length(unique(dpt_binary$x)) == 1 |
      length(unique(dpt_binary$y)) == 1 |
      nrow(dpt_binary) == 0) {
    return (NULL)
  }

  dpt_rast_below <- dpt_binary %>%
    select(x, y, starts_with("dpt")) %>%
    mutate(across(starts_with("dpt_"), ~ ifelse(.x == "below", 1, NA))) %>%
    as.matrix() %>%
    rast(type = "xyz", crs = "+proj=lonlat")

  dpt_rast_above <- dpt_binary %>%
    select(x, y, starts_with("dpt")) %>%
    mutate(across(starts_with("dpt_"), ~ ifelse(.x == "above", 1, NA))) %>%
    as.matrix() %>%
    rast(type = "xyz", crs = "+proj=lonlat")

  dpt_rast_area <-  dpt_binary %>%
    select(x, y, cellArea_km2) %>%
    as.matrix() %>%
    rast(type = "xyz", crs = "+proj=lonlat")

  total_area_from_rast <- exact_extract(dpt_rast_area, ard, fun = "sum")
  dpt_below_from_rast <- exact_extract(dpt_rast_below * dpt_rast_area, ard, fun = "sum")
  dpt_above_from_rast <- exact_extract(dpt_rast_above * dpt_rast_area, ard, fun = "sum")

  if (nrow(ard) > 1) { # multipolygons reaching over 180 degrees longitude
    total_area_from_rast <- sum(total_area_from_rast)
    dpt_below_from_rast <- colSums(dpt_below_from_rast)
    dpt_above_from_rast <- colSums(dpt_above_from_rast)
  }

  dpt_areas_below <- tibble(timestep = names(dpt_rast_below),
                            dptClass = "below",
                            areaInClass = as.numeric(dpt_below_from_rast))
  dpt_areas_above <- tibble(timestep = names(dpt_rast_above),
                            dptClass = "above",
                            areaInClass = as.numeric(dpt_above_from_rast))

  ret <- list(dpt_areas_below, dpt_areas_above, total_area_from_rast) %>%
    setNames(c("below", "above", "total_area"))

  return (ret)

}

.aggregate_dpts <- function(cfg, dpts_df, ice_df, ardiv_id, dpts_period,
                            ardiv_df = NULL, ardiv_polys = NULL) {
  # exclusion of ice cells
  if (!is.null(ice_df)) {
    dpts_df <- dpts_df %>%
      left_join(ice_df, by = c("x", "y")) %>%
      filter(is.na(isIce)) %>%
      select(-isIce)
  }

  # filtering of cells within a specific area
  if (!is.null(ardiv_df) | !is.null(ardiv_polys)) {

    area_meta <- paste0(cfg$areal_division, "_id", ardiv_id)
    # ardivs in 1:1 coordinate grid with dpts; (default if both given)
    if (!is.null(ardiv_df)) {
      dpts_df <- dpts_df %>%
        left_join(ardiv_df, by = c("x", "y")) %>%
        filter(id == ardiv_id)
    } else { # ardivs as arbitrary polygons
      ard <- ardiv_polys %>%
        filter(id == ardiv_id)
      bb <- st_bbox(ard)
      dpts_df <- dpts_df %>%
        filter(x >= floor(bb$xmin) & x <= ceiling(bb$xmax) &
               y >= floor(bb$ymin) & y <= ceiling(bb$ymax))
    }
  } else {
    area_meta <- "global"
  }

  # binary departures
  dpt_binary <- dpts_df %>%
    mutate(across(starts_with("dpt_"), ~ .classify_binary_dpt(.x)))

  # arbitrary shape areas if and only if gridded areal divisions not given
  if (!is.null(ardiv_polys) & is.null(ardiv_df)) {

    rast_area_dpt <- .get_dpts_from_rast(dpt_binary, ard)
    if (is.null(rast_area_dpt)) { return (NULL) } # small area not covered by enough raster

    timesteps <- dpts_df %>%
      select(starts_with("dpt")) %>%
      colnames() %>%
      tibble(timestep = .)

    area_dpt <- bind_rows(left_join(timesteps, rast_area_dpt$below, by = "timestep"),
                          left_join(timesteps, rast_area_dpt$above, by = "timestep"))

    if (length(rast_area_dpt) > 3) {
      area_dpt <- area_dpt %>%
        left_join(rast_area_dpt[[4]], by = "timestep") # wshs temp
    }

    area_not_dpt <- area_dpt %>%
      group_by(timestep) %>%
      summarise(areaInClass = rast_area_dpt$total_area - sum(areaInClass)) %>%
      mutate(dptClass = NA)

    area_dpt <- area_dpt %>%
      bind_rows(area_not_dpt) %>%
      mutate(landShareInClass = areaInClass / rast_area_dpt$total_area,
             dptClass = factor(dptClass, levels = c("above", "below", NA))) %>%
      arrange(timestep, dptClass)

  } else { # gridded areal divisions or global

    dpt_binary_long <- dpt_binary %>%
      select(x, y, cellArea_km2, starts_with("dpt_")) %>%
      pivot_longer(-c(x, y, cellArea_km2),
                   names_to = "timestep", values_to = "dptClass")

    # land area shares (%)
    area_dpt <- dpt_binary_long %>%
      group_by(timestep, dptClass) %>%
      summarise(areaInClass = sum(cellArea_km2)) %>%
      mutate(landShareInClass = areaInClass / sum(areaInClass),
             dptClass = factor(dptClass, levels = c("above", "below", NA))) %>%
      ungroup()

  }

  # ensure all dpt classes being present in data
  area_dpt_full <- area_dpt %>%
    expand(timestep, dptClass) %>%
    left_join(area_dpt, by = c("timestep" = "timestep",
                               "dptClass" = "dptClass")) %>%
    mutate(landShareInClass = ifelse(is.na(landShareInClass), 0, landShareInClass))

  # meta
  area_dpt_full <- area_dpt_full %>%
    mutate(variable = cfg$variable,
           timestep = str_replace(timestep, paste0("dpt_", cfg$variable), ""),
           period = dpts_period,
           impactmodel = unique(dpts_df$impactmodel),
           forcing = unique(dpts_df$forcing),
           area = area_meta) %>%
    tidyr::separate(timestep, into = c("year", "month")) %>%
    mutate(year = as.numeric(year),
           month = as.numeric(month)) %>%
    arrange(year) %>%
    select(-areaInClass) %>%
    select(variable, impactmodel, forcing, period, area, everything())

  return (area_dpt_full)

}

# public function
get_area_aggregates <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("computing percentages of land area with local deviations...",
        cfg$log_file, append = TRUE)
  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_")

  # run for a given areal division or global if nothing is given
  ardiv_lbl <- ifelse(cfg$areal_division == "", "global", cfg$areal_division)
  ardiv_df <- NULL
  ardiv_polys <- NULL
  ardiv_ids <- c(NA)

  if (ardiv_lbl != "global") {

    ardiv_rast_exists <- file.exists(paste0("data/ardiv/", cfg$areal_division, ".tif"))
    ardiv_poly_exists <- file.exists(paste0("data/ardiv/", cfg$areal_division, ".gpkg"))

    if (!ardiv_rast_exists & !ardiv_poly_exists) {
      stop(paste0("neither polys nor rast for areal division ",
                  cfg$areal_division, " found"))
    }
    if (ardiv_rast_exists) {
      ardiv_df <- rast(paste0("data/ardiv/", cfg$areal_division, ".tif")) %>%
        as.data.frame(xy = TRUE) %>%
        as_tibble()
      ardiv_ids <- unique(ardiv_df$id)
    }
    if (ardiv_poly_exists) {
      ardiv_polys <- read_sf(paste0("data/ardiv/", cfg$areal_division, ".gpkg"))
      ardiv_ids <- unique(ardiv_polys$id)
    }
  }

  write(paste0("processing areal division ", ardiv_lbl, "..."),
        cfg$log_file, append = TRUE)

  # check that enough directories exist
  skip <- vector(length = length(cfg$period_label))
  for (p in 1:length(cfg$period_label)) {
    dpt_fldr <- paste0("Data/", cfg$variable, "/", cfg$period_label[p],
                       "/", cfg$variable, "_", ardiv_lbl, "_departure_aggregates/")
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
  i_out <- i_out %>% paste(collapse = "|")
  f_out <- f_out %>% paste(collapse = "|")
  dpt_files <- list.files(paste0("Data/", cfg$variable),
                          recursive = TRUE, full.names = TRUE)
  dpt_files <- dpt_files[grepl(paste0(cfg$variable, "_departures"), dpt_files) &
                         grepl(f_out, dpt_files) &
                         grepl(i_out, dpt_files)]

  # compute and save percentages of land area with local deviations
  for (p in 1:length(cfg$period_label)) {

    if (skip[p]) { next() }

    # ice areas to be excluded
    ice_path <- paste0("Data/", cfg$variable, "/", cfg$period_label[p], "/",
                       cfg$variable, "_departures/", cfg$variable,
                       "_ice_cells_", cfg$period_label[p], "_", cfg$grid_label,
                       ".rds")
    if (!file.exists(ice_path)) {
      warning("Ice has not been built for this config")
      ice_df <- NULL
    } else {
      ice_df <- readRDS(ice_path)
    }

    for (m in 1:12) {

      month_files <- dpt_files[grepl(month.name[m], dpt_files) &
                               grepl(cfg$period_label[p], dpt_files)]
      month_data <- month_files %>%
        lapply(FUN = readRDS)

      ardiv_aggregates <- vector(mode = "list", length = length(ardiv_ids))
      for (a in 1:length(ardiv_ids)) {
        month_aggregates <- month_data %>%
          lapply(FUN = .aggregate_dpts, cfg = cfg, ice_df = ice_df,
                 ardiv_id = ardiv_ids[a], dpts_period = cfg$period_label[p],
                 ardiv_df = ardiv_df, ardiv_polys = ardiv_polys)
        ardiv_aggregates[[a]] <- month_aggregates
      }

      out_id <- ifelse(!cfg$areal_division == "", cfg$areal_division, "global")
      out_files <- month_files %>%
        lapply(FUN = function(x, cfg) {
          str_replace_all(x, paste0(cfg$variable, "_departures"),
                          paste0(cfg$variable, "_", out_id, "_departure_aggregates"))},
          cfg = cfg)

      for (i in 1:length(month_files)) {

        all_ardiv_aggregates <- vector(mode = "list", length = length(ardiv_aggregates))
        for (j in 1:length(ardiv_aggregates)) {
          all_ardiv_aggregates[[j]] <- ardiv_aggregates[[j]][[i]]
        }

        # dimensions should be 3 classes x n years per period
        data_coverage <- 3 * (cfg$period_end[p] - cfg$period_begin[p] + 1)
        check_coverage <- all_ardiv_aggregates %>%
          lapply(nrow) %>%
          unlist()
        not_covered <- which(check_coverage != data_coverage)

        # log areas which are not covered by departures data
        if (length(not_covered) > 0) {
          write(paste0("data not covering area id ", ardiv_ids[not_covered]),
                cfg$log_file, append = TRUE)
        }
        all_ardiv_aggregates <- all_ardiv_aggregates[check_coverage == data_coverage] %>%
          bind_rows()

        saveRDS(all_ardiv_aggregates, out_files[[i]])
        write(paste0("saved ", out_files[[i]]), cfg$log_file, append = TRUE)

      }
    }
  }
}
