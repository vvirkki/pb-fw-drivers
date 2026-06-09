# ------------------------------------------------------------------------------
#
# 01_prepare_monthly_data
#
# Read ISIMIP NetCDF rasters and transform them to R data frames describing
# monthly values (columns) in grid cells (rows).
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

options(dplyr.summarise.inform = FALSE)

.check_outdir <- function(out_path) {
  is_filled <- lapply(month.name, grepl, x = list.files(out_path)) %>%
    lapply(any) %>%
    unlist() %>%
    all()
  return (is_filled)
}

# fill value is mean of all non-NA values
.fill_NA_values <- function(x) {
  xvals <- x[3:length(x)]
  xvals <- ifelse(is.na(xvals), mean(xvals, na.rm = TRUE), xvals)
  ret <- c(x[1], x[2], xvals) %>%
    setNames(names(x))
  return (ret)
}

.detect_cpt <- function(x) {
  xvals <- x[3:length(x)]
  if (var(xvals) < 1e-6) {
    cp <- length(xvals) # minor variance; changepoint set to the end of data
  } else {
    cp <- cpt.meanvar(xvals, method = "AMOC", class = FALSE)["cpt"]
  }
  ret <- c(x[1], x[2], cp) %>%
    setNames(c("x", "y", "cpt"))
  return (ret)
}

# ISIMIP netcdfs -> monthly full-timeseries cell data frames -------------------
.prepare_df <- function(cfg, period, impactmodel, forcing) {

  # setup
  period_no <- which(cfg$period_label == period)
  ystart <- cfg$period_begin[period_no]
  yend <- cfg$period_end[period_no]

  raw_files <- paste0(cfg$raw_data_root, "/", cfg$variable, "/", period, "/raw") %>%
    list.files(full.names = TRUE)
  files_in <- raw_files[grepl(impactmodel, raw_files)]
  files_in <- files_in[grepl(forcing, files_in)]

  if (length(files_in) == 0) {

    ensmem_not_exist <- paste0("ensemble member ", impactmodel, " ", forcing, " ", period, " missing")
    write(ensmem_not_exist, cfg$log_file, append = TRUE)
    stop(ensmem_not_exist)

  }

  if (all(grepl("daily", files_in))) { # daily data

    data_chunks <- files_in %>%
      lapply(., .daily_to_monthly, nm_prefix = cfg$variable)
    data_rast <- rast(data_chunks)

  } else { # monthly data

    data_rast <- rast(files_in)

    if (any(grepl("depth", varnames(data_rast)))) {
      stop("depth layer found")
      data_rast <- data_rast[[-1]]
    }

    # temporal extent determined by config params
    nm <- paste0(cfg$variable,
                 sort(rep(seq(ystart, yend, 1), 12)), ".",
                 str_pad(seq(1,12,1), 2, side = "left", pad = "0"))

    # for LPJmL rootmoist, take sum of all depth layers (kg/m3)
    if (impactmodel == "lpjml5-7-10-fire" & cfg$variable == "rootmoist") {

      nminit <- names(data_rast)
      ts <- nminit %>%
        str_split("_") %>%
        lapply(function(x) { return(x[[3]]) }) %>%
        unlist() %>%
        unique() %>%
        as.numeric()

      sum_rast <- vector(mode = "list", length = length(ts))
      for (i in 1:length(ts)) {
        sel <- paste0("rootmoist_depth=", c(0.1, 0.35, 0.75, 1.5, 2.5), "_", ts[i])
        sum_rast[[i]] <- sum(data_rast[[sel]])
      }
      data_rast <- rast(sum_rast)

    }

    if (forcing == "20crv3-era5") {
      data_rast <- data_rast[[1:length(nm)]] # limit years to 2019 (2021 in data)
    }
    names(data_rast) <- nm

  }

  cell_areas <- data_rast[[1]] %>%
    cellSize(unit = "km") %>%
    setNames(c("cellArea_km2")) %>%
    as.data.frame(xy = TRUE)

  impactmodel_out <- str_replace_all(impactmodel, "-", "_")
  forcing_out <- str_replace_all(forcing, "-", "_")
  mids <- str_pad(seq(1,12,1), 2, side = "left", pad = "0")
  data_nm <- names(data_rast)

  # prepare output
  out_fldr <- paste0("Data/", cfg$variable, "/", period, "/", cfg$variable,
                     "_monthly",
                     ifelse(impactmodel != "",
                            paste0("/", impactmodel_out), ""),
                     paste0("/", forcing_out))

  if (!.check_outdir(out_fldr)) {
    dir.create(out_fldr, recursive = TRUE)
  } else {
    write(paste0("outputs in ", out_fldr, " are complete, skipping..."),
          cfg$log_file, append = TRUE)
    return ()
  }

  for (m in 1:length(mids)) {

    sel <- paste0(".", mids[m])
    sel_layers <- data_nm[grepl(sel, data_nm, fixed = TRUE)]

    # rast to data frame
    mnth_df <- data_rast %>%
      subset(sel_layers) %>%
      as.data.frame(xy = TRUE, cells = TRUE) %>%
      mutate(impactmodel = impactmodel_out,
             forcing = forcing_out) %>%
      left_join(cell_areas, by = c("x" = "x", "y" = "y"))

    # negative and NA diagnostics
    infilling_diagnostics <- mnth_df %>%
      filter(if_any(starts_with(cfg$variable), ~ is.na(.x) | .x < 0)) %>%
      mutate(month = month.name[m]) %>%
      select(impactmodel, forcing, month, cell, cellArea_km2, x, y, everything())

    # convert negative values to NA
    mnth_df <- mnth_df %>%
      mutate(across(starts_with(cfg$variable), ~ ifelse(.x < 0, NA, .x)))

    # check for NA values, fill with mean of non-NA values and write to log if any
    na_values <- mnth_df %>%
      filter(if_any(starts_with(cfg$variable), ~ is.na(.x)))

    if (nrow(na_values) > 0) {

      nvals_mat <- na_values %>%
        select(starts_with(cfg$variable)) %>%
        as.matrix()

      na_filled_values <- mnth_df %>%
        filter(cell %in% na_values$cell) %>%
        select(x, y, starts_with(cfg$variable)) %>%
        as.matrix() %>%
        apply(MARGIN = 1, FUN = .fill_NA_values, simplify = FALSE) %>%
        bind_rows()

      non_filled_values <- mnth_df %>%
        filter(!cell %in% na_values$cell) %>%
        select(x, y, starts_with(cfg$variable))

      mnth_df <- mnth_df %>%
        select(-starts_with(cfg$variable)) %>%
        left_join(bind_rows(na_filled_values, non_filled_values),
                  by = c("x" = "x", "y" = "y"))

      write(paste0("negative or NA monthly ", cfg$variable, " in ",
                   nrow(na_values), " cells in ",
                   length(which(is.na(nvals_mat))), " values"),
            cfg$log_file, append = TRUE)

    }

    # detect changepoints
    cpts <- mnth_df %>%
      select(x, y, starts_with(cfg$variable)) %>%
      as.matrix() %>%
      apply(MARGIN = 1, FUN = .detect_cpt, simplify = FALSE) %>%
      bind_rows()

    mnth_df <- mnth_df %>%
      left_join(cpts, by = c("x" = "x", "y" = "y")) %>%
      select(impactmodel, forcing, cell, cellArea_km2, x, y, cpt, everything())

    nrow_neg <- nrow(mnth_df %>% filter(if_any(starts_with(cfg$variable), ~ . < 0)))
    if (nrow_neg > 0) {
      stop("negative values still found in data preparation")
    }
    out_file <- paste0(cfg$variable, "_monthly_", period, "_", month.name[m], "_",
                       ifelse(!is.na(impactmodel_out),
                              paste0(impactmodel_out, "_"), ""),
                       forcing_out, ".rds")
    out <- paste0(out_fldr, "/", out_file)
    saveRDS(mnth_df, out)
    write(paste0("saved ", out), cfg$log_file, append = TRUE)

    if (nrow(infilling_diagnostics) > 0) {
      saveRDS(infilling_diagnostics,
              paste0("logs/", str_replace(out_file, "_monthly_", "_infilled_")))
    }

  }
}

# public function
prepare_monthly_data <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("preparing data...", cfg$log_file, append = TRUE)
  periods <- cfg$period_label
  impactmodels <- cfg$impactmodel
  forcings <- cfg$forcing

  for (i in 1:length(periods)) {
    for (j in 1:length(impactmodels)) {
      for (k in 1:length(forcings)) {
        .prepare_df(cfg, periods[i], impactmodels[j], forcings[k])
      }
    }
  }
}
