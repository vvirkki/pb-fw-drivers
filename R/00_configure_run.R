# ------------------------------------------------------------------------------
#
# 00_configure_run
#
# Read a configuration file describing a run setup of GHMs and GCMs together
# with parameters on how the analysis should be run.
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

configure_run <- function(config_path) {

  # read config file; return a list of options to be used in analysis
  cfg_data <- readLines(config_path) %>%
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

  cfg$tempdir <- tempdir()
  if (!dir.exists("logs")) { dir.create("logs") }
  cfg$log_file <- paste0("logs/log ", Sys.time(), " ", cfg$variable, " ",
                         paste0(paste(cfg$impactmodel, collapse = "_"), " ",
                                paste(cfg$forcing, collapse = "_"), ".txt"))

  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("opening log...", cfg$log_file, append = TRUE)

  return (cfg)

}
