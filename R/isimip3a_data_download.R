# ------------------------------------------------------------------------------
#
# isimip3a_data_download
#
# Download data from ISIMIP data portal using a list of URLs, with variable to
# be downloaded and download directory specified from cmd; impact models and
# forcing scenarios prescribed here as constant.
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

library(curl)
library(doParallel)
library(stringr)
library(dplyr)
library(tidyr)

# ISIMIP 3a bulk download URLs -------------------------------------------------

cmd_args <- commandArgs(trailingOnly = TRUE)
dvar <- cmd_args[1]
dl_loc <- cmd_args[2]

dl_dir <- paste0(dl_loc, "/", dvar)
if (!dir.exists(dl_dir)) { dir.create(dl_dir) }

ghm <- c("CLASSIC", "CWatM", "H08", "HydroPy", "JULES-ES-VN6P3", "JULES-W2",
         "LPJmL5-7-10-fire", "MIROC-INTEG-LAND", "ORCHIDEE-MICT",
         "SSiB4-TRIFFID-Fire", "VISIT", "WaterGAP2-2e",
         "WEB-DHM-SG") # + "JULES-W2-DDM30"

clim <- c("20CRv3", "20CRv3-ERA5", "20CRv3-W5E5", "GSWP3-W5E5")

exp <- c("obsclim_histsoc_default",
         "counterclim_histsoc_default",
         "obsclim_1901soc_default",
         "counterclim_1901soc_default")

for (i in 1:length(exp)) {
  exp_dir <- paste0(dl_dir, "/", exp[i], "/raw")
  if (!dir.exists(exp_dir)) { dir.create(exp_dir, recursive = TRUE) }
}

urls <- c()
for (i in 1:length(ghm)) {
  for (j in 1:length(clim)) {
    for (k in 1:length(exp)) {

      ending <- case_when(
        clim[j] == "20CRv3" ~ "1901_2015.nc",
        clim[j] == "20CRv3-ERA5" ~ "1901_2021.nc",
        clim[j] == "20CRv3-W5E5" ~ "1901_2019.nc",
        clim[j] == "GSWP3-W5E5" ~ "1901_2019.nc"
      )
      newurl <- paste0("https://files.isimip.org/ISIMIP3a/OutputData/water_global/",
                       ghm[i], "/", tolower(clim[j]), "/historical/", tolower(ghm[i]),
                       "_", tolower(clim[j]), "_", exp[k], "_", dvar, "_", ending)
      urls <- c(urls, newurl)

    }
  }
}

# download ---------------------------------------------------------------------

for (i in 1:length(urls)) {

  out <- urls[i] %>%
    str_split("/") %>%
    unlist() %>%
    tail(1)

  out_loc <- paste0(dl_dir, "/",
                    exp[which(lapply(exp, grepl, x = out) %>% unlist())], "/raw/",
                    out)
  h <- new_handle()
  handle_setopt(h, ssl_verifyhost = 0, ssl_verifypeer = 0)

  if (!file.exists(out_loc) & RCurl::url.exists(urls[i])){

    message(paste0("downloading: ", out))
    curl_download(urls[i], out_loc, handle = h)

  } else if (file.exists(out_loc)) {

    message(paste0("exists already: ", out))

  } else if (!RCurl::url.exists(urls[i])) {

    message(paste0("not in repository: ", out))

  }

}

# ISIMIP 3a mapping of downloaded data files -----------------------------------

vars <- c("dis", "rootmoist")
exp_crf <- c("obsclim", "counterclim")
exp_dhf <- c("histsoc", "1901soc")
exp_sens <- c("default")
ghm_lower <- tolower(ghm)
clim_lower <- tolower(clim)

ll <- list.files(dl_dir, full.names = TRUE, recursive = TRUE)
meta_list <- vector(mode = "list", length = length(ll))
for (i in 1:length(ll)) {

  data_file <- ll[i] %>%
    str_split("/") %>%
    unlist() %>%
    tail(1) %>%
    str_split("_") %>%
    unlist()

  meta <- c(vars[which(vars %in% data_file)],
            ghm_lower[which(ghm_lower %in% data_file)],
            clim_lower[which(clim_lower %in% data_file)],
            exp_crf[which(exp_crf %in% data_file)],
            exp_dhf[which(exp_dhf %in% data_file)],
            exp_sens[which(exp_sens %in% data_file)]) %>%
    setNames(c("variable", "impactmodel", "forcing", "exp_crf", "exp_dhf", "exp_sens"))
  meta_list[[i]] <- meta

}

meta_tbl <- meta_list %>%
  bind_rows() %>%
  mutate(experiment = paste(exp_crf, exp_dhf, exp_sens, sep = "_")) %>%
  select(-starts_with("exp_")) %>%
  arrange(variable, experiment, impactmodel, forcing) %>%
  mutate(exists = 1)

meta_exhaustive <- expand_grid(variable = vars,
                               impactmodel = ghm_lower,
                               forcing = clim_lower,
                               experiment = exp) %>%
  left_join(meta_tbl, by = c("variable", "impactmodel", "forcing", "experiment")) %>%
  pivot_wider(id_cols = c(variable, impactmodel, forcing),
              names_from = experiment,
              values_from = exists) %>%
  filter(!if_all(-c(variable, impactmodel, forcing), ~ is.na(.x)))

write.table(meta_exhaustive, paste0(dl_dir, "/isimip3a_raw_table_", Sys.time(), ".csv"),
            sep = ";", row.names = FALSE)

# rename folder
dummy <- file.rename(dl_dir, str_replace(dl_dir, "_global_monthly", ""))
