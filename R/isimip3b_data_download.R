# ------------------------------------------------------------------------------
#
# isimip3b_data_download
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

# ISIMIP 3b bulk download URLs -------------------------------------------------

cmd_args <- commandArgs(trailingOnly = TRUE)
dvar <- cmd_args[1]
dl_loc <- cmd_args[2]

dl_dir <- paste0(dl_loc, "/", dvar)
if (!dir.exists(dl_dir)) { dir.create(dl_dir) }

ghm <- c("CLASSIC", "CWatM", "H08",  "JULES-ES-VN6P3", "JULES-W2",
         "MIROC-INTEG-LAND", "VISIT", "WaterGAP2-2e", "WEB-DHM-SG")

gcm <- c("GFDL-ESM4", "UKESM1-0-LL", "MPI-ESM1-2-HR", "IPSL-CM6A-LR", "MRI-ESM2-0")

exp <- c("historical_histsoc_default",
         "picontrol_1850soc_default",
         "ssp126_2015soc_default",
         "ssp370_2015soc_default",
         "ssp585_2015soc_default")

for (i in 1:length(exp)) {
  exp_dir <- paste0(dl_dir, "/", exp[i], "/raw")
  if (!dir.exists(exp_dir)) { dir.create(exp_dir, recursive = TRUE) }
}

urls <- c()
for (i in 1:length(ghm)) {
  for (j in 1:length(gcm)) {
    for (k in 1:length(exp)) {

      period <- case_when(
        grepl("ssp", exp[k]) ~ "future",
        TRUE ~ "historical"
      )
      ending <- case_when(
        grepl("ssp", exp[k]) ~ "2015_2100.nc",
        grepl("historical", exp[k]) ~ "1850_2014.nc",
        grepl("picontrol", exp[k]) ~ "1850_2014.nc"
      )

      newurl <- paste0("https://files.isimip.org/ISIMIP3b/OutputData/water_global/",
                       ghm[i], "/", tolower(gcm[j]), "/", period, "/", tolower(ghm[i]),
                       "_", tolower(gcm[j]), "_w5e5_", exp[k], "_", dvar, "_", ending)
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

# ISIMIP 3b mapping of downloaded data files -----------------------------------

vars <- c("dis", "rootmoist")
exp_crf <- c("historical", "picontrol", "ssp126", "ssp370", "ssp585")
exp_dhf <- c("histsoc", "1850soc", "2015soc")
exp_sens <- c("default")
ghm_lower <- tolower(ghm)
gcm_lower <- tolower(gcm)

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
            gcm_lower[which(gcm_lower %in% data_file)],
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
                               forcing = gcm_lower,
                               experiment = exp) %>%
  left_join(meta_tbl, by = c("variable", "impactmodel", "forcing", "experiment")) %>%
  pivot_wider(id_cols = c(variable, impactmodel, forcing),
              names_from = experiment,
              values_from = exists) %>%
  filter(!if_all(-c(variable, impactmodel, forcing), ~ is.na(.x)))

write.table(meta_exhaustive, paste0(dl_dir, "/isimip3b_raw_table_",
                                    Sys.time(), ".csv"), sep = ";", row.names = FALSE)

# rename folder
dummy <- file.rename(dl_dir, str_replace(dl_dir, "_global_monthly", ""))
