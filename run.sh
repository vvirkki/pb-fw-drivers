#!/bin/bash

# ------------------------------------------------------------------------------
### blue water analysis begins -------------------------------------------------
# ------------------------------------------------------------------------------

# prepare data and intermediate outputs ----------------------------------------

Rscript R/main_dataprep.R dis 12
Rscript R/main_outputs.R dis

Rscript R/supplementary_outputs.R dis

Rscript R/future_dataprep.R dis 12
Rscript R/future_outputs.R dis

# prepare tables and raw graphics for figures besides Fig. S3-S4 ---------------

Rscript R/plot01_global_timeseries.R dis main_ensemble figures_tables_timeseries_main
Rscript R/plot01_global_timeseries.R dis supplementary_ensemble figures_tables_timeseries_figS12

Rscript R/plot02_global_maps.R dis hybas4 main_ensemble figures_tables_maps_main

# preparing tables and raw graphics ends ---------------------------------------

# move supplementary outputs around to compose Fig. S3 and Fig. S4 -------------

from_dir="output/dis_supplementary_single_ensemble_members"
to_dir="output/dis_supplementary_ensemble"
plot_dir="figures_tables_maps_figS3S4"
ensmems=(
  "watergap2-2e_20crv3-era5"
  "watergap2-2e_20crv3-w5e5"
  "watergap2-2e_gswp3-w5e5"
)

for em in "${ensmems[@]}"; do

  # pick up individual ensemble member outputs one by one
  created_em_dirs=()
  while IFS= read -r src_em_dir; do

      scenario_dir=$(basename "$(dirname "$src_em_dir")")
      dest_em_dir="$to_dir/$scenario_dir/$em"
      mkdir -p "$dest_em_dir"
      cp -a -- "$src_em_dir"/. "$dest_em_dir"/
      created_em_dirs+=("$dest_em_dir")

  done < <(find "$from_dir" -mindepth 2 -maxdepth 2 -type d -name "$em")

  # create and collate maps for an individual ensemble member
  Rscript R/plot02_global_maps.R dis hybas4 supplementary_ensemble figures_tables_maps_figS3S4

  em_plot_dir="$plot_dir/${em}"
  mkdir -p "$em_plot_dir"
  find "$plot_dir" -maxdepth 1 -type f -exec mv -- {} "$em_plot_dir"/ \;

  # clean up before moving to the next ensemble member
  for d in "${created_em_dirs[@]}"; do
      rm -rf -- "$d"
  done

done

# composing Fig. S3-S4 ends ----------------------------------------------------
# blue water analysis ends -----------------------------------------------------

# ------------------------------------------------------------------------------
### green water analysis begins ------------------------------------------------
# ------------------------------------------------------------------------------

# prepare data and intermediate outputs ----------------------------------------

Rscript R/main_dataprep.R rootmoist 8
Rscript R/main_outputs.R rootmoist

Rscript R/supplementary_outputs.R rootmoist

Rscript R/future_dataprep.R rootmoist 8
Rscript R/future_outputs.R rootmoist

# prepare tables and raw graphics for figures besides Fig. S3-S4 ---------------

Rscript R/plot01_global_timeseries.R rootmoist main_ensemble figures_tables_timeseries_main
Rscript R/plot01_global_timeseries.R rootmoist supplementary_ensemble figures_tables_timeseries_figS12

Rscript R/plot02_global_maps.R rootmoist hybas4 main_ensemble figures_tables_maps_main

# preparing tables and raw graphics ends ---------------------------------------

# move supplementary outputs around to compose Fig. S3 and Fig. S4 -------------

from_dir="output/rootmoist_supplementary_single_ensemble_members"
to_dir="output/rootmoist_supplementary_ensemble"
plot_dir="figures_tables_maps_figS3S4"
ensmems=(
  "miroc-integ-land_20crv3-era5"
  "miroc-integ-land_20crv3-w5e5"
)

for em in "${ensmems[@]}"; do

  # pick up individual ensemble member outputs one by one
  created_em_dirs=()
  while IFS= read -r src_em_dir; do

      scenario_dir=$(basename "$(dirname "$src_em_dir")")
      dest_em_dir="$to_dir/$scenario_dir/$em"
      mkdir -p "$dest_em_dir"
      cp -a -- "$src_em_dir"/. "$dest_em_dir"/
      created_em_dirs+=("$dest_em_dir")

  done < <(find "$from_dir" -mindepth 2 -maxdepth 2 -type d -name "$em")

  # create and collate maps for an individual ensemble member
  Rscript R/plot02_global_maps.R rootmoist hybas4 supplementary_ensemble figures_tables_maps_figS3S4

  em_plot_dir="$plot_dir/${em}"
  mkdir -p "$em_plot_dir"
  find "$plot_dir" -maxdepth 1 -type f -exec mv -- {} "$em_plot_dir"/ \;

  # clean up before moving to the next ensemble member
  for d in "${created_em_dirs[@]}"; do
      rm -rf -- "$d"
  done

done

# composing Fig. S3-S4 ends ----------------------------------------------------
# green water analysis ends ----------------------------------------------------

# ------------------------------------------------------------------------------
### joint blue and green water analysis begins ---------------------------------
# ------------------------------------------------------------------------------

Rscript R/plot03_percentile_rank_figure.R hybas4 figures_tables_maps_main figures_tables_maps_main

### joint blue and green water analysis ends -----------------------------------

# ------------------------------------------------------------------------------
### final data exports for Planetary Health Check ------------------------------
# ------------------------------------------------------------------------------

Rscript R/temp01_planetary_health_check_exports.R
