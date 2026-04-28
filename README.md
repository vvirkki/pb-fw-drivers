## pb-fw-drivers

### Code repository for

### "Regionally divergent drivers behind transgressions of the freshwater change planetary boundary"

Vili Virkki, Lauren Seaby Andersen, Sofie te Wierik, Dieter Gerten, Miina Porkka

**Accepted for publication in *Nature Communications***\
ADD LINK

**Output data available in a Zenodo repository**\
<https://doi.org/10.5281/zenodo.19663530>

**Corresponding author & repository author**\
Vili Virkki (vili.virkki@uef.fi)

**Please cite the *Nature Communications* publication if using code from this repository or data from the Zenodo repository in another publication.**

### Repository structure and workflow

The complete workflow, including data preparation, analysis, and production of final outputs, is described in the bash script `run.sh`, calling various R scripts (summarised below) in the intended, reproducible order. However, acquiring hydrological simulation data and other input data are not described in `run.sh`. The data preparation and output scripts additionally perform process logging to plain text files that are created automatically.

#### Data

Upon execution, the R scripts read from `Data`, where hydrological simulation data from ISIMIP and other required data should be stored. Each placeholder folder in `Data` has an associated readme explaining the contents required for the analysis.

`Data/ardiv` contains polygons for two levels (level 3 and level 4) of [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins), which are used for the analyses. Additionally, `Data/ensemble_selection` contains csv files that are used to select ensemble members to be used for the analysis. These two folders and their contents should not be altered.

The R scripts will otherwise create their required folder structures and add intermediate outputs within `Data` and into a new folder `output` that will be created in the repository root. Final outputs and graphics are written to additional new folders.

Hydrological simulation data from ISIMIP are not given with this repository, but information on acquiring them is given below.

#### R

`00_` to `06_`: prepare data based on ISIMIP NetCDFs and perform analysis as described in the manuscript. Omitting `02_` from the file naming sequence is intentional.

`output01_`: write intermediate output csv files depicting global and regional land area with local deviations.

`plot00_`: helper functions for plotting and final outputs.

`plot01_`: plot global time series of land area with local deviations.

`plot02_`: plot maps of regional land area with local deviations and scenario contributions.

`plot03_`: plot synthesising analysis of scenario contributions on streamflow and soil moisture deviations.

`main_`, `future_`, and `supplementary_`: control workflow for preparing data and intermediate outputs. Common properties include:

-   assume working directory in the repository root

-   silently overwrite intermediate outputs in `Data`

-   call required libraries (`libraries.R`) and other scripts within

`temp01_`: export data relevant to [Planetary Health Check](https://www.planetaryhealthcheck.org).

`temp02_`: diagnose the infilling of missing values in data preparation.

`isimip_3a_` and `isimip_3b_`: helper scripts for downloading hydrological simulation data from the ISIMIP repository.

#### configs

This folder contains plaintext files that describe parameters used in the analysis. The format of the files follows a two-piece convention (separated by semicolon) in which the first piece denotes the parameter name and the second piece denotes the parameter value.

The config files describing analysis parameters (tagged with `main`) are read by `R/00_configure_run.R` and a list of parameters is then passed to data preparation (`01_` to `06_`) and output (`output01_`) scripts.

#### Other files

`environment.txt`: Description of the computing environment, including operating system and versions of R and R libraries, with which the analysis has been run.

### Hydrological simulation data acquisition

Hydrological simulation data should be acquired from [ISIMIP](https://data.isimip.org/search/) and placed into a folder in `Data` with the following structure. Scripts `R/isimip_3a_` and `R/isimip_3b_` should help in this task.

```         
CONFIG_RAW_DATA_ROOT
└── CONFIG_VARIABLE
    ├── CONFIG_SCENARIO_LABEL_1
    │   └── raw
    │       └── add files here
    └── CONFIG_SCENARIO_LABEL_2
    │   └── raw
    │       └── add files here
    │
    │   ...
    │
    └── CONFIG_SCENARIO_LABEL_N
        └── raw
            └── add files here
```

### Changelog and versioning

Should the repository be updated with new versions, main changes will be briefly summarised here. Commits describing new versions are marked respectively with Git tags.

#### v1.0.0.

Version used in producing the results shown in the *Nature Communications* manuscript (ADD LINK).

### License

Attribution 4.0 International (CC BY 4.0)\
<https://creativecommons.org/licenses/by/4.0/>
