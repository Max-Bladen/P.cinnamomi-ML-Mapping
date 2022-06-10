# Structure of repository

-   Report
    -   Figures
        -   contains all figures found in `Main Report.Rmd` and
            `Supplementary Figures.Rmd`. PNG format.
    -   Tables
        -   contains all information shown in the tables of
            `Main Report.Rmd`. `.rda` format.
    -   `*.Rmd`s and associated `*pdf`s for the main report and
        supplementary figures
-   Analysis
    -   Data
        -   Labelled.Data
            -   data in original format (`.xlsx`) and intermediary
                format (`Long.Lat.Class.csv`). Refer to
                `Analysis/Processing.R` for the derivation of the
                latter.
        -   Spatial.Data
            -   clipped `.TIF` files for the remote sensing derived
                features. Additionally, *TPI*, *Aspect* and *Slope* were
                derived through QGIS and are found here.
        -   Sun.Data
            -   measurements of solar elevation and azimuth. Refer to
                `Analysis/Feature.Engineering.R` for *PRR* calculation
                which uses this data.
        -   three `.rda` files used to transfer data between the `.R`
            files in `Analysis`.
    -   Internals
        -   two `.R` files contain internally built functions and
            package dependencies
    -   Metrics
        -   metrics objects of all eight models (best seven models by
            type plus the final model). `.rda` files.
    -   Four `.R` scripts used for all analysis, model building and
        visualisations.

# Reproducibility

To reproduce, work through `.R` files in `Analysis` in the order:

> `Feature.Engineering.R` -> `Processing.R` -> `Unsupervised.Analysis.R`
> -> `Supervised.Analysis.R`

Ensure you set your working directory to where-ever the
`C:/...P. cinnamomi ML Mapping/Analysis` folder is located at the start
of each script. Large parts of `Supervised.Analysis.R` are commented
out. This script loads previously generated output (5-fold, 100-repeat)
rather than running RCV method as it can take a lot of time. Methods
used to produce this output is commented out but can easily be rerun by
uncommenting these sections.
