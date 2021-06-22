# VitB12 - DNAm MR power calculations

Data needed to run the analyses:

1. EWAS results (`MA.b12.maternal.main.model.GSM.no.MARBLES.csv` and `MA.b12.newborn.main.model.GSM.csv`)
2. `aries-antenatal-variances.RData`
3. `b12-distribution-edit.xlsx`
4. `vitb12-gwas-data.xlsx` 

## Analyses workflow

0. Check packages to install (can use code below or just look in the scripts)
1. Put all the data (listed above) in a folder named "data"
2. Run through [power-calculations.R](scripts/power-calculations.R)
3. Run `bash make-report.sh`

### Install packages needed

`
if (!require(attachment)) install.packages("attachment")
pkgs_from_script <- att_from_rscripts(path = "scripts")
pkgs_from_report <- att_from_rmds(path = "report", inside_rmd = TRUE)
pkgs <- c(pkgs_from_script, pkgs_from_report)
lapply(pkgs, function(pkg) {if (!require(pkg, character.only = T)) install.packages(pkg)})
`
