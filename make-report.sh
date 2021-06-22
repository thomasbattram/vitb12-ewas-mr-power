#!/bin/bash

Rscript -e "require(knitr); require(markdown); require(bookdown); rmarkdown::render('report/report.Rmd')"