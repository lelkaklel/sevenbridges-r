#!/usr/local/bin/Rscript
library(rmarkdown)
render("~/Code/svnrepos/bioc-devel/sevenbridges/vignettes/bioc-workflow.Rmd",
       output_dir = ".")
