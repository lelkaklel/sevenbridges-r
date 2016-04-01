## ----include=FALSE-------------------------------------------------------
knitr::opts_chunk$set(eval = FALSE)
library(sevenbridges)

## ------------------------------------------------------------------------
if(!require("devtools", quietly = TRUE)){
    install.packages("devtools") 
}
## Install from github for development version  
source("http://bioconductor.org/biocLite.R")
devtools::install_github("tengfei/sevenbridges", build_vignettes=TRUE, 
  repos=BiocInstaller::biocinstallRepos(),
  dependencies=TRUE)

## ------------------------------------------------------------------------
browseVignettes(package = 'sevenbridges')

## ------------------------------------------------------------------------
a <- Auth(token = "fake_token", url = "https://cgc-api.sbgenomics.com/v2/")

## ------------------------------------------------------------------------
a <- Auth(user = "tengfei", platform = "cgc")

## ------------------------------------------------------------------------
(b <- a$billing())
## a single billing group is showing

## ------------------------------------------------------------------------
(p <- a$project_new("hackathon", 
                    billing_group_id = b$id, 
                    description = "This project is for CGC hackathon"))

## ------------------------------------------------------------------------
## p$delete()

## ---- comment = "", eval = TRUE------------------------------------------
rbx <- Tool(id = "runif", 
            label = "runif",
            hints = requirements(docker(pull = "rocker/r-base"), 
                                 cpu(1), mem(2000)), 
            baseCommand = "Rscript -e 'runif(100)'", 
            stdout = "output.txt",
            outputs = output(id = "random", glob = "*.txt"))
rbx$toJSON(pretty = TRUE)

## ---- eval = TRUE--------------------------------------------------------
## provide scripts
## Make a new script file
fd <- fileDef(name = "runif.R",
              content = "sed.seed(1)
                   runif(100)")
rbx <- Tool(id = "runif", 
            label = "runif",
            hints = requirements(docker(pull = "rocker/r-base"), 
                                 cpu(1), mem(2000)),
            requirements = requirements(fd),
            baseCommand = "Rscript runif.R", ## run script you created.
            stdout = "output.txt",
            outputs = output(id = "random", glob = "*.txt"))   

## ---- eval = TRUE--------------------------------------------------------
fl <- system.file("docker/rnaseqGene", "Dockerfile", 
                  package = "sevenbridges")

## ----comment='', eval = TRUE---------------------------------------------
cat(readLines(fl), sep = '\n')

## ---- eval = TRUE--------------------------------------------------------
fl <- system.file("docker/rnaseqGene/src", "performDE.R", 
                  package = "sevenbridges")

## ----comment='', eval = TRUE---------------------------------------------
cat(readLines(fl), sep = '\n')

## ---- eval = TRUE--------------------------------------------------------
fl <- system.file("docker/rnaseqGene/report", "rnaseqGene.Rmd", 
                  package = "sevenbridges")

## ----comment='', eval = TRUE---------------------------------------------
cat(readLines(fl, n = 50), sep = '\n')

## ---- eval = TRUE--------------------------------------------------------

rbx <- Tool(id = "rnaseqGene", 
            label = "rnaseqgene",
            description = "A RNA-seq Differiencial Expression Flow and Report",
            hints = requirements(docker(pull = "tengfei/rnaseqgene"), cpu(1), mem(2000)), 
            baseCommand = "performDE.R", 
            inputs = list(
                input(
                    id = "bamfiles", label = "bam files",
                    description = "a list of bam files",
                    type = "File...",  ## or type = ItemArray("File")
                    prefix = "--bamfiles",
                    itemSeparator = ","
                ), 
                input(
                    id = "design", label = "design matrix",
                    type = "File",
                    prefix = "--design"
                ),
                input(
                    id = "gtffile", label =  "gene feature files",
                    type = "File",
                    prefix = "--gtffile"
                ),
                input(
                    id = "format", label =  "report foramt html or pdf",
                    type = enum("format", c("pdf", "html")),
                    prefix = "--format"
                )
            ),
            outputs = list(
                output(id = "report", label = "report", 
                       description = "A reproducible report created by Rmarkdown",
                       glob = Expression(engine = "#cwl-js-engine",
                                         script = "x = $job[['inputs']][['format']];
                                                  if(x == 'undefined' || x == null){
                                                   x = 'html';
                                                    };
                                                  'rnaseqGene.' +  x")),
                output(id = "heatmap", label = "heatmap", 
                       description = "A heatmap plot to show the Euclidean distance between samples",
                       glob = "heatmap.pdf"),
                output(id = "count", label = "count", 
                       description = "Reads counts matrix",
                       glob = "count.csv"),
                output(id = "de", label = "Differential expression table", 
                       description = "Differential expression table",
                       glob = "de.csv")
                ))

## ----comment = "", eval = TRUE-------------------------------------------
rbx
rbx$toJSON(pretty = TRUE)
rbx$toJSON()
## or write to external file
## fl <- "~/Downloads/rnaseqGene.json"
## write(rbx$toJSON(pretty = TRUE), fl)

## ------------------------------------------------------------------------
## add App you just created 
(rna.app <- p$app_add("rnaseqgene", rbx))

## ------------------------------------------------------------------------
fl <- system.file("extdata", "sample1.fastq", package = "sevenbridges")
(p <- a$project(id = "tengfei/quickstart"))
## by default load .meta for the file
p$upload(fl, overwrite = TRUE)
## pass metadata
p$upload(fl, overwrite = TRUE, metadata = list(library_id = "testid2", platform = "Illumina x11"))
## rename
p$upload(fl, overwrite = TRUE, name = "sample_new_name.fastq", 
         metadata = list(library_id = "new_id"))

## upload folder
dir.ext <- system.file("extdata", package = "sevenbridges")
p$upload(dir.ext, overwrite = TRUE)

## upload file list
fls <- list.files(dir.ext, recursive = TRUE, full.names = TRUE)
p$upload(fls, overwrite = TRUE)

## ---- comment = "", eval = TRUE------------------------------------------
download.fl <- system.file("extdata/download.txt", package = "sevenbridges")
cat(readLines(download.fl), sep = '\n')

## ------------------------------------------------------------------------
## get file id you need as inout
(bamfiles.in <- p$file(".bam"))
(design.in <- p$file("sample_table.csv"))
(gtf.in <- p$file("Homo_sapiens.GRCh37.75_subset.gtf"))

## ------------------------------------------------------------------------
## add a new Task
(tsk <- p$task_add(name = "RNA DE report new", 
                   description = "RNA DE analysis report", 
                   app = rna.app$id,
                   inputs = list(bamfiles = bamfiles.in, 
                                 design = design.in,
                                 gtffile = gtf.in)))

## don't forget to run a draft task
tsk$run()

## ------------------------------------------------------------------------
## monitor the task
## tsk$monitor()

## ------------------------------------------------------------------------
setTaskHook("completed", function(){
    tsk$download("~/Downloads")
})
tsk$monitor()

## ------------------------------------------------------------------------
tsk$download("~/Downloads")

## ----comment = "", eval = TRUE-------------------------------------------
fl <- system.file("docker/reporttool/rabix/reporttool.json", 
                  package = "sevenbridges")
cat(readLines(fl), sep= "\n")

## ------------------------------------------------------------------------
## directly add json file
p <- a$project(id = "tengfei/hackathon")
(report.app <- p$app_add("report-tool", fl))

## ----comment = "", eval = TRUE-------------------------------------------
fl <- system.file("docker/reporttool/Dockerfile", 
                  package = "sevenbridges")
cat(readLines(fl), sep= "\n")

## ----comment = "", eval = FALSE------------------------------------------
fl <- system.file("docker/reporttool/src/report.R", 
                  package = "sevenbriges")
cat(readLines(fl), sep= "\n")

## ----comment = "", eval = FALSE------------------------------------------
fl <- system.file("docker/reporttool/rabix/generator.R", 
                  package = "sevenbriges")
cat(readLines(fl), sep= "\n")

## ------------------------------------------------------------------------
browseVignettes("sevenbridges")

