# atac

Analyst: Alexandre Pelletier


## Setup

``` r
renv::use_python()
install.packages("reticulate")
install.packages("here")

system(paste(here::here("renv/python/virtualenvs/renv-python-3.7.3/bin/python"), "-m ensurepip"))
reticulate::py_install("scipy")


use_python <- function(type = "virtualenv", path = here::here()) {
  project_directory <- here::here()
  working_directory <- gsub("PROJECT", "DATATMP", project_directory)
  dir.create(
    path = file.path(working_directory, "python"),
    recursive = TRUE, showWarnings = FALSE, mode = "0775"
  )
  file.symlink(
    from = file.path(working_directory, "python"),
    to = file.path(project_directory, "renv", "python")
  )
  renv::use_python(type = type)
  system(paste(here::here("renv/python/virtualenvs/renv-python-3.7.3/bin/python"), "-m ensurepip"))
}

renv::snapshot()

```

<!--
## Design

``` bash
nohup Rscript scripts/01-design.R > logs/01.log &
```

## Quality Control

``` bash
nohup Rscript -e 'rmarkdown::render(input = here::here("scripts", "02-qc.Rmd"), output_file = "QC.html", output_dir = here::here("reports"), encoding = "UTF-8")' > logs/02.log &
```

## Statistical Analyses

``` bash
nohup Rscript scripts/03-analysis.R > logs/03.log &
```

## Meeting Slides

### 2021-02-22

``` bash
nohup Rscript -e 'rmarkdown::render(input = here::here("scripts", "20210222-meeting.Rmd"), output_dir = here::here("reports"))' > logs/20210222-meeting.log &
```
-->
