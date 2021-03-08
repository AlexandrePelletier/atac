library(data.table)
options(repos = c(CRAN="https://pbil.univ-lyon1.fr/CRAN/"))

renv::use_python()
install.packages("reticulate")
install.packages("here")
reticulate::py_install("scipy")

system(paste(here::here("renv/python/virtualenvs/renv-python-3.7.3/bin/python"), "-m ensurepip"))


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
