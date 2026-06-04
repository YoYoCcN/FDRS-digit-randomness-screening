# Install required R packages for FDRS

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 1000)

cran_pkgs <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "scales",
  "ranger",
  "glmnet",
  "e1071",
  "isotree",
  "pROC",
  "ggrepel"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  } else {
    message("Package already installed: ", pkg)
  }
}

message("All required packages have been checked.")
