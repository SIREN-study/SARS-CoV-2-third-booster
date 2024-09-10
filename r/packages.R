# Running this R file will ensure all of the required packages and dependencies are installed.

install.packages("pak")

pak::pkg_install(
  c(
    "broom", "dplyr", "finalfit", "forcats", "future",
    "future.callr", "ggnewscale", "ggplot2", "ggpp",
    "gt", "here", "janitor", "lubridate", "msm",
    "patchwork", "readstata13", "readxl", "scales",
    "survival", "tidyr"
  )
)
