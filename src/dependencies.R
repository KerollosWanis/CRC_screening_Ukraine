requiredPackages = c(
  'tidyverse',
  'data.table',
  'pbapply'
)
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}