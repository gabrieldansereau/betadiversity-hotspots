install.packages("renv")

library(renv)

renv::init()

renv::status()
renv::diagnostics()
renv::history()
renv::activate()
renv::dependencies()
renv::project()
renv::activate()
renv::load()
renv::settings()

# Enable automatic snapshot (more similar to Julia project manager)
getOption("renv.config.auto.snapshot")
options(renv.config.auto.snapshot = TRUE)
# add to ./.Rprofile to enable by default (in project)

# Test adding a package
# install.packages("randomForest")
# library(randomForest)
# renv::status()
# renv::remove("randomForest")

# Tweak settings
renv::settings
# options(renv.settings.snapshot.type = "explicit") 
options(renv.settings.snapshot.type = "all")
# add to ./.Rprofile