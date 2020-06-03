# Load packages
source(file.path("src", "required.R"))

# Prepare data
subset_qc <- TRUE
source(here("src", "02a_training_data-preparation.R"))

# Check data
spa
env
spe
spa_full
env_full
spe_full

# Remove site from full range sets
env_full <- select(env_full, -site)
spe_full <- select(spe_full, -site)

# Export data
write_tsv(spe, here("..", "embarcadero-demo", "data", "spe.csv"))
write_tsv(spa, here("..", "embarcadero-demo", "data", "spa.csv"))
write_tsv(env, here("..", "embarcadero-demo", "data", "env.csv"))
write_tsv(spe_full, here("..", "embarcadero-demo", "data", "spe_full.csv"))
write_tsv(spa_full, here("..", "embarcadero-demo", "data", "spa_full.csv"))
write_tsv(env_full, here("..", "embarcadero-demo", "data", "env_full.csv"))

