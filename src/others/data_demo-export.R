# Load packages
source(file.path("src", "required.R"))

# Prepare data
subset_qc <- TRUE
source(here("src", "02a_training_data-preparation.R"))

# Check data
spe
vars
spe_full
vars_full

# Export data
write_tsv(spe,  here("..", "embarcadero-demo", "data", "spe.csv"))
write_tsv(vars, here("..", "embarcadero-demo", "data", "env.csv"))
write_tsv(spe_full,  here("..", "embarcadero-demo", "data", "spe_full.csv"))
write_tsv(vars_full, here("..", "embarcadero-demo", "data", "env_full.csv"))
