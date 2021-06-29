## Variable selection ####

# To run after main bart script (where data is prepared)

warning("Testing variable selection for 3 species only. This can take a while to run at full scale.")
pred_plot

# Select fewer species
sort(colSums(spe)/nrow(spe), decreasing = TRUE)
spe_sel <- c("sp1", "sp6", "sp17")
vars_sel <- map(spe_sel, ~ NULL)
names(vars_sel) <- spe_sel

# Run variable selection
# variable_selection <- TRUE
if (exists("variable_selection") && isTRUE(variable_selection)) {
  tictoc::tic("total")
  for(sp in spe_sel){
    tictoc::tic(sp)
    set.seed(42)
    message(paste0("Variable selection for ", sp, " (", which(sp == spe_sel), "/", length(spe_sel)), ")")
    # Save plot to png
    png(here("fig", "bart", paste0("x_bart_vars-select_", sp, ".png")))
    step_vars <- variable.step(
      y.data = spe[[sp]],
      x.data = env[xnames],
      iter = 50
    )
    dev.off()
    # Save variables to list
    vars_sel[[sp]] <- step_vars
    tictoc::toc()
  }
  tictoc::toc()
  vars_sel
}
