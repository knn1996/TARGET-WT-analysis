library(MOFA2)
library(ggplot2)

mofa_input <- readRDS("data/mofa_input.rds")

N_FACTORS <- 12
SEEDS <- c(42, 101, 2024, 7, 555)
VAR_THRESHOLD <- 0.01
dir.create("models", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

build_and_run <- function(seed, outfile) {
  obj <- create_mofa(mofa_input)

  data_opts <- get_default_data_options(obj)
  data_opts$scale_views <- TRUE

  model_opts <- get_default_model_options(obj)
  model_opts$num_factors <- N_FACTORS
  model_opts$likelihoods <- c(expression = "gaussian", methylation = "gaussian")

  train_opts <- get_default_training_options(obj)
  train_opts$convergence_mode <- "slow"
  train_opts$drop_factor_threshold <- VAR_THRESHOLD
  train_opts$seed <- seed

  obj <- prepare_mofa(obj,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )

  run_mofa(obj, outfile = outfile, use_basilisk = TRUE)
}

model_main <- build_and_run(SEEDS[1], "models/mofa_wt_seed42.hdf5")
saveRDS(model_main, "models/mofa_wt_main.rds")

pdf("results/mofa_data_overview.pdf", width = 7, height = 4)
print(plot_data_overview(model_main))
dev.off()

var_exp <- get_variance_explained(model_main)
write.csv(var_exp$r2_per_factor$group1, "results/mofa_variance_per_factor.csv")
write.csv(var_exp$r2_total$group1, "results/mofa_variance_total.csv")
print(var_exp$r2_per_factor$group1)
print(var_exp$r2_total$group1)

pdf("results/mofa_variance_explained.pdf", width = 6, height = 5)
print(plot_variance_explained(model_main, max_r2 = 15))
print(plot_variance_explained(model_main, plot_total = TRUE)[[2]])
dev.off()

pdf("results/mofa_factor_correlation.pdf", width = 5, height = 5)
plot_factor_cor(model_main)
dev.off()

Z_main <- get_factors(model_main)$group1

stability <- list()
for (s in SEEDS[-1]) {
  m <- build_and_run(s, file.path("models", paste0("mofa_wt_seed", s, ".hdf5")))
  Z <- get_factors(m)$group1
  cm <- abs(cor(Z_main, Z, use = "pairwise.complete.obs"))
  stability[[as.character(s)]] <- apply(cm, 1, max)
}

stability_df <- as.data.frame(do.call(cbind, stability))
stability_df$min_abs_cor <- apply(stability_df, 1, min)
stability_df$stable <- stability_df$min_abs_cor >= 0.8
write.csv(stability_df, "results/mofa_factor_stability.csv")
print(stability_df)
