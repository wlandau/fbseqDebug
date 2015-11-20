#' @include compare_points.R
NULL

#' @title Function \code{simulated_mcmc}
#' @description Simulate data and try to recapture true parameter values
#' @export
#' @param prior alternate prior to try for the betas
#' @param diag can be "gelman" or "none"
simulated_mcmc = function(prior = "Laplace", diag = "gelman"){
  libraries = 16
  genes = 3.5e4

  dir = paste0("sim_", prior, "_", diag)
  if(!file.exists(dir)) dir.create(dir)
  setwd(dir)

  s = scenario_heterosis_model(genes = genes)
  saveRDS(s, paste0("scenario_", prior, ".rds"))

  configs = Configs(diag = diag, priors = prior)
  chain = Chain(s, configs)
  chain = fbseq(chain)

  saveRDS(chain, paste0("chain_", prior, ".rds"))
  comp = compare_points(chain, s@supplement$truth)
  write.table(comp, file = paste0("mcmc_vs_truth_", prior, ".txt"))
  flat = mcmc(flatten(chain))
  pdf(paste0("trace_", prior, ".pdf"))
  plot(flat, density = F)
  dev.off()
  pdf(paste0("density_", prior, ".pdf"))
  plot(flat, trace = F)
  dev.off()
  ggsave(filename = paste0("volcano_", prior, ".pdf"), plot = volcano(chain))

  setwd("..")
}
