#' @include compare_points.R
NULL

#' @title Function \code{simulated_mcmc}
#' @description Simulate data and try to recapture true parameter values
#' @export
#' @param prior alternate prior to try for the betas
simulated_mcmc = function(prior = "normal"){
  libraries = 16
  genes = 3e4

  dir = paste0("sim_", prior)
  if(!file.exists(dir)) dir.create(dir)
  setwd(dir)

  s = scenario_heterosis_model(genes = genes)
  saveRDS(s, paste0("scenario_", prior, ".rds"))

  configs = Configs(priors = prior, burnin = 1e4, thin = 10, iterations = 1e3)
  chain = Chain(s, configs, starts = Starts(h = 0))
  saveRDS(chain, paste0("chain_begin_", prior, ".rds"))

  chain = fbseq(chain, backend = "OpenMP", threads = 4)
  saveRDS(chain, paste0("chain_", prior, ".rds"))
  comp = compare_points(chain, s@supplement$truth)
  write.table(comp, file = paste0("mcmc_vs_truth_", prior, ".txt"))
  flat = mcmc(mcmc_samples(chain))
  pdf(paste0("trace_", prior, ".pdf"))
  plot(flat, density = F)
  dev.off()
  pdf(paste0("density_", prior, ".pdf"))
  plot(flat, trace = F)
  dev.off()
  ggsave(filename = paste0("volcano_", prior, ".pdf"), plot = volcano(chain))

  setwd("..")
}
