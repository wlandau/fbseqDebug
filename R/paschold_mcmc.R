#' @title Function \code{paschold_mcmc}
#' @description Run mcmc on paschold data and save output
#' @export
#' @param prior alternate prior to try for the betas
paschold_mcmc = function(prior = "normal"){
  data(paschold)

  dir = paste0("paschold_", prior)
  if(!file.exists(dir)) dir.create(dir)
  setwd(dir)

  configs = Configs(priors = prior)
  chain = Chain(paschold, configs)
  saveRDS(chain, paste0("chain_begin_", prior, ".rds"))

  chain = fbseq(chain, backend = "OpenMP", threads = 4)
  saveRDS(chain, paste0("chain_", prior, ".rds"))
  flat = mcmc(mcmc_samples(chain))
  pdf(paste0("trace_", prior, ".pdf"))
  plot(flat, density = F)
  dev.off()
  pdf(paste0("density_", prior, ".pdf"))
  plot(flat, trace = F)
  dev.off()
  pl = volcano(chain)
  ggsave(filename = paste0("volcano_", prior, ".pdf"), plot = pl)

  setwd("..")
}
