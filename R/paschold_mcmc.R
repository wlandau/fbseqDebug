#' @title Function \code{paschold_mcmc}
#' @description Run mcmc on paschold data and save output
#' @export
#' @param prior alternate prior to try for the betas
#' @param diag can be "gelman" or "none"
paschold_mcmc = function(prior = "normal", diag = "gelman"){
  data(paschold)

  dir = paste0("paschold_", prior, "_", diag)
  if(!file.exists(dir)) dir.create(dir)
  setwd(dir)

  configs = Configs(diag = diag, priors = prior)
  chain = Chain(paschold, configs)
  chain = fbseq(chain)

  saveRDS(chain, paste0("chain_", prior, ".rds"))
  flat = mcmc(flatten(chain))
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
