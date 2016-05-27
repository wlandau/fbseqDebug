#' @title Function \code{paschold_mcmc}
#' @description Run mcmc on paschold data and save output
#' @export
#' @param prior alternate prior to try for the betas
paschold_mcmc = function(prior = "normal"){
  data(paschold)
  paschold@counts = paschold@counts[1:30,]

  dir = paste0("paschold_", prior)
  if(!file.exists(dir)) dir.create(dir)
  setwd(dir)

  configs = Configs(priors = prior, iterations = 1000, burnin = 100, thin = 1)
  chain = Chain(paschold, configs)
  saveRDS(chain, paste0("chain_begin_", prior, ".rds"))

  chain = fbseq(chain, backend = "OpenMP", threads = 4)
  saveRDS(chain, paste0("chain_", prior, ".rds"))

  setwd("..")
}
