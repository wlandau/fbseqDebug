#' @title Function \code{paschold_mcmc}
#' @description Run mcmc on paschold data and save output
#' @export
#' @param priors alternate prior to try for phi, alpha, and delta
#' @param diag can be "geweke", "gelman", or "none"
paschold_mcmc = function(priors = c("normal", alternate_priors()), diag = "none"){
  data(paschold_counts)
  data(paschold_group)
  counts = get("paschold_counts")
  group = get("paschold_group")  
  features = dim(counts)[1]

  for(prior in priors){
    dir = paste0("paschold_", prior, "_", diag)
    if(!file.exists(dir)) dir.create(dir)
    setwd(dir)

    configs = Configs(diag = diag, max_attempts = 10, phiPrior = prior, alpPrior = prior, delPrior = prior)
    chain = Chain(counts, group, configs)
    chain = fbseq(chain)

    saveRDS(chain, paste0("chain_", prior, ".rds"))
    flat = mcmc(flatten(chain))
    pdf(paste0("trace_", prior, ".pdf"))
    plot(flat, density = F)
    dev.off()
    pdf(paste0("density_", prior, ".pdf"))
    plot(flat, trace = F)
    dev.off()

    setwd("..")
  }
}