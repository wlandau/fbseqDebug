#' @title Function \code{paschold_mcmc}
#' @description Run mcmc on paschold data and save output
#' @export
#' @param priors alternate prior to try for phi, alpha, and delta
#' @param diag can be "geweke", "gelman", or "none"
paschold_mcmc = function(priors = c("normal", alternate_priors()), diag = "gelman"){
  data(paschold)

  for(prior in priors){
    dir = paste0("paschold_", prior, "_", diag)
    if(!file.exists(dir)) dir.create(dir)
    setwd(dir)

    gs = sample.int(dim(paschold@counts)[1], 2)
    ns = sample.int(dim(paschold@counts)[2], 2)

    configs = Configs(diag = diag, max_attempts = 10, priors = c("normal", rep(prior, 4)), nchains_diag = 4, 
      genes_return = gs, genes_return_epsilon = gs, libraries_return_epsilon = ns)

    if(prior == "horseshoe"){
      configs@thetaStart[2:5] = 0
      configs@effects_update_theta = 1
    }

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
}
