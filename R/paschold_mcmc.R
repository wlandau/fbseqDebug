#' @title Function \code{paschold_mcmc}
#' @description Run mcmc on paschold data and save output
#' @export
#' @param priors alternate prior to try for phi, alpha, and delta
#' @param diag can be "geweke", "gelman", or "none"
#' @param constrain_theta set the constrain_theta slot int the \code{fbseq::Chain} object
paschold_mcmc = function(priors = c("normal", alternate_priors()), diag = "gelman", constrain_theta = c(T, F)){
  data(paschold)

  for(prior in priors) for(con in constrain_theta){
    dir = paste0("paschold_", prior, "_", diag, "_", con)
    if(!file.exists(dir)) dir.create(dir)
    setwd(dir)

    configs = Configs(diag = diag, priors = prior, constrain_theta = con)
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
