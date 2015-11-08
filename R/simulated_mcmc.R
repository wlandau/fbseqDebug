#' @include compare_points.R
NULL

#' @title Function \code{simulated_mcmc}
#' @description Simulate data and try to recapture true parameter values
#' @export
#' @param priors alternate prior to try for phi, alpha, and delta
#' @param diag can be "geweke" or "gelman"
#' @param constrain_theta set the constrain_theta slot int the \code{fbseq::Chain} object
simulated_mcmc = function(priors = c("normal", alternate_priors()), diag = "gelman", constrain_theta = c(T, F)){
  libraries = 16
  genes = 3.5e4
  for(prior in priors) for(con in constrain_theta){
    dir = paste0("sim_", prior, "_", diag, "_", con)
    if(!file.exists(dir)) dir.create(dir)
    setwd(dir)

    s = scenario_heterosis_model(genes = genes, libraries = libraries)
    saveRDS(s, paste0("scenario_", prior, ".rds"))

    configs = Configs(diag = diag, priors = prior, constrain_theta = con)
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
}
