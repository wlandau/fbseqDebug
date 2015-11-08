#' @include compare_points.R
NULL

#' @title Function \code{simulated_mcmc}
#' @description Simulate data and try to recapture true parameter values
#' @export
#' @param priors alternate prior to try for phi, alpha, and delta
#' @param diag can be "geweke" or "gelman"
simulated_mcmc = function(priors = c("normal", alternate_priors()), diag = "gelman"){
  libraries = 16
  genes = 3.5e4
  for(prior in priors){
    dir = paste0("sim_", prior, "_", diag)
    if(!file.exists(dir)) dir.create(dir)
    setwd(dir)

    s = scenario_heterosis_model(genes = genes, libraries = libraries)
    saveRDS(s, paste0("scenario_", prior, ".rds"))

    gs = sample.int(genes, 2)
    ns = sample.int(libraries, 2)

    configs = Configs(diag = diag, max_attempts = 10, priors = prior, nchains_diag = 4,
      genes_return = gs, genes_return_epsilon = gs, libraries_return_epsilon = ns,

      burnin = 1e7
#      effects_update_theta = 1, 
    )
    chain = Chain(s, configs)
#    chain@thetaStart[2:5] = 0
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
