#' @title Function \code{credible_debug}
#' @description Run mcmc on paschold data and check crude credible intervals against actual quantile
#' credible intervals
#' @export
#' @param priors alternate prior to try for phi, alpha, and delta
credible_debug = function(priors = c("normal", special_beta_priors())){
  data(paschold)

  for(prior in priors){
    dir = paste0("credible_", prior)
    if(!file.exists(dir)) dir.create(dir)
    setwd(dir)

    N = dim(paschold@counts)[2]
    G = dim(paschold@counts)[1]

    configs = Configs(priors = prior, iterations = 1e4, thin = 1, burnin = 1e3,
      genes_return = sample.int(G, 12), libraries_return = sample.int(N, 12),
      genes_return_epsilon = sample.int(G, 4), libraries_return_epsilon = sample.int(N, 3))
    chain = Chain(paschold, configs)
    out = fbseq(chain)
    saveRDS(out, "out.rds")

    samples = mcmc_samples(out)
    estimates = estimates(out)

    for(name in colnames(samples)){
      x = samples[,name]
      realmean = mean(x)
      qlow = quantile(x, 0.025)
      qhigh = quantile(x, 0.975)
      if(length(unique(x)) > 3)
        x = x[x > quantile(x, 0.0025) & x < quantile(x, 0.9975)]
      df = data.frame(x = x)

      pl = ggplot(data = df) + 
        geom_histogram(mapping = aes(x = x, y = ..density..), binwidth = diff(range(x))/30, alpha = 0.6) +
        labs(title = paste0(name, "\n")) + 
        xlab("\nMCMC libraries") +
        ylab("Density\n") +
        geom_vline(xintercept = estimates[name, "mean"], color = "green", linetype = 6) +
        geom_vline(xintercept = realmean, color = "green") +
        geom_vline(xintercept = estimates[name, "lower_ci_0.95"], color = "blue", linetype = 6) +
        geom_vline(xintercept = qlow, color = "blue") +
        geom_vline(xintercept = estimates[name, "upper_ci_0.95"], color = "red", linetype = 6) +
        geom_vline(xintercept = qhigh, color = "red") +
        theme(panel.background = element_rect(fill='white'),
                   panel.border = element_rect(color="black", fill = NA),
                   panel.grid.major = element_line(color="lightgray"),
                   panel.grid.minor = element_blank())

      ggsave(filename = paste0(name, ".pdf"), plot = pl, width = 8, height = 8, dpi = 1600)
    }

    setwd("..")
  }
}