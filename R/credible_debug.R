#' @title Function \code{credible_debug}
#' @description Run mcmc on paschold data and check crude credible intervals against actual quantile
#' credible intervals
#' @export
#' @param priors alternate prior to try for phi, alpha, and delta
#' @param diag can be "geweke", "gelman", or "none"
credible_debug = function(priors = c("normal", alternate_priors()), diag = "none"){
  data(paschold_counts)
  data(paschold_group)
  counts = get("paschold_counts")
  group = get("paschold_group")  

  for(prior in priors){
    dir = paste0("credible_", prior)
    if(!file.exists(dir)) dir.create(dir)
    setwd(dir)

    N = dim(counts)[2]
    G = dim(counts)[1]

    configs = Configs(diag = "none", ess = 0, phiPrior = prior, alpPrior = prior, delPrior = prior, 
      features_return = sample.int(G, 12), samples_return = sample.int(N, 12),
      features_return_eps = sample.int(G, 4), samples_return_eps = sample.int(N, 3))
    chain = Chain(counts, group, configs)
    chain = fbseq(chain)
    saveRDS(chain, "chain.rds")

    flat = flatten(chain)
    est = estimates(chain)$params

    for(name in colnames(flat)){
      x = flat[,name]
      qlow = quantile(x, 0.025)
      qhigh = quantile(x, 0.975)
      x = x[x > quantile(x, 0.0025) & x < quantile(x, 0.9975)]
      df = data.frame(x = x)

      pl = ggplot(data = df) + 
        geom_histogram(mapping = aes(x = x, y = ..density..), binwidth = diff(range(x))/30, alpha = 0.6) +
        labs(title = paste0(name, "\n")) + 
        xlab("\nMCMC samples") +
        ylab("Density\n") +
        geom_vline(xintercept = est[name, "lower_ci_0.95"], color = "blue", linetype = 6) +
        geom_vline(xintercept = qlow, color = "green") +
        geom_vline(xintercept = est[name, "upper_ci_0.95"], color = "red", linetype = 6) +
        geom_vline(xintercept = qhigh, color = "black") +
        theme(panel.background = element_rect(fill='white'),
                   panel.border = element_rect(color="black", fill = NA),
                   panel.grid.major = element_line(color="lightgray"),
                   panel.grid.minor = element_blank())

      ggsave(filename = paste0(name, ".pdf"), plot = pl, width = 8, height = 8, dpi = 1600)
    }

    setwd("..")
  }
}