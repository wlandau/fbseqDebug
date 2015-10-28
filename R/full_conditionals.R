#' @include full_conditionals_check.R full_conditionals_utils.R
NULL

#' @title Function \code{sample_full_conditionals}
#' @description Get MCMC samples using the \code{fbseq} package.
#' @export
#' @param dir directory for output files
#' @param counts RNA-seq count data matrix
#' @param group Experimental design. A vector of integers,
#' one for each RNA-seq sample/library, denoting the genetic
#' variety of that sample. You must use 1 for parent 1, 2 for the hybrid,
#' and 3 for parent 2.
#' @param starts A \code{fbseq::Starts} object of MCMC starting values.
#' @param priors name of prior distributions on the phi's, alpha's, and delta's.
sample_full_conditionals = function(dir, counts, group, starts = Starts(), priors = "Laplace"){
  stopifnot(priors %in% alternate_priors())
  runtimes = NULL
  vars = parameters()
  for(v in vars){
    print(v)
    N = dim(counts)[2]
    G = dim(counts)[1]
    ns = sample.int(N, 12)
    gs = sample.int(G, 12)
    nse = sample.int(N, 3)
    gse = sample.int(G, 4)

    configs = Configs(diag = "none", ess = 0, burnin = 1e3, iterations = 1e4, thin = 0, returns = v, updates = v, 
                                 samples_return = ns, features_return = gs, samples_return_eps = nse, features_return_eps = gse,
                                 phiPrior = priors, alpPrior = priors, delPrior = priors)
    chain = Chain(counts, group, configs, starts)
    file = paste0(dir, "chains/", v, ".rds")
    t0 = proc.time()
    chain = fbseq(chain)
    t1 = proc.time() - t0
    print(t1)
    runtimes = rbind(runtimes, t1)
    saveRDS(chain, file)
  }
  rownames(runtimes) = vars
  print(runtimes)
}

#' @title Function \code{plot_full_conditionals}
#' @description Plot MCMC samples against their full conditional densities and make traceplots.
#' @export
#' @param dir directory for output files
plot_full_conditionals = function(dir){
  setwd(paste0(dir, "plots"))
  for(name in parameters()){
    file = paste0("../chains/", name, ".rds")
    if(file.exists(file)){
      chain = readRDS(file)
      get(paste0(name, "_check"))(chain)
    }
  }
  setwd("../..")
}

#' @title Function \code{full_conditionals_paschold}
#' @description Get MCMC samples from Paschold data.
#' @export
#' @param priors name of prior distributions on the phi's, alpha's, and delta's.
full_conditionals_paschold = function(priors = "Laplace"){
  stopifnot(priors %in% alternate_priors())
  dir = paste0(priors, "_full_conditionals_paschold/")
  make_dirs(dir)  
  data(paschold_counts)
  data(paschold_group)
  counts = get("paschold_counts")
  group = get("paschold_group")
  sample_full_conditionals(dir, counts, group, priors = priors)
}

#' @title Function \code{full_conditionals_simulated}
#' @description Get MCMC samples from simulated data.
#' @export
#' @param priors name of prior distributions on the phi's, alpha's, and delta's.
full_conditionals_simulated = function(priors = "Laplace"){
  stopifnot(priors %in% alternate_priors())
  dir = paste0(priors, "_full_conditionals_simulated/")
  make_dirs(dir)
  gen = generate_data(phiPrior = priors, alpPrior = priors, delPrior = priors)
  sample_full_conditionals(dir, gen$counts, gen$group, starts = gen$truth, priors = priors)
}

#' @title Function \code{full_conditionals}
#' @description Check MCMC samples against their full conditionals.
#' @export
#' @param priors names of prior distributions to try for phi, alpha, and delta.
full_conditionals = function(priors = alternate_priors()){
  for(prior in priors){
    stopifnot(prior %in% alternate_priors())
    full_conditionals_paschold(priors = prior)
    plot_full_conditionals(paste0(prior, "_full_conditionals_paschold/"))
    full_conditionals_simulated(priors = prior)
    plot_full_conditionals(paste0(prior, "_full_conditionals_simulated/"))
  }
}
