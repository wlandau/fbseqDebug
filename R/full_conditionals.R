#' @include full_conditionals_check.R full_conditionals_utils.R
NULL

#' @title Function \code{sample_full_conditionals}
#' @description Get MCMC libraries using the \code{fbseq} package.
#' @export
#' @param dir directory for output files
#' @param counts RNA-seq count data matrix
#' @param design Design matrix. See the \code{fbseq} package vignette for details.
#' @param starts A \code{fbseq::Starts} object of MCMC starting values.
#' @param priors name of prior distributions on the phi's, alpha's, and delta's.
sample_full_conditionals = function(dir, counts, design, starts = Starts(), priors = "Laplace"){
  stopifnot(priors %in% alternate_priors())
  runtimes = NULL

    L = dim(design)[2]
    N = dim(counts)[2]
    G = dim(counts)[1]
    ns = sample.int(N, 12)
    gs = sample.int(G, 12)
    nse = sample.int(N, 3)
    gse = sample.int(G, 4)

   configs = Configs(diag = "none", ess = 0, burnin = 1e3, iterations = 1e4, thin = 0, priors = priors,
                                libraries_return = ns, genes_return = gs, libraries_return_epsilon = nse, genes_return_epsilon  = gse,
                                parameter_sets_return = "beta", parameter_sets_update = "beta")


  for(l in 1:L){
    v = paste0("beta_", l)
    print(v)
    configs@effects_update = l
    chain = Chain(counts, design, configs, starts)
    chain@xiStart = runif(length(chain@xiStart), 0.5, 1.5)
    file = paste0(dir, "chains/", v, ".rds")
    t0 = proc.time()
    chain = fbseq(chain)
    t1 = proc.time() - t0
    print(t1)
    runtimes = rbind(runtimes, t1)
    saveRDS(chain, file)
  }

  vars = setdiff(parameters(), "beta")
  for(v in vars){
    print(v)
    configs@parameter_sets_return = configs@parameter_sets_update = v
    chain = Chain(counts, design, configs, starts)
    chain@xiStart = runif(length(chain@xiStart), 0.5, 1.5)
    file = paste0(dir, "chains/", v, ".rds")
    t0 = proc.time()
    chain = fbseq(chain)
    t1 = proc.time() - t0
    print(t1)
    runtimes = rbind(runtimes, t1)
    saveRDS(chain, file)
  }
  rownames(runtimes) = c(paste0("beta_", 1:L), vars)
  print(runtimes)
}

#' @title Function \code{plot_full_conditionals}
#' @description Plot MCMC libraries against their full conditional densities and make traceplots.
#' @export
#' @param dir directory for output files
plot_full_conditionals = function(dir){
  setwd(paste0(dir, "plots"))
  for(file in list.files("../chains/")){
    path = paste0("../chains/", file)
    name = gsub(".rds", "", file)
    if(grepl("beta", name)){
      if(file.exists(path)){
        chain = readRDS(path)
        beta_check(chain)
      }
      next
    }

    if(file.exists(path)){
      chain = readRDS(path)
      get(paste0(name, "_check"))(chain)
    }
  }
  setwd("../..")
}

#' @title Function \code{full_conditionals_paschold}
#' @description Get MCMC libraries from Paschold data.
#' @export
#' @param priors name of prior distributions on the phi's, alpha's, and delta's.
full_conditionals_paschold = function(priors = "Laplace"){
  stopifnot(priors %in% alternate_priors())
  dir = paste0(priors, "_full_conditionals_paschold/")
  make_dirs(dir)  
  data(paschold)
  counts = get("paschold_counts")
  design = get("paschold_design")
  chain = Chain(counts, design)
  chain@piStart = c(0.43, 0, 1)
  sample_full_conditionals(dir, counts, design, starts = Starts(chain), priors = priors)
}

#' @title Function \code{full_conditionals_simulated}
#' @description Get MCMC libraries from simulated data.
#' @export
#' @param priors name of prior distributions on the phi's, alpha's, and delta's.
full_conditionals_simulated = function(priors = "Laplace"){
  stopifnot(priors %in% alternate_priors())
  dir = paste0(priors, "_full_conditionals_simulated/")
  make_dirs(dir)
  gen = generate_data()
  sample_full_conditionals(dir, gen$counts, gen$design, starts = gen$truth, priors = priors)
}

#' @title Function \code{full_conditionals}
#' @description Check MCMC libraries against their full conditionals.
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
