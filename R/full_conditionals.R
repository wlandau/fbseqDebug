#' @include full_conditionals_check.R full_conditionals_utils.R
NULL

#' @title Function \code{sample_full_conditionals}
#' @description Get MCMC libraries using the \code{fbseq} package.
#' @export
#' @param dir directory for output files
#' @param scenario a \code{fbseq::Scenario} object.
#' @param starts A \code{fbseq::Starts} object of MCMC starting values.
#' @param priors name of prior distributions on the phi's, alpha's, and delta's.
sample_full_conditionals = function(dir, scenario, starts = Starts(), priors = "Laplace"){
  stopifnot(priors %in% special_beta_priors())
  runtimes = NULL

  L = dim(scenario@design)[2]
  N = dim(scenario@counts)[2]
  G = dim(scenario@counts)[1]
  ns = sample.int(N, 12)
  gs = sample.int(G, 12)
  nse = sample.int(N, 3)
  gse = sample.int(G, 4)

  configs = Configs(burnin = 1e3, iterations = 1e4, thin = 1, priors = "normal",
                               libraries_return = ns, genes_return = gs, libraries_return_epsilon = nse, 
                               genes_return_epsilon  = gse,
                               parameter_sets_return = "beta", parameter_sets_update = "beta")

  starts@xi = runif(G*L, 0.5, 1.5)

  for(l in 1:L){
    v = paste0("beta_", l)
    print(v)
    configs@effects_update_beta = l
    chain = Chain(scenario, configs, starts)
    file = paste0(dir, "chains/", v, ".rds")
    t0 = proc.time()
    chain = fbseq(chain, backend = "OpenMP", threads = 4, additional_chains = 0)
    t1 = proc.time() - t0
    print(t1)
    runtimes = rbind(runtimes, t1)
    saveRDS(chain, file)
  }

  vars = setdiff(parameters(), "beta")
  for(v in vars){
    print(v)
    configs@parameter_sets_return = configs@parameter_sets_update = v
    configs@priors = ifelse(v == "xi", priors, "normal")
    chain = Chain(scenario, configs, starts)
    file = paste0(dir, "chains/", v, ".rds")
    t0 = proc.time()
    chain = fbseq(chain, backend = "OpenMP", threads = 4, additional_chains = 0)
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
  stopifnot(priors %in% special_beta_priors())
  dir = paste0(priors, "_full_conditionals_paschold/")
  make_dirs(dir)  
  data(paschold)
  paschold@counts = paschold@counts[1:50,]
  sample_full_conditionals(dir, paschold, priors = priors)
}

#' @title Function \code{full_conditionals_simulated}
#' @description Get MCMC libraries from simulated data.
#' @export
#' @param priors name of prior distributions on the phi's, alpha's, and delta's.
full_conditionals_simulated = function(priors = "Laplace"){
  stopifnot(priors %in% special_beta_priors())
  dir = paste0(priors, "_full_conditionals_simulated/")
  make_dirs(dir)
  s = scenario_heterosis_model(genes = 50,
    truth = Starts(h = c(-0.5, 0, 0.5), nu = 2.812224, tau = 0.006780517, sigmaSquared = c(1, 
        0.03724313, 0.0324207, 0.0006229287, 0.06410533), theta = c(3, 
        -0.005734982, -0.02541216, -0.004763663, -0.06341044)))
  sample_full_conditionals(dir, s, starts = s@supplement$truth, priors = priors)
}

#' @title Function \code{full_conditionals}
#' @description Check MCMC libraries against their full conditionals.
#' @export
#' @param priors names of prior distributions to try for phi, alpha, and delta.
full_conditionals = function(priors = special_beta_priors()){
  for(prior in priors){
    stopifnot(prior %in% special_beta_priors())
    full_conditionals_paschold(priors = prior)
    plot_full_conditionals(paste0(prior, "_full_conditionals_paschold/"))
    full_conditionals_simulated(priors = prior)
    plot_full_conditionals(paste0(prior, "_full_conditionals_simulated/"))
  }
}
