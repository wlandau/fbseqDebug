#' @include full_conditionals_utils.R
NULL

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
beta_check = function(chain){
  Z = inits(chain)
  l = chain@effects_update
  name = paste0("beta_", l)
  bpm = matrix(chain@betaPostMean, nrow = chain@G)
  bpmsq = matrix(chain@betaPostMeanSquare, nrow = chain@G)

  for(v in colnames(Z$flat)[grep(name, colnames(Z$flat))]){
    x = as.numeric(Z$flat[,v])
    g = as.integer(strsplit(gsub(paste0(name, "_"), "", v), split = "_")[[1]])

    A = sum(Z$y[g, ]*Z$design[,l])
    B = 1/(2*Z$s@sigmaSquared[l] * Z$xi[g, l])
    C = Z$s@theta[l]

    lkern = Vectorize(function(x){
      A*x - B* (x - C)^2 - sum(exp(Z$design[,l] * x) * exp( Z$epsilon[g, ] + Z$design[,-l] %*% Z$beta[g, -l]))
    }, "x")

    plotfc(x, lkern, v, bpm[g, l], sqrt(bpmsq[g, l]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
epsilon_check = function(chain){
  Z = inits(chain)
  attach(Z, warn.conflicts = F)
  epm = matrix(chain@epsilonPostMean, nrow = G)
  epmsq = matrix(chain@epsilonPostMeanSquare, nrow = G)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    ind = as.integer(strsplit(gsub(paste0(name, "_"), "", v), split = "_")[[1]])
    n = ind[1]
    g = ind[2]

    A = y[g, n]
    B = 1/(2*s@rho[n]*s@gamma[g])
    C = chain@h[n]
    D = exp(sum(design[n,] * Z$beta[g,]))

    lkern = function(x){A*x - B*(x-C)^2 - D*exp(x)}
    plotfc(x, lkern, v, epm[g, n], sqrt(epmsq[g, n]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
gamma_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    g = as.integer(gsub(paste0(name, "_"), "", v))

    shape = (N + s@nuGamma[1])/2
    scale = (s@nuGamma[1]*s@tauGamma[1] + sum((epsilon[g,] - chain@h)^2/s@rho))/2

    lkern = function(x){ldig(x, shape, scale)}
    plotfc(x, lkern, v, chain@gammaPostMean[g], sqrt(chain@gammaPostMeanSquare[g]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
nuGamma_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  lkern = function(x){
    G*(-lgamma(x/2) + (x/2) * log(x*s@tauGamma[1]/2)) - 
    (x/2) * sum(log(s@gamma) + s@tauGamma[1]/s@gamma)
  }
  plotfc(as.numeric(flat), lkern, name, chain@nuGammaPostMean, sqrt(chain@nuGammaPostMeanSquare))
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
nuRho_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  lkern = function(x){
    N*(-lgamma(x/2) + (x/2) * log(x*s@tauRho[1]/2)) - 
    (x/2) * sum(log(s@rho) + s@tauRho[1]/s@rho)
  }
  plotfc(as.numeric(flat), lkern, name, chain@nuRhoPostMean, sqrt(chain@nuRhoPostMeanSquare))
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
rho_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    n = as.integer(gsub(paste0(name, "_"), "", v))
   
    shape = (G + s@nuRho[1])/2
    scale = (s@nuRho[1]*s@tauRho[1] + sum((epsilon[,n] - chain@h[n])^2/s@gamma))/2

    lkern = function(x){ldig(x, shape, scale)}
    plotfc(x, lkern, v, chain@rhoPostMean[n], sqrt(chain@rhoPostMeanSquare[n]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
sigmaSquared_check = function(chain){
  Z = inits(chain)
  attach(Z, warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    l = as.integer(gsub(paste0(name, "_"), "", v))
    shape = (G - 1)/2
    scale = 0.5 * sum((Z$beta[,l] - s@theta[l])^2/xi[,l])
    lkern = function(x){ldig(x, shape, scale)}
    plotfc(x, lkern, v, chain@sigmaSquaredPostMean[l], sqrt(chain@sigmaSquaredPostMeanSquare[l]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
tauGamma_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  shape = s@aGamma + G*s@nuGamma/2
  rate = s@bGamma + (s@nuGamma/2) * sum(1/s@gamma)
  lkern = function(x){dgamma(x, shape = shape, rate = rate, log = T)}
  plotfc(as.numeric(flat), lkern, name, chain@tauGammaPostMean, sqrt(chain@tauGammaPostMeanSquare))
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
tauRho_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  shape = s@aRho + N*s@nuRho/2
  rate = s@bRho + (s@nuRho/2) * sum(1/s@rho)
  lkern = function(x){dgamma(x, shape = shape, rate = rate, log = T)}
  plotfc(as.numeric(flat), lkern, name, chain@tauRhoPostMean, sqrt(chain@tauRhoPostMeanSquare))
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
theta_check = function(chain){
  Z = inits(chain)
  attach(Z, warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    l = as.integer(gsub(paste0(name, "_"), "", v))
    A = (1/(s@c[l]^2) + (1/s@sigmaSquared[l]) * sum(1/xi[,l]))/2
    B = (1/s@sigmaSquared[l]) * sum(Z$beta[,l]/xi[,l])
    lkern = function(x){dnorm(x, mean = B/(2*A), sd = sqrt(1/(2*A)), log = T)}
    plotfc(x, lkern, v, chain@thetaPostMean[l], sqrt(chain@thetaPostMeanSquare[l]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
xi_check = function(chain){
  Z = inits(chain)
  attach(Z, warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    ind = as.integer(strsplit(gsub(paste0(name, "_"), "", v), split = "_")[[1]])
    l = ind[1]
    g = ind[2]
   
    xipm = matrix(chain@xiPostMean, nrow = chain@G)
    xipmsq = matrix(chain@xiPostMeanSquare, nrow = chain@G)

    prior = alternate_priors()[chain@priors[l]]
    z = (Z$beta[g, l] - s@theta[l])^2/(2 * s@sigmaSquared[l])
    if(prior == "Laplace") {
      a = z
      b = s@k[l]
      lkern = function(x){-0.5*log(x) - a/x - b*x}
    } else if (prior == "t") {
      a = s@k[l] + 1.5
      b = z + s@r[l]
      lkern = function(x){-a*log(x) - b/x}
    } else if (prior == "horseshoe") {
      lkern = function(x){-log(x*(1 + x)) - z/x}
    } else {
      lkern = function(x){ifelse(x > 0.9 & x < 1.1, 1, 0)}
    }

    plotfc(x, lkern, v, xipm[g, l], sqrt(xipmsq[g, l]))
  }
}
