#' @include full_conditionals_utils.R
NULL

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
alp_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    g = as.integer(gsub(paste0(name, "_"), "", v))
    
    A = sum(y[g , group != 2]) - sum(y[g , group == 2])
    B = s@theAlp
    C = 2*s@sigAlp^2 * s@xiAlp[g]
    D = exp(s@phi[g] - s@del[g]) * sum(exp(eps[g, group == 1])) + 
          exp(s@phi[g] + s@del[g]) * sum(exp(eps[g, group == 3]))
    E = exp(s@phi[g] + s@del[g]) * sum(exp(eps[g, group == 2]))

    lkern = function(x){A*x - (x - B)^2/C - D*exp(x) - E*exp(-x)}
    plotfc(x, lkern, v, chain@alpPostMean[g], sqrt(chain@alpPostMeanSq[g]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
del_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    g = as.integer(gsub(paste0(name, "_"), "", v))
    
    A = sum(y[g , group != 1]) - sum(y[g , group == 1])
    B = s@theDel
    C = 2*s@sigDel^2 * s@xiDel[g]
    D = exp(s@phi[g] - s@alp[g]) * sum(exp(eps[g, group == 2])) + 
          exp(s@phi[g] + s@alp[g]) * sum(exp(eps[g, group == 3]))
    E = exp(s@phi[g] + s@alp[g]) * sum(exp(eps[g, group == 1]))

    lkern = function(x){A*x - (x - B)^2/C - D*exp(x) - E*exp(-x)}
    plotfc(x, lkern, v, chain@delPostMean[g], sqrt(chain@delPostMeanSq[g]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
eps_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  epm = matrix(chain@epsPostMean, nrow = G)
  epms = matrix(chain@epsPostMeanSq, nrow = G)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    ind = as.integer(strsplit(gsub(paste0(name, "_"), "", v), split = "_")[[1]])
    n = ind[1]
    g = ind[2]

    A = y[g, n]
    B = 0
    C = 2*s@rho[n]^2*s@gam[g]^2
    D = exp(eta(n, g, s, group))
    E = 0

    lkern = function(x){A*x - (x - B)^2/C - D*exp(x) - E*exp(-x)}
    plotfc(x, lkern, v, epm[g, n], sqrt(epms[g, n]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
gam_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    g = as.integer(gsub(paste0(name, "_"), "", v))
   
    x = x^2
    shape = (N + s@nuGam[1])/2
    scale = (s@nuGam[1]*s@tauGam[1]^2 + sum((eps[g,]/s@rho)^2))/2

    lkern = function(x){ldig(x, shape, scale)}
    plotfc(x, lkern, v, chain@gamPostMean[g]^2, chain@gamPostMeanSq[g])
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
nuGam_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  lkern = function(x){G*(
    -lgamma(x/2) + (x/2) * log((x*s@tauGam[1]^2)/2) - (x/G) * sum(log(s@gam) + (s@tauGam[1]^2/2) * (1/s@gam^2))
  )}
  plotfc(as.numeric(flat), lkern, name, chain@nuGamPostMean, sqrt(chain@nuGamPostMeanSq))
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
nuRho_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  lkern = function(x){N*(
    -lgamma(x/2) + (x/2) * log((x*s@tauRho[1]^2)/2) - (x/N) * sum(log(s@rho) + (s@tauRho[1]^2/2) * (1/s@rho^2))
  )}
  plotfc(as.numeric(flat), lkern, name, chain@nuRhoPostMean, sqrt(chain@nuRhoPostMeanSq))
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
phi_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    g = as.integer(gsub(paste0(name, "_"), "", v))
    
    A = N*mean(y[g,])
    B = s@thePhi
    C = 2*s@sigPhi^2 * s@xiPhi[g]
    D = exp(s@alp[g] - s@del[g]) * sum(exp(eps[g, group == 1])) +
          exp(- s@alp[g] + s@del[g]) * sum(exp(eps[g, group == 2])) +
          exp(s@alp[g] + s@del[g]) * sum(exp(eps[g, group == 3]))
    E = 0

    lkern = function(x){A*x - (x - B)^2/C - D*exp(x) - E*exp(-x)}
    plotfc(x, lkern, v, chain@phiPostMean[g], sqrt(chain@phiPostMeanSq[g]))
  }
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
   
    x = x^2
    shape = (G + s@nuRho[1])/2
    scale = (s@nuRho[1]*s@tauRho[1]^2 + sum((eps[,n]/s@gam)^2))/2

    lkern = function(x){ldig(x, shape, scale)}
    plotfc(x, lkern, v, chain@rhoPostMean[n]^2, chain@rhoPostMeanSq[n])
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
sigAlp_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  shape = (G - 1)/2
  scale = 0.5 * sum((s@alp - s@theAlp)^2/s@xiAlp)
  lkern = function(x){ldig(x, shape, scale)}
  plotfc(as.numeric(flat)^2, lkern, name, chain@sigAlpPostMean^2, chain@sigAlpPostMeanSq)
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
sigDel_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  shape = (G - 1)/2
  scale = 0.5 * sum((s@del - s@theDel)^2/s@xiDel)
  lkern = function(x){ldig(x, shape, scale)}
  plotfc(as.numeric(flat)^2, lkern, name, chain@sigDelPostMean^2, chain@sigDelPostMeanSq)
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
sigPhi_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  shape = (G - 1)/2
  scale = 0.5 * sum((s@phi - s@thePhi)^2/s@xiPhi)
  lkern = function(x){ldig(x, shape, scale)}
  plotfc(as.numeric(flat)^2, lkern, name, chain@sigPhiPostMean^2, chain@sigPhiPostMeanSq)
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
tauGam_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  shape = s@aGam + G*s@nuGam/2
  rate = s@bGam + (s@nuGam/2) * sum(1/s@gam^2)
  lkern = function(x){dgamma(x, shape = shape, rate = rate, log = T)}
  plotfc(as.numeric(flat)^2, lkern, name, chain@tauGamPostMean^2, chain@tauGamPostMeanSq)
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
tauRho_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  shape = s@aRho + N*s@nuRho/2
  rate = s@bRho + (s@nuRho/2) * sum(1/s@rho^2)
  lkern = function(x){dgamma(x, shape = shape, rate = rate, log = T)}
  plotfc(as.numeric(flat)^2, lkern, name, chain@tauRhoPostMean^2, chain@tauRhoPostMeanSq)
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
theAlp_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  x = as.numeric(flat)
  denom = s@cAlp^2*sum(1/s@xiAlp) + s@sigAlp^2
  numm = s@cAlp^2 * sum(s@alp/s@xiAlp)
  numv = s@cAlp^2 * s@sigAlp^2
  lkern = function(x){dnorm(x, mean = numm/denom, sd = sqrt(numv/denom), log = T)}
  plotfc(as.numeric(flat), lkern, name, chain@theAlpPostMean, sqrt(chain@theAlpPostMeanSq))
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
theDel_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  x = as.numeric(flat)
  denom = s@cDel^2*sum(1/s@xiDel) + s@sigDel^2
  numm = s@cDel^2 * sum(s@del/s@xiDel)
  numv = s@cDel^2 * s@sigDel^2
  lkern = function(x){dnorm(x, mean = numm/denom, sd = sqrt(numv/denom), log = T)}
  plotfc(as.numeric(flat), lkern, name, chain@theDelPostMean, sqrt(chain@theDelPostMeanSq))
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
thePhi_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  x = as.numeric(flat)
  denom = s@cPhi^2*sum(1/s@xiPhi) + s@sigPhi^2
  numm = s@cPhi^2 * sum(s@phi/s@xiPhi)
  numv = s@cPhi^2 * s@sigPhi^2
  lkern = function(x){dnorm(x, mean = numm/denom, sd = sqrt(numv/denom), log = T)}
  plotfc(as.numeric(flat), lkern, name, chain@thePhiPostMean, sqrt(chain@thePhiPostMeanSq))
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
xiPhi_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    g = as.integer(gsub(paste0(name, "_"), "", v))

    prior = alternate_priors()[chain@phiPrior]
    z = (s@phi[g] - s@thePhi)^2/(2 * s@sigPhi^2)
    if(prior == "Laplace") {
      a = z
      b = s@kPhi
      lkern = function(x){-0.5*log(x) - a/x - b*x}
    } else if (prior == "t") {
      a = s@kPhi + 1.5
      b = z + s@rPhi
      lkern = function(x){-a*log(x) - b/x}
    } else if (prior == "horseshoe") {
      lkern = function(x){-log(x*(1 + x)) - z/x}
    } else {
      lkern = function(x){ifelse(x > 0.9 & x < 1.1, 1, 0)}
    }
    
    plotfc(x, lkern, v, chain@xiPhiPostMean[g], sqrt(chain@xiPhiPostMeanSq[g]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
xiAlp_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    g = as.integer(gsub(paste0(name, "_"), "", v))
   
    prior = alternate_priors()[chain@phiPrior]
    z = (s@alp[g] - s@theAlp)^2/(2 * s@sigAlp^2)
    if(prior == "Laplace") {
      a = z
      b = s@kAlp
      lkern = function(x){-0.5*log(x) - a/x - b*x}
    } else if (prior == "t") {
      a = s@kAlp + 1.5
      b = z + s@rAlp
      lkern = function(x){-a*log(x) - b/x}
    } else if (prior == "horseshoe") {
      lkern = function(x){-log(x*(1 + x)) - z/x}
    } else {
      lkern = function(x){ifelse(x > 0.9 & x < 1.1, 1, 0)}
    }

    plotfc(x, lkern, v, chain@xiAlpPostMean[g], sqrt(chain@xiAlpPostMeanSq[g]))
  }
}

#' @title \code{*_check} functions
#' @description Check MCMC parameter samples against the true full conditional.
#' @export
#' @param chain a \code{fbseq::Chain} object
xiDel_check = function(chain){
  attach(inits(chain), warn.conflicts = F)
  for(v in colnames(flat)){
    x = as.numeric(flat[,v])
    g = as.integer(gsub(paste0(name, "_"), "", v))
   
    prior = alternate_priors()[chain@phiPrior]
    z = (s@del[g] - s@theDel)^2/(2 * s@sigDel^2)
    if(prior == "Laplace") {
      a = z
      b = s@kDel
      lkern = function(x){-0.5*log(x) - a/x - b*x}
    } else if (prior == "t") {
      a = s@kDel + 1.5
      b = z + s@rDel
      lkern = function(x){-a*log(x) - b/x}
    } else if (prior == "horseshoe") {
      lkern = function(x){-log(x*(1 + x)) - z/x}
    } else {
      lkern = function(x){ifelse(x > 0.9 & x < 1.1, 1, 0)}
    }

    plotfc(x, lkern, v, chain@xiDelPostMean[g], sqrt(chain@xiDelPostMeanSq[g]))
  }
}
