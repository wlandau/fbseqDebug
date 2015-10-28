#' @title Function \code{compare_points}
#' @description Compare MCMC parameter samples of a single \code{Chain}
#' with a named vector of point estimates. This is useful, for example, in
#' simulation studies that begin with a set of "true" parameters. We can see if
#' our 95\% credible intervals recapture the truth.
#' @export
#' 
#' @return A data frame with parameters as rows and the following columns.
#' \describe{
#'   \item{\code{pointsIn95ci}}{a matrix of logical values, each \code{TRUE}
#' if and only if the parameter's starting value lies within the 95-percent credible interval
#' estimated from the MCMC samples.}
#'   \item{\code{meansMinuspoints}}{mean of each parameter's MCMC samples (excluding burnin and thinning)
#' minus that parameter's starting value.}
#'   \item{\code{mediansMinuspoints}}{The same as \code{meansMinuspoints}, but with medians}
#' }
#' @param chain \code{Chain} object
#' @param points either a \code{Starts} object or a named numeric vector with 
#' names in \code{names(mcmc(chain))}. \code{Starts} objects will be automatically
#' flattened into named numeric vectors.
compare_points = function(chain, points){

  if(class(points) == "Starts")
    points = flatten(points)

  if(is.null(names(points)))
    stop("names(points) must be nonempty and comprise elements of colnames(mcmc(chain)).")

  parms = flatten(chain)
  ns = intersect(colnames(parms), names(points))
  parms = as.matrix(parms[,ns], ncol = length(ns))
  points = points[ns]
  
  fs = list(
    means = mean,
    medians = median,
    q025 = function(x) quantile(x, 0.025),
    q975 = function(x) quantile(x, 0.975))

  s = lapply(fs, function(f) apply(parms, 2, f))

  data.frame(pointsIn95ci = s$q025 < points & points < s$q975,
                    meansMinusPoints = s$means - points,
                    mediansMinusPoints = s$medians - points)
}
