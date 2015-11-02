#' @title Function \code{inits}
#' @description preparation for the \code{check_*} functions
#' @export
#' @return a list of things the \code{check_*} functions need
#' @param chain a \code{fbseq::Chain} object.
inits = function(chain){
  s = Starts(chain)
  list(chain = chain,
        beta = matrix(s@beta, nrow = chain@G),
        delta = matrix(s@delta, nrow = chain@G),
        epsilon = matrix(s@epsilon, nrow = chain@G),
        xi = matrix(s@xi, ncol = chain@L),
        flat = as.matrix(flatten(chain)),
        G = chain@G,
        design = matrix(chain@design, ncol = chain@L),
        N = chain@N,
        L = chain@L,
        name = names(which(chain@parameter_sets_update == 1)),
        s = s,
        y = matrix(chain@counts, ncol = chain@N))
}

#' @title Function \code{plotfc}
#' @description Plots samples from a full conditional: a histogram overlayed with the true 
#' target, plus a traceplot.
#' @export
#' @param x vector of MCMC samples
#' @param lkern log kernel of the true full conditional
#' @param name name of the parameter
#' @param postmean posterior mean
#' @param postmeansq posterior mean of squares
plotfc = function(x, lkern, name, postmean = NULL, postmeansq = NULL){
  if(!is.null(postmean) & !is.null(postmeansq)) postmeansq = postmeansq * sign(postmean)
  x = as.numeric(x)
  if(length(unique(x)) > 3) x = x[x > quantile(x, 0.0025) & x < quantile(x, 0.9975)]
  m = mean(x)
  xs = seq(min(x), max(x), length.out = 1e3)
  area = trapz(xs, exp(lkern(xs) - lkern(m)))
  plotfn = function(x){exp(lkern(x) - lkern(m))/area}
  df = data.frame(x = x)
  print(paste("plotting", name))

  pl = ggplot(data = df) + 
    geom_histogram(mapping = aes(x = x, y = ..density..), binwidth = diff(range(xs))/30, alpha = 0.6) +
    stat_function(fun = plotfn) + 
    labs(title = paste0(name, "\n")) + 
    xlab("\nMCMC samples") +
    ylab("Density\n") +
    theme(panel.background = element_rect(fill='white'),
               panel.border = element_rect(color="black", fill = NA),
               panel.grid.major = element_line(color="lightgray"),
               panel.grid.minor = element_blank())

  if(!is.null(postmean)) pl = pl + geom_vline(xintercept = postmean, color = "blue")
  if(!is.null(postmeansq)) pl = pl + geom_vline(xintercept = postmeansq, color = "green", linetype = 6)
  ggsave(filename = paste0("hist_", name, ".pdf"), plot = pl, width = 8, height = 8, dpi = 1600)


  t = floor(seq(from = 1, to = length(x), length.out = 1e3))
  y = x[t]

  pl = qplot(x = t, y = y, geom = "line") + 
    labs(title = paste0(name, "\n")) + 
    xlab("\nMCMC iteration") +
    ylab("MCMC sample\n") +
    theme(panel.background = element_rect(fill='white'),
               panel.border = element_rect(color="black", fill = NA),
               panel.grid.major = element_line(color="lightgray"),
               panel.grid.minor = element_blank())

  if(!is.null(postmean)) pl = pl + geom_hline(yintercept = postmean, color = "blue")
  if(!is.null(postmeansq)) pl = pl + geom_hline(yintercept = postmeansq, color = "green", linetype = 6)
  ggsave(filename = paste0("trace_", name, ".pdf"), plot = pl, width = 8, height = 8, dpi = 1600)
}

#' @title Function \code{ldig}
#' @description log inverse gamma density
#' @export
#' @return log inverse gamma density
#' @param x value
#' @param a shape
#' @param b scale
ldig = function(x, a, b){
  a*log(b) - lgamma(a) - (a + 1)*log(x) - b/x
}

#' @title Function \code{make_dirs}
#' @description Makes directories for checking MCMC samples against their full conditionals.
#' @export
#' @param dir directory for output files
make_dirs = function(dir){
  for(d in paste0(dir, c("", "chains", "plots"))) 
    if(!file.exists(d))
      dir.create(d)
}
