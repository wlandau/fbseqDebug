#' @title Function \code{inits}
#' @description preparation for the \code{check_*} functions
#' @export
#' @return a list of things the \code{check_*} functions need
#' @param chain a \code{fbseq::Chain} object.
inits = function(chain){
  group = chain@group
  s = Starts(chain)
  y = matrix(chain@counts, ncol = length(group))

  list(chain = chain,
        eps = matrix(s@eps, nrow = dim(y)[1]),
        flat = as.matrix(flatten(chain)),
        G = dim(y)[1],
        group = group,
        N = length(group),
        name = names(which(chain@updates == 1)),
        s = s,
        y = y)
}

#' @title Function \code{eta}
#' @description Phi_g - alpha_g, phi_g + delta_g, or phi_g + alpha_g, depending on the library \code{n}. 
#' @export
#' @return Phi_g - alpha_g, phi_g + delta_g, or phi_g + alpha_g, depending on the library \code{n}. 
#' @param n library
#' @param g gene
#' @param s a \code{Starts} object
#' @param group Experimental design. A vector of integers,
#' one for each RNA-seq sample/library, denoting the genetic
#' variety of that sample. You must use 1 for parent 1, 2 for parent 2,
#' and 3 for the hybrid.
eta = function(n, g, s, group){
  if(group[n] == 1) return(s@phi[g] + s@alp[g] - s@del[g])
  if(group[n] == 2) return(s@phi[g] - s@alp[g] + s@del[g])
  if(group[n] == 3) return(s@phi[g] + s@alp[g] + s@del[g])
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
  x = x[x > quantile(x, 0.0025) & x < quantile(x, 0.9975)]
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
