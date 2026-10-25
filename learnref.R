.learnref_sources <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
if (!length(.learnref_sources)) stop("Load this file with source().")
.learnref_root <- dirname(normalizePath(tail(.learnref_sources, 1L)[[1L]], mustWork=TRUE))
rm(.learnref_sources)

source(file.path(.learnref_root, "certified.R"), local=TRUE)
source(file.path(.learnref_root, "fixed_gate.R"), local=TRUE)

learnref_cache <- new.env()

learnref_default_windows <- function() {
  g <- expand.grid(center = c(.3, .5, .7), width = c(.2, .4))
  t(apply(g, 1, function(z) c(z[1] - z[2] / 2, z[1], z[1] + z[2] / 2)))
}

learnref_test <- function(groups, cap, m = NULL, windows = learnref_default_windows(),
                          rule = c('combined', 'joint', 'median', 'trimmed'), alpha = .05, total_tests = NULL,
                          allocation=c('recommended','study')) {
  rule <- match.arg(rule)
  allocation <- match.arg(allocation)
  if (is.list(groups)) {
    mats <- lapply(groups, as.matrix); m <- nrow(mats[[1]])
    stopifnot(all(vapply(mats, nrow, 0) == m), length(unique(vapply(mats, ncol, 0))) == 1)
    x <- do.call(rbind, mats)
  } else { x <- as.matrix(groups); stopifnot(!is.null(m)) }
  K <- nrow(x) %/% m
  stopifnot(K * m == nrow(x), K >= 3, cap >= 0, cap < ceiling((K - 1) / 2), all(is.finite(x)))
  if (is.null(total_tests)) total_tests <- K
  values <- lapply(seq_len(ncol(x)),function(j)matrix(x[,j],m,K))
  fit <- lr_scan(values,windows,cap,total_tests=total_tests,alpha=alpha,lambdas=.5,cache=learnref_cache,
                 allocation=allocation,mass_split=if(rule=='joint').5 else 0)
  col <- c(combined = 'combined_050', joint = 'joint_050', median = 'median', trimmed = 'trimmed')[[rule]]
  out <- data.frame(group = seq_len(K), log_p = fit$log_p[, col], p = fit$p[, col])
  if (total_tests == K) out$selected <- lr_select(out$p, alpha)
  attr(out,'numerical_diagnostics') <- fit$diagnostics
  out
}

