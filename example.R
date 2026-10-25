.example_args <- grep("^--file=", commandArgs(FALSE), value=TRUE)
.example_sources <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
.example_file <- if (length(.example_sources)) tail(.example_sources, 1L)[[1L]] else sub("^--file=", "", .example_args[[1L]])
source(file.path(dirname(normalizePath(.example_file)), "learnref.R"), local=TRUE)
rm(.example_args, .example_sources, .example_file)

set.seed(20261018)
K <- 40
m <- 2048
groups <- lapply(seq_len(K), function(k) {
  u <- runif(m)
  if (k <= 4) {
    i <- u > .4 & u < .6
    u[i] <- .4 + .2 * ((u[i] - .4) / .2)^2
  }
  qexp(u)
})
fit <- learnref_test(groups, cap=8)
stopifnot(identical(which(fit$selected), 1:4),
  attr(fit, "numerical_diagnostics")$failures == 0,
  attr(fit, "numerical_diagnostics")$probability_failures == 0)
print(fit)

gate_counts <- t(vapply(groups, function(x) c(low=sum(x <= 1), high=sum(x > 1)), numeric(2)))
gate_intervals <- learnref_gate_intervals(gate_counts, m=m, cap=8)
print(gate_intervals)
print(learnref_gate_ratio(gate_intervals, "low", "high"))
