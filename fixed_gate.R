learnref_gate_intervals <- function(counts, m, cap, error=.05,
                                    pilot_error=0, cache=new.env()) {
  counts <- as.matrix(counts)
  K <- nrow(counts); H <- ncol(counts)
  stopifnot(K >= 1L, H >= 1L, length(m)==1L, is.finite(m),
            m >= 1, m<=.Machine$integer.max, m==floor(m), all(is.finite(counts)),
            all(counts==floor(counts)), all(counts>=0), all(counts<=m),
            length(error)==1L, is.finite(error), error>0, error<1,
            length(pilot_error)==1L, is.finite(pilot_error),
            pilot_error>=0, pilot_error<=1)
  kg <- ceiling(K/2)
  med <- apply(counts,2L,function(x) sort.int(x,partial=kg)[[kg]])
  gates <- colnames(counts)
  if (is.null(gates)) gates <- paste0('gate_',seq_len(H))
  stopifnot(!anyDuplicated(gates))
  ans <- data.frame(gate=gates,median_count=as.integer(med),
                    lower=rep(0,H),upper=rep(1,H),stringsAsFactors=FALSE)
  level <- lr_sub(error,pilot_error,FALSE)
  admissible <- length(cap)==1L && is.finite(cap) && cap==floor(cap) &&
                cap>=0 && cap<K/2 && level>0
  status <- if (admissible) 'certified' else 'full_interval'
  eta <- NA_real_
  if (admissible) {
    eta <- lr_div(level,2*H,FALSE)
    ci <- tryCatch(lr_limits(as.integer(med),as.integer(m),as.integer(K),
                             as.integer(cap),eta,cache),error=function(e) NULL)
    if (is.null(ci)) status <- 'full_interval' else {
      ok <- is.finite(ci[,1L]) & is.finite(ci[,4L]) &
            ci[,1L]>=0 & ci[,4L]<=1 & ci[,1L]<=ci[,4L]
      ans$lower[ok] <- ci[ok,1L]
      ans$upper[ok] <- ci[ok,4L]
      if (!all(ok)) status <- 'partial_full_interval'
    }
  }
  attr(ans,'calibration') <- list(K=K,m=m,cap=cap,gates=H,
      total_error=error,pilot_error=pilot_error,count_error=max(0,level),
      error_per_tail=eta,status=status,
      denominator='All original m observations in every group')
  ans
}

learnref_gate_ratio <- function(intervals,numerator,complement) {
  locate <- function(which_gate) {
    if (is.character(which_gate)) return(match(which_gate,intervals$gate))
    if (!(is.integer(which_gate) || is.double(which_gate)) ||
        length(which_gate)!=1L || !is.finite(which_gate) ||
        which_gate!=floor(which_gate))
      stop('Gate indices must be integer-valued numbers or gate names.')
    which_gate
  }
  g <- locate(numerator); b <- locate(complement)
  stopifnot(length(g)==1L,length(b)==1L,!is.na(g),!is.na(b),
            g>=1L,g<=nrow(intervals),b>=1L,b<=nrow(intervals),g!=b)
  gl <- intervals$lower[[g]]; gu <- intervals$upper[[g]]
  bl <- intervals$lower[[b]]; bu <- intervals$upper[[b]]
  if (!all(is.finite(c(gl,gu,bl,bu))) || gl<0 || bl<0 ||
      gu>1 || bu>1 || gl>gu || bl>bu) return(c(lower=0,upper=1))
  if (lr_add(gl,bl,FALSE)>1) return(c(lower=0,upper=1))
  den_lower_endpoint <- lr_add(gl,bu,TRUE)
  den_upper_endpoint <- lr_add(gu,bl,FALSE)
  lower <- if (den_lower_endpoint>0) lr_div(gl,den_lower_endpoint,FALSE) else 0
  upper <- if (den_upper_endpoint>0) lr_div(gu,den_upper_endpoint,TRUE) else 1
  c(lower=max(0,min(1,lower)),upper=max(0,min(1,upper)))
}

