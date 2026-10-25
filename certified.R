.lr_root <- .learnref_root
.lr_lib <- file.path(.lr_root, paste0("certified_probability", .Platform$dynlib.ext))
.lr_native <- file.path(.lr_root, c("probability.cpp", "certificates.cpp"))
if (!file.exists(.lr_lib) || any(file.mtime(.lr_lib) < file.mtime(.lr_native))) {
  .lr_build <- function() {
    old <- Sys.getenv(c("PKG_CPPFLAGS", "PKG_LIBS"), unset=NA_character_)
    on.exit({
      Sys.unsetenv(names(old)[is.na(old)])
      if (any(!is.na(old))) do.call(Sys.setenv, as.list(old[!is.na(old)]))
    })
    inc <- Sys.getenv("MPFR_INCLUDE", if (file.exists("/opt/homebrew/include/mpfr.h")) "/opt/homebrew/include" else "")
    lib <- Sys.getenv("MPFR_LIB", if (file.exists("/opt/homebrew/lib/libmpfr.dylib")) "/opt/homebrew/lib" else "")
    cppflags <- paste(Sys.getenv("PKG_CPPFLAGS"), if (nzchar(inc)) paste0("-I", shQuote(inc)) else "")
    ldflags <- paste(Sys.getenv("PKG_LIBS"), if (nzchar(lib)) paste0("-L", shQuote(lib)) else "", "-lmpfr -lgmp")
    Sys.setenv(PKG_CPPFLAGS=cppflags, PKG_LIBS=ldflags)
    unlink(file.path(.lr_root, "probability.o"))
    status <- system2(file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", shQuote(file.path(.lr_root, "probability.cpp")), "-o", shQuote(.lr_lib)))
    if (status != 0L) stop("Native compilation failed. Install MPFR/GMP and an R-compatible C++ compiler.")
  }
  .lr_build()
  rm(.lr_build)
}
if (!is.loaded("lr_beta")) dyn.load(.lr_lib)

reset_gap_audit<-function(){gap_audit<<-new.env();gap_audit$n<-gap_audit$zero<-gap_audit$failed<-0L;gap_audit$max_gap<-0;gap_audit$max_guard<-0;gap_audit$max_iterations<-0L;gap_audit$max_bits<-0L;gap_audit$by_J<-integer(64)}

reset_gap_audit()

box_certificate<-function(lo,hi,a,tolerance=1e-7,max_bits=1024L){
 stopifnot(length(lo)==length(a),length(hi)==length(a),tolerance>0)
 z<-.Call('certified_box',as.double(lo),as.double(hi),as.double(a),as.double(tolerance),as.integer(max_bits))
 list(value=z[1],upper=z[2],gap=z[3],bits=as.integer(z[4]),iterations=as.integer(z[5]),failed=as.logical(z[6]),guard=0)
}

box_min_discrepancy<-function(lo,hi,a){
 z<-box_certificate(lo,hi,a);gap_audit$n<-gap_audit$n+1L;gap_audit$by_J[length(a)]<-gap_audit$by_J[length(a)]+1L
 gap_audit$zero<-gap_audit$zero+as.integer(z$upper==0);gap_audit$failed<-gap_audit$failed+as.integer(z$failed)
 if(is.finite(z$gap))gap_audit$max_gap<-max(gap_audit$max_gap,z$gap)
 gap_audit$max_iterations<-max(gap_audit$max_iterations,z$iterations);gap_audit$max_bits<-max(gap_audit$max_bits,z$bits);z$value
}

subset_min_discrepancy<-function(R,a,keep,low,high,evals=300L){
 T<-sum(a);q<-a/T;lo<-a+low;hi<-a+high;box<-box_min_discrepancy(lo,hi,a)
 l<-sum(lo);u<-sum(hi)
 for(i in seq_len(60)){mid<-(l+u)/2;if(sum(pmin(hi,pmax(lo,mid*q)))>=mid)l<-mid else u<-mid}
 z<-pmin(hi,pmax(lo,l*q));z<-z/sum(z);h0<-log(q/z);h0<-h0-h0[length(h0)]
 total_low<-T+sum(sort(rowSums(R))[seq_len(keep)]);best<- -Inf;best_h<-h0
 proposal<-function(h){hh<-c(h,0);sr<-subset_ratio_max(R,a,keep,exp(hh),total_low);v<-sum(a*hh)-T*log(sr$rho);if(is.finite(v)&&v>best){best<<-v;best_h<<-hh};v}
 proposal(h0[-length(h0)])
 if(length(a)==2)optimize(function(h)-proposal(h),h0[1]+c(-5,5),tol=1e-10) else if(length(a)>2)optim(h0[-length(h0)],function(h)-proposal(h),method='Nelder-Mead',control=list(maxit=evals,reltol=1e-10))
 certificate<-function(h){b<-exp(h);if(any(!is.finite(b))||any(b<=0))return(0);sr<-subset_ratio_max(R,a,keep,b,total_low);if(!is.finite(sr$rho)||sr$rho<=0)return(0);.Call('certified_subset',matrix(as.double(R),nrow(R)),as.double(a),as.integer(keep),as.double(b),as.double(sr$rho))}
 max(0,box,certificate(h0),certificate(best_h))
}

subset_ratio_max <- function(R, a, keep, eh, tot_low) {
  rs <- rowSums(R); ya <- as.vector(R %*% eh); base_num <- sum(eh * a); base_den <- sum(a)
  sel <- order(ya / rs, decreasing = TRUE)[seq_len(keep)]
  rho <- (base_num + sum(ya[sel])) / (base_den + sum(rs[sel]))
  for (it in 1:50) {
    sc <- ya - rho * rs
    sel <- order(sc, decreasing = TRUE)[seq_len(keep)]
    new <- (base_num + sum(ya[sel])) / (base_den + sum(rs[sel]))
    if (new <= rho * (1 + 1e-13)) break
    rho <- new
  }
  sc <- ya - rho * rs
  F <- base_num - rho * base_den + sum(sort(sc, decreasing = TRUE)[seq_len(keep)])
  list(rho = rho + max(0, F) / tot_low, sel = sel)
}

lr_op <- function(x,y=0,op=1,upper=TRUE) .Call('lr_arithmetic',as.double(x),as.double(y),as.integer(op))[,if(upper)2 else 1]

lr_add <- function(x,y,upper=TRUE) lr_op(x,y,1,upper)

lr_sub <- function(x,y,upper=TRUE) lr_op(x,y,2,upper)

lr_mul <- function(x,y,upper=TRUE) lr_op(x,y,3,upper)

lr_div <- function(x,y,upper=TRUE) lr_op(x,y,4,upper)

lr_log <- function(x,upper=TRUE) lr_op(x,0,5,upper)

lr_tail <- function(x,a) .Call('lr_psi',as.double(x),as.double(a))

lr_thresholds <- function(K,alpha) .Call('lr_by_thresholds',as.integer(K),as.double(alpha))

lr_select <- function(p,alpha=.05) {
  th<-lr_thresholds(length(p),alpha)[,1];o<-order(p);ii<-which(p[o]<=th)
  if(!length(ii))return(rep(FALSE,length(p)))
  p<=th[max(ii)]
}

lr_qbeta <- function(pl,pu,a,b,proposal=qbeta((pl+pu)/2,a,b)) {
  n<-max(length(pl),length(pu),length(a),length(b),length(proposal))
  .Call('lr_beta',as.integer(rep(a,length.out=n)),as.integer(rep(b,length.out=n)),
        as.double(rep(pl,length.out=n)),as.double(rep(pu,length.out=n)),as.double(rep(proposal,length.out=n)))
}

lr_cdf <- function(k,n,p,upper=FALSE) .Call('lr_binomial',as.integer(n),as.integer(k),as.double(p),as.logical(upper))

lr_limits <- function(y,m,N,cap,eta,cache) {
  key<-paste('limits',m,N,cap,sprintf('%a',eta));k0<-ceiling(N/2);d<-N-cap
  if(!exists(key,cache,inherits=FALSE)) {
    nl<-lr_qbeta(eta,eta,k0-cap,d-k0+cap+1)
    nr<-lr_qbeta(eta,eta,N-k0+1-cap,k0)
    assign(key,list(limits=matrix(NA_real_,m+1,4),nl=nl,nr=nr),cache)
  }
  obj<-get(key,cache,inherits=FALSE);idx<-unique(y[is.na(obj$limits[y+1,1])])
  if(length(idx)) {
    z<-matrix(NA_real_,length(idx),4)
    z[idx==0,1:2]<-0;z[idx==m,3:4]<-1
    j<-which(idx>0)
    if(length(j))z[j,1:2]<-lr_qbeta(obj$nr[1],obj$nr[2],idx[j],m-idx[j]+1)
    j<-which(idx<m)
    if(length(j)) {
      pl<-lr_sub(1,obj$nl[2],FALSE);pu<-lr_sub(1,obj$nl[1],TRUE)
      prop<-qbeta(mean(obj$nl),idx[j]+1,m-idx[j],lower.tail=FALSE)
      z[j,3:4]<-lr_qbeta(pl,pu,idx[j]+1,m-idx[j],prop)
    }
    obj$limits[idx+1,]<-z;assign(key,obj,cache)
  }
  obj$limits[y+1,,drop=FALSE]
}

lr_low_keep <- function(p,m,eps) {
  k<-as.integer(qbinom(eps,m,p))-1L
  
  
  repeat {
    bad<-which(k>=0 & lr_cdf(k,m,p)[,2]>eps)
    if(!length(bad))break
    k[bad]<-k[bad]-1L
  }
  cand<-k+1L;ok<-cand<m & lr_cdf(cand,m,p)[,2]<=eps;k[ok]<-cand[ok]
  pmax(0L,k+1L)
}

lr_screen <- function(y,m,N,cap,eta,eps,cache) {
  key<-paste('screen',m,N,cap,sprintf('%a',eta),sprintf('%a',eps))
  if(!exists(key,cache,inherits=FALSE))assign(key,matrix(NA_real_,m+1,2),cache)
  tab<-get(key,cache,inherits=FALSE);idx<-unique(y[is.na(tab[y+1,1])])
  if(length(idx)) {
    ci<-lr_limits(idx,m,N,cap,eta,cache)
    tab[idx+1,]<-cbind(lr_low_keep(ci[,1],m,eps),m-lr_low_keep(lr_sub(1,ci[,4],FALSE),m,eps))
    assign(key,tab,cache)
  }
  tab[y+1,,drop=FALSE]
}

lr_median_box <- function(ci,a,tol) {
  outer<-box_certificate(ci[,1],ci[,4],a,tol)
  
  
  if(any(ci[,2]>ci[,3]))return(c(value=0,gap=Inf,failed=1))
  inner<-box_certificate(ci[,2],ci[,3],a,tol)
  gap<-lr_sub(inner$upper,outer$value)
  fail<-outer$failed||inner$failed||!is.finite(gap)||gap>tol
  c(value=if(fail)0 else outer$value,gap=gap,failed=as.integer(fail))
}

lr_binary_limits <- function(ci) rbind(ci,c(lr_sub(1,ci[4],FALSE),lr_sub(1,ci[3],TRUE),lr_sub(1,ci[2],FALSE),lr_sub(1,ci[1],TRUE)))

lr_scan <- function(values,windows,cap,weights=NULL,total_tests=NULL,alpha=.05,
                    lambdas=c(.5,.1),cache=new.env(),tolerance=1e-7,
                    allocation=c('study','recommended'),mass_split=0) {
  allocation<-match.arg(allocation)
  stopifnot(length(mass_split)==1,mass_split>=0,mass_split<1)
  if(is.matrix(values))values<-list(values)
  C<-length(values);m<-nrow(values[[1]]);K<-ncol(values[[1]]);N<-K-1;d<-N-cap;k0<-ceiling(N/2)
  stopifnot(K>=3,cap==floor(cap),cap>=0,cap<N/2,m==floor(m),m*N<2^52,
            all(vapply(values,function(x)all(dim(x)==c(m,K))&&all(is.finite(x)),FALSE)),
            all(lambdas>0&lambdas<1),alpha>0,alpha<1,tolerance>0)
  if(is.null(total_tests))total_tests<-K
  stopifnot(total_tests>=K)
  windows<-as.matrix(windows);V<-nrow(windows);J<-ncol(windows)-1L
  ranks<-lapply(seq_len(V),function(v)ceiling(m*windows[v,]));stopifnot(all(vapply(ranks,function(r)min(r)>=1&&max(r)<=m&&all(diff(r)>0),FALSE)))
  if(is.null(weights))weights<-rep(1,V)
  stopifnot(length(weights)==V,all(is.finite(weights)&weights>0))
  
  denom<-0;for(w in weights)denom<-lr_add(denom,w)
  omega<-lr_div(weights,lr_mul(C,denom),FALSE)
  cuts<-lr_thresholds(total_tests,alpha);u1<-cuts[1,1]
  beta<-lr_div(lr_mul(u1,omega,FALSE),2,FALSE)
  eta<-1e-12;eps<-lr_div(1e-9,2*N*C*V*(J+2),FALSE)
  if(allocation=='recommended') {
    
    shape_share<-lr_sub(1,mass_split,FALSE)
    beta<-lr_mul(beta,shape_share,FALSE)
    
    screen_limit<-lr_div(lr_mul(lr_mul(lr_sub(1,max(lambdas),FALSE),shape_share,FALSE),u1,FALSE),4,FALSE)
    eta<-min(eta,lr_div(screen_limit,8*C*V*(J+2),FALSE))
    eps<-min(eps,lr_div(screen_limit,8*N*C*V*(J+2),FALSE))
  }
  allowance<-lr_mul(C*V*(J+2),lr_add(lr_mul(2,eta),lr_mul(2*N,eps)))
  allowance_lo<-lr_mul(C*V*(J+2),lr_add(lr_mul(2,eta,FALSE),lr_mul(2*N,eps,FALSE),FALSE),FALSE)
  cols<-c('median','trimmed',as.vector(rbind(paste0('combined_',sprintf('%03d',round(lambdas*100))),paste0('mass_',sprintf('%03d',round(lambdas*100))),paste0('joint_',sprintf('%03d',round(lambdas*100))))))
  pmed<-pmedlo<-rep(1,K);ptr<-ptrlo<-pmtr<-pmtrlo<-rep(1,K);pm<-pmlo<-pmm<-pmmlo<-matrix(1,K,length(lambdas))
  failures<-0L;maxgap<-0;removed<-matrix(0,K,C*V);gap_cases<-0L
  maxcut<-lr_mul(cuts[total_tests,2],lr_sub(1,min(lambdas),TRUE))
  ar<-sort(unique(unlist(ranks)));ix<-lapply(ranks,function(r)match(r,ar))
  for(ch in seq_len(C)) {
    vals<-apply(values[[ch]],2,sort);ends<-vals[ar,,drop=FALSE];cum<-array(0L,c(length(ar),K,K))
    for(l in seq_len(K))cum[,,l]<-matrix(findInterval(as.vector(ends),vals[,l]),length(ar),K)
    for(v in seq_len(V)) {
      a<-diff(ranks[[v]]);aa<-c(sum(a),m+1-sum(a));idx<-ix[[v]]
      bc<-lr_mul(lambdas,beta[v],FALSE)
      for(k in seq_len(K)) {
        cc<-cum[idx,k,-k,drop=FALSE];dim(cc)<-c(length(idx),N)
        Y<-t(cc[-1,,drop=FALSE]-cc[-length(idx),,drop=FALSE]);outside<-cbind(cc[1,],m-cc[nrow(cc),]);ys<-cbind(Y,outside)
        med<-apply(ys,2,function(y)sort(y,partial=k0)[k0]);scr<-lr_screen(med,m,N,cap,eta,eps,cache)
        bad<-rowSums(ys<rep(scr[,1],each=N)|ys>rep(scr[,2],each=N))>0
        removed[k,(ch-1)*V+v]<-sum(bad)
        for(j in seq_len(length(lambdas)+1L)) {
          bt<-if(j==1)beta[v] else bc[j-1]
          ci<-lr_limits(med[seq_len(J)],m,N,cap,lr_div(bt,2*J,FALSE),cache)
          z<-lr_median_box(ci,a,tolerance);failures<-failures+z['failed'];maxgap<-max(maxgap,z['gap']);gap_cases<-gap_cases+1L
          tail<-lr_tail(z['value'],a)
          pp<-min(1,lr_div(lr_add(tail[,2],bt),omega[v]));plo<-min(1,lr_div(lr_add(tail[,1],bt,FALSE),omega[v],FALSE))
          if(j==1){pmed[k]<-min(pmed[k],pp);pmedlo[k]<-min(pmedlo[k],plo)}
          else{pm[k,j-1]<-min(pm[k,j-1],pp);pmlo[k,j-1]<-min(pmlo[k,j-1],plo)}
        }
        mid<-sort(rowSums(Y),partial=k0)[k0]
        for(j in seq_along(lambdas)) {
          ci<-lr_limits(mid,m,N,cap,lr_div(bc[j],2,FALSE),cache)
          z<-lr_median_box(lr_binary_limits(as.numeric(ci)),aa,tolerance)
          failures<-failures+z['failed'];maxgap<-max(maxgap,z['gap']);gap_cases<-gap_cases+1L
          tail<-lr_tail(z['value'],aa)
          pmm[k,j]<-min(pmm[k,j],1,lr_div(lr_add(tail[,2],bc[j]),omega[v]))
          pmmlo[k,j]<-min(pmmlo[k,j],1,lr_div(lr_add(tail[,1],bc[j],FALSE),omega[v],FALSE))
        }
        if(sum(bad)>cap)next
        R<-Y[!bad,,drop=FALSE];nr<-nrow(R);t<-nr-d
        srt<-apply(R,2,sort);if(is.null(dim(srt)))srt<-matrix(srt,nrow=nr)
        low<-colSums(srt[seq_len(d),,drop=FALSE]);high<-colSums(srt[(nr-d+1):nr,,drop=FALSE])
        z<-box_certificate(a+low,a+high,a,tolerance);x<-z$value
        failures<-failures+z$failed;maxgap<-max(maxgap,z$gap);gap_cases<-gap_cases+1L
        if(t>0) {
          sel<-order(abs(rowSums(R)-median(rowSums(R))))[seq_len(d)];vec<-a+colSums(R[sel,,drop=FALSE])
          feasible<-box_certificate(vec,vec,a,tolerance)
          possible<-lr_div(lr_tail(feasible$upper,a)[,1],omega[v],FALSE)<=maxcut
          if(possible)x<-max(x,subset_min_discrepancy(R,a,d,low,high))
        }
        tail<-lr_tail(x,a)
        ptr[k]<-min(ptr[k],1,lr_div(tail[,2],omega[v]));ptrlo[k]<-min(ptrlo[k],1,lr_div(tail[,1],omega[v],FALSE))
        rr<-sort(rowSums(R));ll<-sum(head(rr,d));uu<-sum(tail(rr,d));nn<-d*m
        z<-box_certificate(aa+c(ll,nn-uu),aa+c(uu,nn-ll),aa,tolerance)
        failures<-failures+z$failed;maxgap<-max(maxgap,z$gap);gap_cases<-gap_cases+1L
        tail<-lr_tail(z$value,aa)
        pmtr[k]<-min(pmtr[k],1,lr_div(tail[,2],omega[v]));pmtrlo[k]<-min(pmtrlo[k],1,lr_div(tail[,1],omega[v],FALSE))
      }
    }
  }
  ptr<-pmin(1,lr_add(ptr,allowance));pmtr<-pmin(1,lr_add(pmtr,allowance))
  ptrlo<-pmin(1,lr_add(ptrlo,allowance_lo,FALSE));pmtrlo<-pmin(1,lr_add(pmtrlo,allowance_lo,FALSE))
  out<-cbind(median=pmed,trimmed=ptr)
  outlo<-cbind(median=pmedlo,trimmed=ptrlo)
  for(j in seq_along(lambdas)) {
    
    sh<-pmin(1,lr_div(ptr,lr_sub(1,lambdas[j],FALSE)),lr_div(pm[,j],lambdas[j]))
    ma<-pmin(1,lr_div(pmtr,lr_sub(1,lambdas[j],FALSE)),lr_div(pmm[,j],lambdas[j]))
    out<-cbind(out,sh,ma,pmin(1,lr_mul(sh,2),lr_mul(ma,2)))
    shlo<-pmin(1,lr_div(ptrlo,lr_sub(1,lambdas[j],TRUE),FALSE),lr_div(pmlo[,j],lambdas[j],FALSE))
    malo<-pmin(1,lr_div(pmtrlo,lr_sub(1,lambdas[j],TRUE),FALSE),lr_div(pmmlo[,j],lambdas[j],FALSE))
    outlo<-cbind(outlo,shlo,malo,pmin(1,lr_mul(shlo,2,FALSE),lr_mul(malo,2,FALSE)))
  }
  colnames(out)<-cols
  probability_gap<-matrix(lr_sub(out,outlo),nrow=K)
  cutoff_gap<-max(lr_sub(cuts[,2],cuts[,1]))
  numerical_budget<-lr_mul(.01,u1,FALSE)
  imprecise<-probability_gap>numerical_budget
  out[imprecise]<-1
  if(cutoff_gap>numerical_budget)out[]<-1
  list(p=out,log_p=matrix(lr_log(out),nrow=K,dimnames=dimnames(out)),removed=removed,
       diagnostics=list(failures=unname(failures),max_discrepancy_gap=unname(maxgap),box_cases=gap_cases,
                        probability_bits=192L,allocation=allocation,mass_split=mass_split,
                        coverage=beta,screen_eta=eta,screen_epsilon=eps,
                        screen_allowance=allowance,min_BY_threshold=cuts[1,1],
                        max_probability_gap=max(probability_gap),probability_failures=sum(imprecise),
                        numerical_budget=numerical_budget,max_BY_rounding_gap=cutoff_gap))
}

