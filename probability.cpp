
#include "certificates.cpp"
#include <limits>

static Interval expi(const Interval& a){Interval z;mpfr_exp(z.lo,a.lo,MPFR_RNDD);mpfr_exp(z.hi,a.hi,MPFR_RNDU);return z;}
static Interval sqrti(const Interval& a){Interval z;mpfr_sqrt(z.lo,a.lo,MPFR_RNDD);mpfr_sqrt(z.hi,a.hi,MPFR_RNDU);return z;}
static Interval unit(Interval z){
  if(mpfr_sgn(z.lo)<0)mpfr_set_zero(z.lo,0);
  if(mpfr_cmp_ui(z.hi,1)>0)mpfr_set_ui(z.hi,1,MPFR_RNDN);
  if(mpfr_cmp_ui(z.lo,1)>0)mpfr_set_ui(z.lo,1,MPFR_RNDN);
  return z;
}

static Interval bin_cdf(int n,int k,double p){
  if(k<0)return Interval(0.0);
  if(k>=n||p==0)return Interval(1.0);
  if(p==1)return Interval(0.0);
  Interval P(p),Q=sub(Interval(1.0),P),term,ratio,sum;
  bool left=(k+1<=n-k);
  if(left){mpfr_pow_ui(term.lo,Q.lo,n,MPFR_RNDD);mpfr_pow_ui(term.hi,Q.hi,n,MPFR_RNDU);ratio=divi(P,Q);}
  else{mpfr_pow_ui(term.lo,P.lo,n,MPFR_RNDD);mpfr_pow_ui(term.hi,P.hi,n,MPFR_RNDU);ratio=divi(Q,P);}
  sum=term;int stop=left?k:n-k-1;
  for(int j=0;j<stop;++j){
    mpfr_mul(term.lo,term.lo,ratio.lo,MPFR_RNDD);mpfr_mul(term.hi,term.hi,ratio.hi,MPFR_RNDU);
    mpfr_mul_ui(term.lo,term.lo,n-j,MPFR_RNDD);mpfr_mul_ui(term.hi,term.hi,n-j,MPFR_RNDU);
    mpfr_div_ui(term.lo,term.lo,j+1,MPFR_RNDD);mpfr_div_ui(term.hi,term.hi,j+1,MPFR_RNDU);
    mpfr_add(sum.lo,sum.lo,term.lo,MPFR_RNDD);mpfr_add(sum.hi,sum.hi,term.hi,MPFR_RNDU);
  }
  return unit(left?sum:sub(Interval(1.0),sum));
}
static Interval beta_cdf(int a,int b,double p){return unit(sub(Interval(1.0),bin_cdf(a+b-1,a-1,p)));}

static Interval beta_quantile(int a,int b,double pl,double pu,double proposal){
  Interval answer;mpfr_set_ui(answer.hi,1,MPFR_RNDN);
  if(pl==0&&pu==0)return Interval(0.0);
  if(pl==1&&pu==1)return Interval(1.0);
  if(!R_FINITE(proposal)||proposal<0||proposal>1)proposal=.5;
  double radius=std::max(8*std::max(proposal-std::nextafter(proposal,0.0),std::nextafter(proposal,1.0)-proposal),std::numeric_limits<double>::denorm_min());
  bool lok=pl==0,hik=pu==1;
  for(int j=0;j<70&&!(lok&&hik);++j){
    double l=std::max(0.0,std::nextafter(proposal-radius,0.0));
    double u=std::min(1.0,std::nextafter(proposal+radius,1.0));
    if(!lok){Interval f=beta_cdf(a,b,l);if(mpfr_cmp_d(f.hi,pl)<=0){mpfr_set_d(answer.lo,l,MPFR_RNDD);lok=true;}}
    if(!hik){Interval f=beta_cdf(a,b,u);if(mpfr_cmp_d(f.lo,pu)>=0){mpfr_set_d(answer.hi,u,MPFR_RNDU);hik=true;}}
    radius*=2;
  }
  return answer;
}
extern "C" SEXP lr_beta(SEXP aa,SEXP bb,SEXP pp,SEXP qq,SEXP proposed){
  bits=192;int n=LENGTH(aa);SEXP out=PROTECT(allocMatrix(REALSXP,n,2));
  for(int i=0;i<n;++i){
    int a=INTEGER(aa)[i],b=INTEGER(bb)[i];double pl=REAL(pp)[i],pu=REAL(qq)[i];
    if(a<1||b<1||pl<0||pu>1||pl>pu)Rf_error("Invalid beta enclosure input");
    Interval z=beta_quantile(a,b,pl,pu,REAL(proposed)[i]);REAL(out)[i]=down(z);REAL(out)[i+n]=up(z);
    if(i%128==0)R_CheckUserInterrupt();
  }
  UNPROTECT(1);return out;
}
extern "C" SEXP lr_binomial(SEXP nn,SEXP kk,SEXP pp,SEXP upper){
  bits=192;int n=LENGTH(kk),m=asInteger(nn);SEXP out=PROTECT(allocMatrix(REALSXP,n,2));
  for(int i=0;i<n;++i){Interval z=bin_cdf(m,INTEGER(kk)[i],REAL(pp)[i]);if(asLogical(upper))z=unit(sub(Interval(1.0),z));REAL(out)[i]=down(z);REAL(out)[i+n]=up(z);}
  UNPROTECT(1);return out;
}
extern "C" SEXP lr_psi(SEXP xx,SEXP aa){
  bits=192;int n=LENGTH(xx),J=LENGTH(aa);Interval T,A,h((J-1)/2.0);
  for(int j=0;j<J;++j){Interval a(REAL(aa)[j]);T=add(T,a);A=add(A,divi(Interval(1.0),mul(Interval(12.0),a)));}
  Interval rem=divi(Interval(1.0),mul(Interval(12.0),T));
  SEXP out=PROTECT(allocMatrix(REALSXP,n,2));
  for(int i=0;i<n;++i){Interval ans(1.0);double x0=REAL(xx)[i];
    if(x0>0&&R_FINITE(x0)){
      Interval x(x0),r=divi(add(h,sqrti(add(mul(h,h),mul(Interval(4.0),mul(A,x))))),mul(Interval(2.0),x));r=unit(r);
      Interval exponent=add(sub(mul(sub(r,Interval(1.0)),x),mul(h,logi(r))),add(divi(A,r),rem));
      ans=unit(expi(exponent));
    }else if(x0==R_PosInf)ans=Interval(0.0);
    REAL(out)[i]=down(ans);REAL(out)[i+n]=up(ans);
  }
  UNPROTECT(1);return out;
}
extern "C" SEXP lr_arithmetic(SEXP xx,SEXP yy,SEXP operation){
  bits=192;int n=std::max(LENGTH(xx),LENGTH(yy)),op=asInteger(operation);SEXP out=PROTECT(allocMatrix(REALSXP,n,2));
  for(int i=0;i<n;++i){Interval x(REAL(xx)[i%LENGTH(xx)]),y(REAL(yy)[i%LENGTH(yy)]),z;
    if(op==1)z=add(x,y);else if(op==2)z=sub(x,y);else if(op==3)z=mul(x,y);else if(op==4)z=divi(x,y);else if(op==5)z=logi(x);else z=expi(x);
    REAL(out)[i]=down(z);REAL(out)[i+n]=up(z);
  }
  UNPROTECT(1);return out;
}
extern "C" SEXP lr_by_thresholds(SEXP kk,SEXP alpha){
  bits=192;int K=asInteger(kk);Interval H;for(int i=1;i<=K;++i)H=add(H,divi(Interval(1.0),Interval((double)i)));
  SEXP out=PROTECT(allocMatrix(REALSXP,K,2));
  for(int i=1;i<=K;++i){Interval z=divi(mul(Interval(asReal(alpha)),Interval((double)i)),mul(Interval((double)K),H));REAL(out)[i-1]=down(z);REAL(out)[i-1+K]=up(z);}
  UNPROTECT(1);return out;
}
