#include <mpfr.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <R.h>
#include <Rinternals.h>

static mpfr_prec_t bits=128;
struct Interval {
  mpfr_t lo,hi;
  Interval(){mpfr_init2(lo,bits);mpfr_init2(hi,bits);mpfr_set_zero(lo,0);mpfr_set_zero(hi,0);}
  explicit Interval(double x):Interval(){mpfr_set_d(lo,x,MPFR_RNDD);mpfr_set_d(hi,x,MPFR_RNDU);}
  explicit Interval(mpfr_srcptr x):Interval(){mpfr_set(lo,x,MPFR_RNDD);mpfr_set(hi,x,MPFR_RNDU);}
  Interval(const Interval& x):Interval(){mpfr_set(lo,x.lo,MPFR_RNDD);mpfr_set(hi,x.hi,MPFR_RNDU);}
  Interval& operator=(const Interval& x){if(this!=&x){mpfr_set(lo,x.lo,MPFR_RNDD);mpfr_set(hi,x.hi,MPFR_RNDU);}return *this;}
  ~Interval(){mpfr_clear(lo);mpfr_clear(hi);}
};
static Interval add(const Interval&a,const Interval&b){Interval z;mpfr_add(z.lo,a.lo,b.lo,MPFR_RNDD);mpfr_add(z.hi,a.hi,b.hi,MPFR_RNDU);return z;}
static Interval sub(const Interval&a,const Interval&b){Interval z;mpfr_sub(z.lo,a.lo,b.hi,MPFR_RNDD);mpfr_sub(z.hi,a.hi,b.lo,MPFR_RNDU);return z;}
static Interval mul(const Interval&a,const Interval&b){
  Interval z;mpfr_t t;mpfr_init2(t,bits);mpfr_set_inf(z.lo,1);mpfr_set_inf(z.hi,-1);
  mpfr_srcptr aa[]={a.lo,a.hi},bb[]={b.lo,b.hi};
  for(auto x:aa)for(auto y:bb){mpfr_mul(t,x,y,MPFR_RNDD);if(mpfr_less_p(t,z.lo))mpfr_set(z.lo,t,MPFR_RNDD);mpfr_mul(t,x,y,MPFR_RNDU);if(mpfr_greater_p(t,z.hi))mpfr_set(z.hi,t,MPFR_RNDU);}mpfr_clear(t);return z;
}
static Interval divi(const Interval&a,const Interval&b){Interval inv;mpfr_ui_div(inv.lo,1,b.hi,MPFR_RNDD);mpfr_ui_div(inv.hi,1,b.lo,MPFR_RNDU);return mul(a,inv);}
static Interval logi(const Interval&a){Interval z;mpfr_log(z.lo,a.lo,MPFR_RNDD);mpfr_log(z.hi,a.hi,MPFR_RNDU);return z;}
static Interval positive(const Interval&a){Interval z(a);if(mpfr_sgn(z.lo)<0)mpfr_set_zero(z.lo,0);if(mpfr_sgn(z.hi)<0)mpfr_set_zero(z.hi,0);return z;}
static Interval point_upper(const Interval&a){return Interval(a.hi);}
static Interval point_lower(const Interval&a){return Interval(a.lo);}
static double down(const Interval&a){return mpfr_get_d(a.lo,MPFR_RNDD);}
static double up(const Interval&a){return mpfr_get_d(a.hi,MPFR_RNDU);}

extern "C" SEXP certified_box(SEXP lower,SEXP upper,SEXP shape,SEXP tolerance,SEXP maxbits){
  const int J=LENGTH(shape);const double tol=asReal(tolerance);double ans[]={0,R_PosInf,R_PosInf,0,0,1};
  for(int j=0;j<J;++j)if(!R_FINITE(REAL(lower)[j])||!R_FINITE(REAL(upper)[j])||!R_FINITE(REAL(shape)[j])||REAL(lower)[j]<0||REAL(upper)[j]<=0||REAL(lower)[j]>REAL(upper)[j]||REAL(shape)[j]<=0)goto done;
  for(bits=128;bits<=asInteger(maxbits);bits*=2){
    std::vector<Interval> a,l,u,q,y;
    Interval T,L,U;
    for(int j=0;j<J;++j){a.emplace_back(REAL(shape)[j]);l.emplace_back(REAL(lower)[j]);u.emplace_back(REAL(upper)[j]);T=add(T,a[j]);L=add(L,l[j]);U=add(U,u[j]);}
    bool contains=true;
    for(int j=0;j<J&&contains;++j)for(int k=0;k<J;++k){Interval x=mul(l[j],a[k]),z=mul(u[k],a[j]);if(mpfr_cmp(x.hi,z.lo)>0){contains=false;break;}}
    if(contains){ans[0]=ans[1]=ans[2]=0;ans[3]=bits;ans[5]=0;break;}

    if(mpfr_zero_p(L.hi)){ans[0]=ans[1]=ans[2]=0;ans[3]=bits;ans[5]=0;break;}
    for(int j=0;j<J;++j)q.push_back(divi(a[j],T));
    bool singleton=true;for(int j=0;j<J;++j)if(mpfr_cmp(l[j].lo,u[j].hi)!=0)singleton=false;
    if(singleton||J==2){
      std::vector<Interval> corner;
      if(singleton)corner=l;
      else {

        Interval leftA=mul(a[0],u[1]),leftB=mul(a[1],l[0]);
        Interval rightA=mul(a[0],l[1]),rightB=mul(a[1],u[0]);
        if(mpfr_cmp(leftA.hi,leftB.lo)<0){corner.push_back(l[0]);corner.push_back(u[1]);}
        else if(mpfr_cmp(rightA.lo,rightB.hi)>0){corner.push_back(u[0]);corner.push_back(l[1]);}
        else continue;
      }
      Interval S,value;for(auto&c:corner)S=add(S,c);
      for(int j=0;j<J;++j)value=add(value,mul(a[j],logi(divi(q[j],divi(corner[j],S)))));
      value=positive(value);ans[0]=down(value);ans[1]=up(value);ans[2]=up(sub(Interval(ans[1]),Interval(ans[0])));ans[3]=bits;ans[4]=0;
      if(R_FINITE(ans[2])&&ans[2]>=0&&ans[2]<=tol){ans[5]=0;break;}
      continue;
    }
    mpfr_t left,right,mid,total,t;mpfr_inits2(bits,left,right,mid,total,t,(mpfr_ptr)0);
    mpfr_set(left,L.lo,MPFR_RNDN);mpfr_set(right,U.hi,MPFR_RNDN);
    for(int it=1;it<=static_cast<int>(bits)+40;++it){
      mpfr_add(mid,left,right,MPFR_RNDN);mpfr_div_2ui(mid,mid,1,MPFR_RNDN);mpfr_set_zero(total,0);y.clear();
      for(int j=0;j<J;++j){mpfr_mul(t,mid,a[j].lo,MPFR_RNDN);mpfr_div(t,t,T.lo,MPFR_RNDN);if(mpfr_less_p(t,l[j].lo))mpfr_set(t,l[j].lo,MPFR_RNDN);if(mpfr_greater_p(t,u[j].hi))mpfr_set(t,u[j].hi,MPFR_RNDN);y.emplace_back(t);mpfr_add(total,total,t,MPFR_RNDN);}
      if(mpfr_cmp(total,mid)>=0)mpfr_set(left,mid,MPFR_RNDN);else mpfr_set(right,mid,MPFR_RNDN);
      if(it<40||it%8)continue;
      Interval S;for(auto& yy:y)S=add(S,yy);
      Interval primal,linear,residual;
      for(int j=0;j<J;++j){
        Interval z=divi(y[j],S);primal=add(primal,mul(a[j],logi(divi(q[j],z))));

        mpfr_mul(t,a[j].lo,S.lo,MPFR_RNDN);mpfr_div(t,t,T.lo,MPFR_RNDN);mpfr_div(t,t,y[j].lo,MPFR_RNDN);Interval b(t);
        linear=add(linear,mul(a[j],logi(b)));Interval diff=sub(b,Interval(1.0));
        residual=add(residual,mul(diff,mpfr_cmp_ui(b.lo,1)>=0?u[j]:l[j]));
      }
      Interval rho=add(Interval(1.0),divi(point_upper(positive(residual)),point_lower(L)));
      Interval dual=positive(sub(linear,mul(T,logi(rho))));primal=positive(primal);

      double dl=down(dual),du=up(primal);Interval gap=sub(Interval(du),Interval(dl));
      ans[0]=dl;ans[1]=du;ans[2]=up(gap);ans[3]=bits;ans[4]=it;
      if(R_FINITE(ans[2])&&ans[2]>=0&&ans[2]<=tol){ans[5]=0;break;}
    }
    mpfr_clears(left,right,mid,total,t,(mpfr_ptr)0);
    if(ans[5]==0)break;
  }
done:
  if(ans[5]){ans[0]=0;ans[2]=ans[1];}
  SEXP out=PROTECT(allocVector(REALSXP,6));for(int i=0;i<6;++i)REAL(out)[i]=ans[i];UNPROTECT(1);return out;
}

extern "C" SEXP certified_subset(SEXP counts,SEXP shape,SEXP nkeep,SEXP weights,SEXP trial){
  bits=192;const int J=LENGTH(shape),N=INTEGER(getAttrib(counts,R_DimSymbol))[0],keep=asInteger(nkeep);
  Interval T,den,base,linear;std::vector<Interval>a,b;std::vector<Interval> scores,totals;
  Interval rho(asReal(trial));
  for(int j=0;j<J;++j){a.emplace_back(REAL(shape)[j]);b.emplace_back(REAL(weights)[j]);T=add(T,a[j]);linear=add(linear,mul(a[j],logi(b[j])));base=add(base,mul(a[j],sub(b[j],rho)));}
  for(int i=0;i<N;++i){Interval score,tot;for(int j=0;j<J;++j){Interval v(REAL(counts)[i+N*j]);tot=add(tot,v);score=add(score,mul(v,sub(b[j],rho)));}scores.push_back(point_upper(score));totals.push_back(point_lower(tot));}
  std::sort(scores.begin(),scores.end(),[](const Interval&x,const Interval&y){return mpfr_greater_p(x.hi,y.hi);});
  std::sort(totals.begin(),totals.end(),[](const Interval&x,const Interval&y){return mpfr_less_p(x.lo,y.lo);});den=T;for(int i=0;i<keep;++i){base=add(base,scores[i]);den=add(den,totals[i]);}
  Interval rhoUpper=add(rho,divi(point_upper(positive(base)),point_lower(den)));
  Interval lower=positive(sub(linear,mul(T,logi(rhoUpper))));return ScalarReal(down(lower));
}
