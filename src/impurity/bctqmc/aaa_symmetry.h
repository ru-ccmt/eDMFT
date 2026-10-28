// Symmetry-constrained fits of independent real scalar bath channels.
// The generic matrix AAA implementation remains the pole initializer.
#ifndef EDMFT_AAA_SYMMETRY_H
#define EDMFT_AAA_SYMMETRY_H
#include "aaa.h"
#include <numeric>
#include <bit>
#include <cstdint>

namespace edmft::aaa {
struct ChannelSymmetry {
  std::vector<unsigned char> even;
  bool active() const {return std::find(even.begin(),even.end(),1)!=even.end();}
};

// Only purely anomalous CIX channels are eligible for automatic parity
// detection. A complex/odd-frequency channel is never forced even in auto.
inline ChannelSymmetry detect_even_channels(const std::vector<MatrixXcd>& f,
    const std::vector<unsigned char>& anomalous,const std::string& mode,
    double tolerance,std::ostream* log=nullptr) {
  if(mode!="auto" && mode!="none" && mode!="even")
    throw std::invalid_argument("aaa_anomalous_symmetry must be auto, none or even");
  if(!(tolerance>=0.) || !std::isfinite(tolerance))
    throw std::invalid_argument("aaa_symmetry_tolerance must be finite and nonnegative");
  ChannelSymmetry s;s.even.assign(anomalous.size(),0);
  if(mode=="none")return s;
  for(int c=0;c<int(anomalous.size());++c)if(anomalous[c]) {
    double scale=0.,imag=0.;for(const auto& x:f){
      if(x.rows()!=int(anomalous.size()) || x.cols()!=x.rows() || !x.allFinite())
        throw std::invalid_argument("invalid AAA symmetry channel data");
      scale=std::max(scale,std::abs(x(c,c)));imag=std::max(imag,std::abs(x(c,c).imag()));
    }
    // No amplitude cutoff: an arbitrarily small nonzero bath is retained.
    const double defect=scale>0. ? imag/scale : 0.;
    if(mode=="even" && defect>tolerance)
      throw std::invalid_argument("anomalous input violates requested even-frequency symmetry");
    s.even[c]=scale>0. && defect<=tolerance;
    if(log)*log<<"AAA anomalous symmetry: channel="<<c<<" even="<<int(s.even[c])
      <<" input_relative_defect="<<defect<<" tolerance="<<tolerance<<std::endl;
  }
  return s;
}

namespace parity_detail {
// Intel's fast floating-point mode may optimize away isfinite/allFinite.
// Check the representation at the nonlinear optimizer boundary instead.
inline bool finite_number(double x) {
  static_assert(sizeof(double)==sizeof(std::uint64_t) && std::numeric_limits<double>::is_iec559);
  return (std::bit_cast<std::uint64_t>(x)&UINT64_C(0x7ff0000000000000))!=UINT64_C(0x7ff0000000000000);
}
inline bool finite_number(Complex x) {return finite_number(x.real()) && finite_number(x.imag());}
template<class Derived> inline bool finite_values(const Eigen::MatrixBase<Derived>& x) {
  for(Eigen::Index i=0;i<x.size();++i)if(!finite_number(x.derived().coeff(i)))return false;
  return true;
}

struct Layout {
  int pairs=0;
  VectorXd x; // positive pair magnitudes, followed by unrestricted normal poles
  int count() const {return x.size()+pairs;}
  VectorXd poles() const {
    VectorXd p(count());for(int k=0;k<pairs;++k){p(2*k)=-x(k);p(2*k+1)=x(k);}
    for(int k=pairs;k<x.size();++k)p(pairs+k)=x(k);return p;
  }
};
inline void validate(const std::vector<Complex>& z,const std::vector<MatrixXcd>& f,
                     const ChannelSymmetry& s) {
  detail::validate_samples(z,f);
  if(f[0].rows()!=f[0].cols() || s.even.size()!=size_t(f[0].rows()))
    throw std::invalid_argument("AAA even symmetry requires independent scalar channels");
  for(const auto& a:f)for(int i=0;i<a.rows();++i)for(int j=0;j<a.cols();++j)
    if(!std::isfinite(a(i,j).real()) || !std::isfinite(a(i,j).imag()) ||
       (i!=j && a(i,j)!=Complex(0.,0.)))
      throw std::invalid_argument("AAA even symmetry requires finite diagonal channel data");
}
inline bool needs_pair(const std::vector<MatrixXcd>& f,const ChannelSymmetry& s) {
  for(int c=0;c<int(s.even.size());++c)if(s.even[c])
    for(const auto& a:f)if(a(c,c)!=Complex(0.,0.))return true;
  return false;
}
inline Layout paired_seed(const VectorXd& seed) {
  std::vector<double> e;bool zero=false;
  for(double p:seed){if(!std::isfinite(p))throw std::invalid_argument("nonfinite AAA pole");if(p==0.)zero=true;else e.push_back(std::abs(p));}
  std::sort(e.begin(),e.end());e.erase(std::unique(e.begin(),e.end()),e.end());
  Layout l;l.pairs=e.size();l.x.resize(e.size()+zero);
  for(int i=0;i<int(e.size());++i)l.x(i)=e[i];if(zero)l.x(e.size())=0.;return l;
}
inline Layout existing_layout(const PoleExpansion& source,const ChannelSymmetry& s) {
  const int n=source.poles.size();std::vector<bool> seen(n,false);std::vector<double> e,q;
  for(int i=0;i<n;++i)if(!seen[i]) {
    const double p=source.poles(i);int mate=-1;
    if(p!=0.)for(int j=i+1;j<n;++j)if(!seen[j] && source.poles(j)==-p){mate=j;break;}
    seen[i]=true;
    if(mate>=0){seen[mate]=true;e.push_back(std::abs(p));}
    else {
      for(int c=0;c<int(s.even.size());++c)if(s.even[c] && source.residues[i](c,c)!=Complex(0.,0.))
        throw std::invalid_argument("unpaired nonzero anomalous AAA residue");
      q.push_back(p);
    }
  }
  Layout l;l.pairs=e.size();l.x.resize(e.size()+q.size());
  for(int i=0;i<int(e.size());++i)l.x(i)=e[i];for(int i=0;i<int(q.size());++i)l.x(e.size()+i)=q[i];return l;
}
struct Fit {
  PoleExpansion expansion;
  std::vector<double> errors;
  double relative=0.,worst=0.,maximum=0.,norm=0.;
};
inline Fit measure(const std::vector<Complex>& z,const std::vector<MatrixXcd>& f,PoleExpansion expansion) {
  if(!finite_values(expansion.poles))throw std::runtime_error("nonfinite constrained AAA poles");
  for(const auto& residue:expansion.residues)
    if(!finite_values(residue))throw std::runtime_error("nonfinite constrained AAA residues");
  Fit fit;fit.expansion=std::move(expansion);const int nc=f[0].rows();
  double error=0.,scale=0.;std::vector<double> e(nc,0.),b(nc,0.);
  for(int i=0;i<int(z.size());++i){const auto v=fit.expansion.evaluate(z[i]);const MatrixXcd d=v-f[i];
    error+=d.squaredNorm();scale+=f[i].squaredNorm();fit.maximum=std::max(fit.maximum,d.norm());
    for(int c=0;c<nc;++c){e[c]+=std::norm(d(c,c));b[c]+=std::norm(f[i](c,c));}}
  fit.norm=std::sqrt(error);fit.relative=std::sqrt(error/std::max(scale,std::numeric_limits<double>::min()));fit.errors.resize(nc);
  for(int c=0;c<nc;++c){fit.errors[c]=b[c]>0. ? std::sqrt(e[c]/b[c]) : -1.;fit.worst=std::max(fit.worst,fit.errors[c]);}
  if(!finite_number(fit.norm) || !finite_number(fit.relative) ||
     !finite_number(fit.maximum) || !finite_number(fit.worst))
    throw std::runtime_error("nonfinite constrained AAA fit");
  return fit;
}
inline Fit fit_layout(const std::vector<Complex>& z,const std::vector<MatrixXcd>& f,
                      const ChannelSymmetry& s,const Layout& l) {
  if(l.pairs<0 || l.pairs>l.x.size())throw std::invalid_argument("invalid constrained AAA layout");
  if(!finite_values(l.x))throw std::runtime_error("nonfinite constrained AAA trial poles");
  // A collapsed pair has a zero anomalous design column. Reject this
  // singular trial before QR; optimize() retains its best finite candidate.
  for(int k=0;k<l.pairs;++k)if(!(l.x(k)>0.))
    throw std::runtime_error("collapsed constrained AAA pole pair");
  const int n=z.size(),nc=f[0].rows();const VectorXd p=l.poles();
  std::vector<int> general,even;
  for(int c=0;c<nc;++c)(s.even[c] ? even : general).push_back(c);
  Fit fit;fit.expansion.poles=p;fit.expansion.residues.assign(p.size(),MatrixXcd::Zero(nc,nc));
  const auto solve=[&](const std::vector<int>& channels,bool tied) {
    if(channels.empty())return;
    const int cols=tied ? l.pairs : p.size();
    MatrixXd a(2*n,cols),b(2*n,channels.size());
    for(int j=0;j<cols;++j)for(int i=0;i<n;++i){
      Complex v=tied ? 1./(z[i]-l.x(j))-1./(z[i]+l.x(j)) : 1./(z[i]-p(j));
      a(i,j)=v.real();a(n+i,j)=v.imag();
    }
    for(int j=0;j<int(channels.size());++j)for(int i=0;i<n;++i){b(i,j)=f[i](channels[j],channels[j]).real();b(n+i,j)=f[i](channels[j],channels[j]).imag();}
    if(cols==0){if(b.squaredNorm()!=0.)throw std::invalid_argument("nonzero even anomalous bath requires at least two poles");return;}
    if(!finite_values(a) || !finite_values(b))
      throw std::runtime_error("nonfinite constrained AAA least-squares input");
    for(int j=0;j<cols;++j)if(a.col(j).cwiseAbs().maxCoeff()==0.)
      throw std::runtime_error("zero constrained AAA design column");
    const MatrixXd c=a.colPivHouseholderQr().solve(b);
    if(!finite_values(c))throw std::runtime_error("nonfinite constrained AAA least-squares solution");
    for(int j=0;j<int(channels.size());++j)for(int k=0;k<cols;++k){
      if(tied){fit.expansion.residues[2*k](channels[j],channels[j])=-c(k,j);fit.expansion.residues[2*k+1](channels[j],channels[j])=c(k,j);}
      else fit.expansion.residues[k](channels[j],channels[j])=c(k,j);
    }
  };
  solve(general,false);solve(even,true);
  return measure(z,f,std::move(fit.expansion));
}
inline bool better(const Fit& a,const Fit& b){return a.worst<b.worst || (a.worst==b.worst && a.relative<b.relative);}
struct Objective {
  const std::vector<Complex>& z;const std::vector<MatrixXcd>& f;
  const ChannelSymmetry& s;int pairs;int evaluations=0;
  VectorXd best_x;
  double best_norm=std::numeric_limits<double>::infinity();
  double operator()(const VectorXd& x,VectorXd& grad) {
    Layout l{pairs,x};const Fit fit=fit_layout(z,f,s,l);grad=VectorXd::Zero(x.size());
    for(int i=0;i<int(z.size());++i){const MatrixXcd d=fit.expansion.evaluate(z[i])-f[i];
      for(int k=0;k<x.size();++k)for(int c=0;c<f[0].rows();++c){
        Complex derivative;
        if(k<pairs)derivative=fit.expansion.residues[2*k+1](c,c)/std::pow(z[i]-x(k),2)-fit.expansion.residues[2*k](c,c)/std::pow(z[i]+x(k),2);
        else derivative=fit.expansion.residues[pairs+k](c,c)/std::pow(z[i]-x(k),2);
        grad(k)+=(std::conj(d(c,c))*derivative).real();
      }}
    if(fit.norm>std::numeric_limits<double>::min())grad/=fit.norm;
    if(!finite_values(grad))throw std::runtime_error("nonfinite constrained AAA gradient");
    if(fit.norm<best_norm){best_norm=fit.norm;best_x=x;}
    ++evaluations;return fit.norm;
  }
};
// Optimize log(E) for symmetric pair magnitudes. This preserves the same
// positive-pole fitting domain while removing its singular E=0 boundary.
struct PositivePairObjective {
  Objective objective;
  VectorXd physical(const VectorXd& variables) const {
    VectorXd x=variables;
    for(int k=0;k<objective.pairs;++k)x(k)=std::exp(variables(k));
    return x;
  }
  double operator()(const VectorXd& variables,VectorXd& grad) {
    const VectorXd x=physical(variables);
    const double value=objective(x,grad);
    for(int k=0;k<objective.pairs;++k)grad(k)*=x(k);
    if(!finite_values(grad))throw std::runtime_error("nonfinite constrained AAA log-pair gradient");
    return value;
  }
};
inline int optimize(const std::vector<Complex>& z,const std::vector<MatrixXcd>& f,
    const ChannelSymmetry& s,Layout& l,double bound,int iterations,double tolerance,
    int* evaluations=nullptr) {
  if(iterations==0 || l.x.size()==0)return 0;
  if(!(bound>0.) || !finite_number(bound))throw std::invalid_argument("invalid AAA refinement bound");
  PositivePairObjective objective{{z,f,s,l.pairs}};LBFGSpp::LBFGSBParam<double> param;
  param.max_iterations=iterations;param.epsilon=tolerance;param.epsilon_rel=tolerance;
  LBFGSpp::LBFGSBSolver<double> solver(param);
  VectorXd lb=VectorXd::Constant(l.x.size(),-bound),ub=VectorXd::Constant(l.x.size(),bound);
  VectorXd variables=l.x;
  for(int k=0;k<l.pairs;++k){
    if(!(l.x(k)>0.) || !finite_number(l.x(k)))throw std::invalid_argument("invalid initial AAA pair");
    variables(k)=std::log(l.x(k));
    lb(k)=std::log(std::numeric_limits<double>::min());ub(k)=std::log(bound);
  }
  double fx=0.;
  const VectorXd initial=l.x;int done;
  try {done=solver.minimize(objective,variables,fx,lb,ub);
    l.x=objective.physical(variables);
    if(!finite_values(l.x))throw std::runtime_error("nonfinite AAA refinement parameters");
  }catch(const std::runtime_error&){
    // Singular/nonfinite trial or failed line search: preserve a finite,
    // symmetry-preserving candidate, including progress made before failure.
    l.x=objective.objective.best_x.size()==initial.size() ? objective.objective.best_x : initial;
    done=-1;
  }
  if(evaluations)*evaluations=objective.objective.evaluations;return done;
}
// Keep partners adjacent: summing all negative poles before all positive
// poles needlessly loses the exact imaginary cancellation on the Matsubara
// axis. Every layout emits (-E,+E) consecutively, then normal-only poles.
inline PoleExpansion paired_order(PoleExpansion a) {return a;}
}

inline double even_symmetry_defect(const PoleExpansion& p,const ChannelSymmetry& s) {
  double defect=0.;
  for(int i=0;i<p.poles.size();++i)for(int c=0;c<int(s.even.size());++c)if(s.even[c]){
    Complex sum=p.residues[i](c,c);
    if(p.poles(i)==0.)defect=std::max(defect,std::abs(sum));
    else {for(int j=0;j<p.poles.size();++j)if(p.poles(j)==-p.poles(i))sum+=p.residues[j](c,c);defect=std::max(defect,std::abs(sum));}}
  return defect;
}

// Analytic projection Phi(z)->[Phi(z)+Phi(-z)]/2. Normal channels retain
// exactly their seed residues and pole locations; added mirrors have zero
// normal residues. This avoids an ill-conditioned refit on a doubled basis.
inline PoleExpansion project_even_channels(const PoleExpansion& source,const ChannelSymmetry& s) {
  if(source.poles.size()==0)return source;
  const auto layout=parity_detail::paired_seed(source.poles);
  PoleExpansion out;out.poles=layout.poles();const int nc=s.even.size();
  out.residues.assign(out.poles.size(),MatrixXcd::Zero(nc,nc));
  for(int k=0;k<source.poles.size();++k){
    int row=-1;for(int j=0;j<out.poles.size();++j)if(out.poles(j)==source.poles(k)){row=j;break;}
    if(row<0)throw std::runtime_error("missing projected AAA pole");
    for(int c=0;c<nc;++c){
      if(!s.even[c])out.residues[row](c,c)+=source.residues[k](c,c);
      else if(source.poles(k)!=0.){
        const Complex half=.5*source.residues[k](c,c);
        out.residues[row](c,c)+=half;out.residues[row^1](c,c)-=half;
      }
    }
  }
  if(even_symmetry_defect(out,s)!=0.)throw std::runtime_error("analytic AAA parity projection failed");
  return out;
}

inline RefinementReport fit_even_channels(const std::vector<Complex>& z,
    const std::vector<MatrixXcd>& f,const VectorXd& seed,const ChannelSymmetry& s,
    double bound,int iterations=200,double tolerance=1e-10,std::ostream* log=nullptr) {
  using namespace parity_detail;validate(z,f,s);
  if(seed.size()==0)return {};
  if(needs_pair(f,s) && paired_seed(seed).pairs==0)
    throw std::invalid_argument("even anomalous bath needs a nonzero pole pair");
  bound=std::max(bound,1.05*seed.cwiseAbs().maxCoeff());
  const auto start=std::chrono::steady_clock::now();
  // Generic AAA/refinement is an initializer, never the consumed bath.
  // Project its rational function exactly, rather than performing another
  // unconstrained least-squares solve on a larger, ill-conditioned pole set.
  RefinementReport report;
  if(iterations>0)report=refine_poles(z,f,seed,-bound,bound,iterations,tolerance);
  else report.expansion=fit_residues(z,f,seed);
  Fit result=measure(z,f,project_even_channels(report.expansion,s));
  report.expansion=std::move(result.expansion);report.relative_rms_error=result.relative;report.max_error=result.maximum;
  report.objective=result.norm;report.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  if(even_symmetry_defect(report.expansion,s)!=0.)throw std::runtime_error("AAA even-pole constraint failed");
  if(log)*log<<"AAA constrained refinement: seed_poles="<<seed.size()<<" poles="<<report.expansion.poles.size()
    <<" stable_seed_projection=yes iterations="<<report.iterations<<" relative_rms="<<report.relative_rms_error
    <<" symmetry_defect=0 seconds="<<report.seconds<<std::endl;
  return report;
}

inline PoleExpansion refit_even_channels(const std::vector<Complex>& z,
    const std::vector<MatrixXcd>& f,const VectorXd& poles,const ChannelSymmetry& s) {
  parity_detail::validate(z,f,s);
  return project_even_channels(fit_residues(z,f,poles),s);
}

inline CompressionReport compress_even_channels(const std::vector<Complex>& z,
    const std::vector<MatrixXcd>& f,const PoleExpansion& source,int target,
    const ChannelSymmetry& s,std::ostream* log=nullptr) {
  using namespace parity_detail;validate(z,f,s);
  if(target<1 || target>source.poles.size())throw std::invalid_argument("invalid constrained AAA pole budget");
  if(even_symmetry_defect(source,s)!=0.)throw std::invalid_argument("input AAA expansion violates even symmetry");
  const bool pair_required=needs_pair(f,s);
  if(pair_required && target<2)throw std::invalid_argument("nonzero even anomalous bath requires at least two AAA poles");
  Layout l=existing_layout(source,s);Fit best=measure(z,f,source);
  if(target==source.poles.size()){
    best.expansion=source;
  }else {
    while(l.count()>target){bool have=false;Layout next;Fit chosen;
      // Each candidate removes a free normal pole, a complete pair, or one
      // member of a pair while setting that pole's anomalous residue to zero.
      const auto score=[&](Layout candidate){
        if(candidate.count()<target || (pair_required && candidate.pairs==0))return;
        Fit trial=fit_layout(z,f,s,candidate);
        if(!have || better(trial,chosen)){have=true;chosen=std::move(trial);next=std::move(candidate);}
      };
      for(int k=0;k<l.x.size();++k){
        Layout a;a.pairs=l.pairs-(k<l.pairs);a.x.resize(l.x.size()-1);
        for(int i=0,j=0;i<l.x.size();++i)if(i!=k)a.x(j++)=l.x(i);score(a);
        if(k<l.pairs)for(int sign:{-1,1}){Layout b=a;b.x.conservativeResize(a.x.size()+1);b.x(a.x.size())=sign*l.x(k);score(b);}
      }
      if(!have)throw std::runtime_error("no symmetry-preserving AAA compression candidate");
      l=std::move(next);best=std::move(chosen);
    }
    // Fixed-pole deletion is only an initializer. Refine the small retained
    // layout with equally normalized physical channels, retaining the better
    // worst-channel result. All original samples remain in the error report.
    std::vector<MatrixXcd> balanced=f;
    for(int c=0;c<f[0].rows();++c){double scale=0.;for(const auto& a:f)scale=std::hypot(scale,std::abs(a(c,c)));
      if(scale>0.)for(auto& a:balanced)a(c,c)/=scale;}
    const Layout initial=l;const double bound=std::max(1.,2.*source.poles.cwiseAbs().maxCoeff());
    // Small budgets are genuinely nonlinear fits. Nearby pair magnitudes
    // give inexpensive alternative starts; a failed search keeps its finite
    // initial candidate rather than corrupting the accepted expansion.
    const std::vector<double> starts=l.count()<=6 ? std::vector<double>{1.,.5,2.} : std::vector<double>{1.};
    int failed=0;
    for(double scale:starts){Layout trial=initial;
      for(int k=0;k<trial.pairs;++k)trial.x(k)=std::min(bound,trial.x(k)*scale);
      if(optimize(z,balanced,s,trial,bound,150,1e-10)<0)++failed;
      Fit refined=fit_layout(z,f,s,trial);if(better(refined,best)){l=std::move(trial);best=std::move(refined);}
    }
    if(log && failed)*log<<"AAA constrained compression: rejected_refinements="<<failed
      <<" finite_constrained_fallback=available"<<std::endl;
  }
  CompressionReport report;report.expansion=paired_order(std::move(best.expansion));report.relative_rms_error=best.relative;
  report.worst_channel_relative_rms=best.worst;report.channel_relative_rms=best.errors;
  if(even_symmetry_defect(report.expansion,s)!=0.)throw std::runtime_error("compressed AAA even-pole constraint failed");
  if(log)*log<<"AAA constrained compression: poles="<<report.expansion.poles.size()<<"/"<<source.poles.size()
    <<" pairs="<<l.pairs<<" worst_channel_relative_rms="<<report.worst_channel_relative_rms<<" symmetry_defect=0"<<std::endl;
  return report;
}
}
#endif
