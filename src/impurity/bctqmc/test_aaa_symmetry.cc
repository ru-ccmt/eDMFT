#include <fstream>
#include <sstream>
#include <iomanip>
#include "aaa.h"
using namespace edmft::aaa;
void require(bool p,const char* why){if(!p)throw std::runtime_error(why);}
void check(const PoleExpansion& p,const ChannelSymmetry& s){
 require(even_symmetry_defect(p,s)==0.,"nonzero pole/residue parity defect");
 for(double w:{.001,.1,.7,3.,31.}){
  auto a=p.evaluate(Complex(0,w)),b=p.evaluate(Complex(0,-w));
  for(int c=0;c<int(s.even.size());++c)if(s.even[c]){
   if(std::abs(a(c,c)-b(c,c))>=2e-13*std::max(1.,std::abs(a(c,c)))){std::cerr<<"frequency="<<w<<" plus="<<a(c,c)<<" minus="<<b(c,c)<<std::endl;}
   require(std::abs(a(c,c)-b(c,c))<2e-13*std::max(1.,std::abs(a(c,c))),"frequency parity failed");
   require(std::abs(a(c,c).imag())<2e-13*std::max(1.,std::abs(a(c,c))),"imaginary even Matsubara value");
  }
 }
 // Direct stable spectral-kernel check at nodes unrelated to the fit grid.
 for(double t:{0.,.17,.5,.83,1.})for(int c=0;c<int(s.even.size());++c)if(s.even[c]){
  double plus=0.,minus=0.,scale=0.;for(int k=0;k<p.poles.size();++k){
   auto kernel=[](double t,double w){return w>=0 ? std::exp(-t*w)/(1.+std::exp(-w)) : std::exp((1.-t)*w)/(1.+std::exp(w));};
   double a=p.residues[k](c,c).real()*kernel(t,35*p.poles(k));
   double b=p.residues[k](c,c).real()*kernel(t,-35*p.poles(k));plus+=a;minus+=b;scale+=std::abs(a)+std::abs(b);
  }require(std::abs(plus+minus)<2e-13*std::max(1.,scale),"imaginary-time reflection failed");
 }
}
int main(int argc,char**argv){try{
 std::cout<<std::setprecision(15);
 std::vector<Complex> z;std::vector<MatrixXcd> f;
 for(int i=0;i<90;++i){z.emplace_back(0,(2*i+1)*.06);MatrixXcd a=MatrixXcd::Zero(2,2);
  a(0,0)=.7/(z.back()+.6)+.2/(z.back()-.6)+.4/(z.back()-2.);
  a(1,1)=.08/(z.back()-.6)-.08/(z.back()+.6);f.push_back(a);}
 auto s=detect_even_channels(f,{0,1},"auto",1e-9,&std::cout);require(s.active()&&!s.even[0]&&s.even[1],"channel metadata failed");
 auto none=detect_even_channels(f,{0,1},"none",1e-9);require(!none.active(),"none failed");
 auto odd=f;for(int i=0;i<int(z.size());++i)odd[i](1,1)=.01*z[i]/(z[i]*z[i]-1.);
 require(!detect_even_channels(odd,{0,1},"auto",1e-9).active(),"odd frequency forced even");
 bool rejected=false;try{detect_even_channels(odd,{0,1},"even",1e-9);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"explicit invalid symmetry accepted");
 auto tiny=f;for(auto& a:tiny)a(1,1)*=1e-20;require(detect_even_channels(tiny,{0,1},"auto",1e-9).active(),"tiny bath dropped");
 auto zero=f;for(auto& a:zero)a(1,1)=0.;require(!detect_even_channels(zero,{0,1},"auto",1e-9).active(),"zero channel activates constraints");
 VectorXd seed(3);seed<<-.4,.9,2.1;
 VectorXd known(3);known<<-.6,.6,2.;
 auto tiny_fit=refit_even_channels(z,tiny,known,s);check(tiny_fit,s);
 double tiny_residue=0.;for(const auto& a:tiny_fit.residues)tiny_residue+=std::abs(a(1,1));
 require(tiny_residue>0.,"tiny nonzero anomalous residue discarded");
 parity_detail::Layout l=parity_detail::paired_seed(seed);parity_detail::Objective objective{z,f,s,l.pairs};VectorXd g;objective(l.x,g);
 for(int k=0;k<g.size();++k){auto a=l.x,b=l.x;a(k)+=1e-6;b(k)-=1e-6;VectorXd q;
  double finite=(objective(a,q)-objective(b,q))/2e-6;
  require(std::abs(finite-g(k))<2e-5*std::max({1.,std::abs(finite),std::abs(g(k))}),"constrained variable-projection gradient failed");}
 parity_detail::Layout mixed;mixed.pairs=1;mixed.x.resize(2);mixed.x<<.4,2.1;
 parity_detail::Objective mixed_objective{z,f,s,mixed.pairs};mixed_objective(mixed.x,g);
 for(int k=0;k<g.size();++k){auto a=mixed.x,b=mixed.x;a(k)+=1e-6;b(k)-=1e-6;VectorXd q;
  double finite=(mixed_objective(a,q)-mixed_objective(b,q))/2e-6;
  require(std::abs(finite-g(k))<2e-5*std::max({1.,std::abs(finite),std::abs(g(k))}),"normal-only pole gradient failed");}
 // The cluster failure first collapses its only pair to E=0, then
 // attempts NaN poles. Both must be rejected before Eigen QR, even with
 // fast floating-point compilation. Normal-only zero poles remain valid.
 for(double bad:{0.,-0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}){
  auto trial=mixed;trial.x(0)=bad;bool rejected_trial=false;
  try{parity_detail::fit_layout(z,f,s,trial);}catch(const std::runtime_error&){rejected_trial=true;}
  require(rejected_trial,"invalid pair reached QR");
 }
 {auto trial=mixed;trial.x(1)=std::numeric_limits<double>::quiet_NaN();bool rejected_trial=false;
  try{parity_detail::fit_layout(z,f,s,trial);}catch(const std::runtime_error&){rejected_trial=true;}
  require(rejected_trial,"nonfinite normal-only pole reached QR");}
 {auto trial=mixed;trial.x(1)=0.;auto fit=parity_detail::fit_layout(z,f,s,trial);
  require(parity_detail::finite_number(fit.norm),"normal-only zero pole was rejected");check(fit.expansion,s);}
 {parity_detail::PositivePairObjective trial{{z,f,s,mixed.pairs}};
  auto variables=mixed.x;variables(0)=std::log(variables(0));VectorXd grad;trial(variables,grad);
  for(int k=0;k<grad.size();++k){auto plus=variables,minus=variables;plus(k)+=1e-6;minus(k)-=1e-6;VectorXd tmp;
   double finite=(trial(plus,tmp)-trial(minus,tmp))/2e-6;
   require(std::abs(finite-grad(k))<2e-5*std::max({1.,std::abs(finite),std::abs(grad(k))}),"log-pair gradient failed");}
  variables(0)=-100.;require(trial.physical(variables)(0)>0.,"log pair reached singular zero");}
 {auto trial=mixed;const auto initial=trial.x;auto large=f;
  for(auto& value:large)value(1,1)*=1e4;
  // Force a numerically singular tiny-pair solve whose required residues
  // overflow, exercising the actual optimizer exception/fallback path.
  require(parity_detail::optimize(z,large,s,trial,2*std::numeric_limits<double>::min(),10,1e-10)<0,"overflowed residue solve accepted");
  require((trial.x-initial).norm()==0.,"optimizer failed to restore finite fallback");
  check(parity_detail::fit_layout(z,large,s,trial).expansion,s);}
 {parity_detail::Objective trial{z,f,s,mixed.pairs};VectorXd grad;
  trial(mixed.x,grad);auto bad=mixed.x;bad(0)=0.;
  try{trial(bad,grad);}catch(const std::runtime_error&){}
  require((trial.best_x-mixed.x).norm()==0.,"invalid trial replaced best finite fit");}
 std::cout<<"PASS collapsed/negative/nonfinite trials, normal-only zero pole, log-pair gradient and finite fallback"<<std::endl;
 auto full=fit_even_channels(z,f,seed,s,10.,120,1e-10,&std::cout);check(full.expansion,s);
 require(full.relative_rms_error<1e-6,"asymmetric normal/even anomalous exact fit failed");
 for(int count:{2,3,4,5}){auto compressed=compress_even_channels(z,f,full.expansion,count,s,&std::cout);check(compressed.expansion,s);require(compressed.expansion.poles.size()==count,"pole budget changed");}
 rejected=false;try{compress_even_channels(z,f,full.expansion,1,s);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"one pole silently destroys nonzero anomalous bath");
 auto exact=compress_even_channels(z,f,full.expansion,full.expansion.poles.size(),s);require((exact.expansion.poles-full.expansion.poles).norm()==0.,"uncapped pole table changed");
 auto a=fit_residues(z,f,seed),b=fit_residues(z,f,seed);require((a.evaluate(z[0])-b.evaluate(z[0])).norm()==0.,"legacy fit unstable");
 std::cout<<"PASS detection, tiny/zero/odd channels, gradient, asymmetric normal data, frequency/time symmetry, compression budgets and one-pole rejection"<<std::endl;
 if(argc>1){std::ifstream in(argv[1]);std::string line;std::vector<Complex> zs;std::vector<MatrixXcd> fs;
  while(std::getline(in,line)){if(line.empty()||line[0]=='#')continue;std::istringstream row(line);double w;row>>w;MatrixXcd a=MatrixXcd::Zero(4,4);
   for(int c=0;c<4;++c){double re,im;row>>re>>im;if(!row)throw std::runtime_error("bad cluster data");a(c,c)={re,im};}zs.emplace_back(0,w);fs.push_back(a);}
  auto sc=detect_even_channels(fs,{0,0,0,1},"auto",1e-9,&std::cout);
  Options opts;opts.tolerance=1e-10;auto raw=fit(zs,fs,opts);auto p=raw.real_poles(1e-7);
  auto accurate=fit_even_channels(zs,fs,p,sc,30.,200,1e-10,&std::cout);check(accurate.expansion,sc);
  {std::ofstream output("cluster_poles_accurate.dat");output<<std::setprecision(17);
   for(int k=0;k<accurate.expansion.poles.size();++k){output<<accurate.expansion.poles(k);for(int c=0;c<4;++c)output<<' '<<accurate.expansion.residues[k](c,c).real();output<<'\n';}}
  const double reference_bound=std::max(30.,1.05*p.cwiseAbs().maxCoeff());
  auto reference=refine_poles(zs,fs,p,-reference_bound,reference_bound,200,1e-10);
  for(double w:{.001,.1,.7,3.,31.}){
   const Complex z(0,w);const MatrixXcd old=reference.expansion.evaluate(z),now=accurate.expansion.evaluate(z);
   for(int c=0;c<3;++c)require(std::abs(now(c,c)-old(c,c))<2e-12*std::max(1.,std::abs(old(c,c))),"normal accurate bath changed during parity projection");
   const Complex projected=.5*(old(3,3)+reference.expansion.evaluate(-z)(3,3));
   require(std::abs(now(3,3)-projected)<2e-12*std::max(1.,std::abs(projected)),"analytic parity projection mismatch");
  }
  {std::ofstream output("cluster_poles_reference.dat");output<<std::setprecision(17);
   for(int k=0;k<reference.expansion.poles.size();++k){output<<reference.expansion.poles(k);for(int c=0;c<4;++c)output<<' '<<reference.expansion.residues[k](c,c).real();output<<'\n';}}
  require(accurate.relative_rms_error<5e-4,"accurate cluster fit lost accuracy");
  for(int count:{6,4,3,2}){auto compressed=compress_even_channels(zs,fs,accurate.expansion,count,sc,&std::cout);check(compressed.expansion,sc);require(compressed.expansion.poles.size()==count,"cluster cap changed");
   std::cout<<"CLUSTER poles="<<count<<" global="<<compressed.relative_rms_error<<" worst="<<compressed.worst_channel_relative_rms<<" frequencies=["<<compressed.expansion.poles.transpose()<<"]"<<std::endl;
   std::ofstream output(std::string("cluster_poles_")+std::to_string(count)+".dat");output<<std::setprecision(17);for(int k=0;k<count;++k){output<<compressed.expansion.poles(k);for(int c=0;c<4;++c)output<<' '<<compressed.expansion.residues[k](c,c).real();output<<'\n';}
  }
  std::cout<<"PASS saved superconducting cluster accurate fit and 6/4/3/2-pole compressed fits"<<std::endl;
 }
 return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
