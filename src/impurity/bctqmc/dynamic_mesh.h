#ifndef _DYNAMIC_MESH
#define _DYNAMIC_MESH

#include <iostream>
#include <fstream>
#include <deque>
#include <set>

#include "util.h"

#include "interpolate.h"
#include <string>

using namespace std;
using namespace blitz;

auto ReadSigma(const string& filename, int Npseudo){
  struct retVals{
    Array<double,1> om;//(nw);
    Array<double,2> App;//(Npseudo,nw);
  };
  ifstream inp(filename);
  retVals r;
  if (!inp.good() ) return r;
  int nw=0;
  string line;
  getline(inp,line);
  auto begin_pos = inp.tellg();
  while(inp){
    double w;
    inp>>w;
    if (!inp) break;
    nw += 1;
    getline(inp,line);
  }
  //cout<<nw<<endl;
  inp.clear();
  inp.seekg(begin_pos);
  r.om.resize(nw);
  r.App.resize(Npseudo,nw);
  for (int i=0; i<nw; ++i){
    inp>>r.om(i);
    for (int j=0; j<Npseudo; ++j) inp>>r.App(j,i);
    getline(inp,line);
  }
  return r;
}

Array<double,1> GiveDh(const Array<double,1>& om){
  int N = om.size();
  Array<double,1> dh(N); 
  dh(0) = 0.5*(om(1)-om(0));
  dh(N-1)= 0.5*(om(N-1)-om(N-2));
  for (int i=1; i<N-1; ++i)
    dh(i) = 0.5*(om(i+1)-om(i-1));
  return dh;
}

auto NewGrid(const Array<double,1>& fw2, const Array<double,1>& w, int N){
  Array<double,1> xn(N); xn=0;
  int nbins = fw2.size();
  int nw=w.size();
  double avgperbin = sum(fw2(Range(1,nw-2)))/(N-1.)*(nw-1.)/(nw-2.);
  //**** redefine the size of each bin  ****
  xn(0) = w(0);
  double prev=w(0), cur=w(0), dx=0;
  double thisbin = 0;
  int ibin = 0;
  for (int newbin=1; newbin<N-1; ++newbin){ // all but the last bin, which is 1.0
    while (thisbin < avgperbin){
      ibin+=1;
      thisbin += fw2(ibin); // adding widetilde{f}(g)
      prev = cur;           // g^old_{l-1}
      cur = w(ibin);     // g^old_{l}
    }

    //cout<<"newbin="<<newbin<<" ibin="<<ibin<<" thisbin="<<thisbin<<" avgperbin="<<avgperbin<<" prev="<<prev<<" cur="<<cur;
    // Explanation is in order : 
    //   g^new_l should be somewhere between g^old_l and g^old_{l-1}, because if we add the last point
    //      we exceeded I/Ng value, while withouth the last point, we do not have enough. Hence need to interpolate.
    //   prev    -- g^{old}_{l-1}
    //   cur     -- g^{old}_l
    //   thisbin -- Sm = f_{l-k}+.... +f_{l-2}+f_{l-1}+f_l
    //   we know that  Sm is just a bit more than we need, i.e., I/N_g, hence we need to compute how much more
    //   using linear interpolation :
    //   g^{new} = g_l - (g_l-g_{l-1}) * (f_{l-k}+....+f_{l-2}+f_{l-1}+f_l - I/N_g)/f_l
    //    clearly
    //         if I/N_g == f_{l-k}+....+f_{l-2}+f_{l-1}+f_l
    //            we will get g^{new} = g_l
    //     and if I/N_g == f_{l-k}+....+f_{l-2}+f_{l-1}
    //            we will get g^{new} = g_{l-1}
    //     and if I/N_g  is between the two possibilities, we will get linear interpolation between
    //     g_{l-1} and g_l
    //
    //double thisbin0 = thisbin;
    thisbin -= avgperbin;   // thisbin <- (f_{l-k}+....+f_{l-2}+f_{l-1}+f_l - I/N_g)
    dx=thisbin/fw2(ibin);    // dx is between 0 and 1
    // cur=g_l and prev=g_{l-1} and dx is the fraction of each we need to take for linear interpolation
    xn(newbin) = cur + (prev-cur)*dx; // linear interpolation between g_l and g_{l-1}
    //if (xn(newbin)-xn(newbin-1)<1e-10){
    //  cout<<"Here is a problem : xn["<<newbin<<"]="<<xn(newbin)<<" cur="<<cur<<" prev="<<prev<<" dx="<<dx<<" thisbin="<<thisbin<<" fw2["<<ibin<<"]="<<fw2(ibin)<<" fw2["<<ibin-1<<"]="<<fw2(ibin-1)<<" w["<<ibin<<"]="<<w(ibin)<<" w["<<ibin-1<<"]="<<w(ibin-1)<<endl;
    //}
    //cout<<" dx="<<dx<<" xn["<<newbin<<"]="<<xn(newbin)<<endl;
    //cout<<newbin<<" ibin="<<ibin<<" thisbin="<<thisbin0<<" prev="<<prev<<" cur="<<cur<<" newg="<<xn(newbin+1)<<endl;
  }
  if (xn(N-2)!=w(nw-1)){
    xn(N-1) = w(nw-1);
  }else{
    xn.resizeAndPreserve(N-1);
  }
  /*
  deque<double> remove;
  for (int i=xn.size()-1; i>1; --i) if (xn(i)-xn(i-1)<1e-12) remove.push_back(i);
  if (remove.size()>0){
    for (int i=0; i<remove.size(); ++i){
      int ir = remove[i];
      for (int j=ir+1; j<xn.size(); ++j) xn(j-1)=xn(j);
    }
    xn.resizeAndPreserve(xn.size()-remove.size());
  }
  */
  return xn;
}

auto FindExtra(const Array<double,1>& fx, const Array<double,1>& w, const Array<double,1>& xn, double small){
  struct retVal{
    deque<pair<int,double>> extra;
    Array<double,1> fn;
  };
  Array<double,1> dhn = GiveDh(xn);
  retVal r;
  r.fn.resize(xn.size());
  int pos=0;
  //for (int i=0; i<xn.size(); ++i) r.fn(i) = fx(InterpLeft(xn(i),pos,w));
  //for (int i=0; i<xn.size(); ++i) r.fn(i) = fx.LinInterp(InterpLeft(xn(i),pos,w));
  for (int i=0; i<xn.size(); ++i) r.fn(i) = linInterp(fx,InterpLeft(xn(i),pos,w));
  double integral = sum(r.fn*dhn);
  for (int i=0; i<xn.size()-1; ++i){
    double xt = 0.5*(xn(i)+xn(i+1));
    if (xn(i+1)-xn(i) < 1e-6*xn(i)*xn(i)) continue; // We don't want to points at the same place.
    //double ft = fx(Interp(xt,w));
    //double ft = fx.LinInterp(Interp(xt,w));
    double ft = linInterp(fx,Interp(xt,w));
    // xn[i+1],xt,xn[i]
    double dint = ft*0.5*(xn(i+1)-xn(i)) - r.fn(i+1)*0.5*(xt-xn(i)) - r.fn(i)*0.5*(xn(i+1)-xt);
    double tintegral = integral + dint;
    double change = abs(dint/integral);
    if (change>small){
      r.extra.push_back( make_pair(i+1,xt) );
      //cout<<"xn["<<i<<"]="<<xn(i)<<" xt="<<xt<<" xn["<<i+1<<"]="<<xn(i+1)<<endl;
    }
  }
  return r;
}

void CombineExtra(const deque<pair<int,double>>& extra, deque<pair<int,double>>& extras){
  set<int> iextras;
  for (int i=0; i<extras.size(); ++i) iextras.insert(extras[i].first);
  //iextras = [extras[i][0] for i in range(len(extras))]
  for (int i=0; i<extra.size(); ++i){
    int ii = extra[i].first;
    //if (not iextras.contains(ii)){
    if (iextras.find(ii)==iextras.end()){
      auto k=extras.begin();
      while(k!=extras.end() and k->first<ii) ++k;
      //if (k>0 and (extras[k-1][0]>ii or extras[k][0]<ii)):
      // print('WARNING extras['+str(k-1)+']=',extras[k-1][0], 'ii=', ii, 'extras['+str(k)+']=',extras[k][0])
      extras.insert(k, extra[i]);
      //cout<<"inserting extra["<<i<<"]=("<<extra[i].first<<","<<extra[i].second<<") at ("<< k->first<<","<<k->second<<"); extras=";
      //for (int j=0; j<extras.size(); ++j)
      //cout<<"("<<extras[j].first<<","<<extras[j].second<<") ";
      //cout<<endl;
    }
  }
}
Array<double,1> InserExtra(const Array<double,1>& xn, const deque<pair<int,double>>& extras){
  Array<double,1> xn2(xn.size()+extras.size());
  int iprevious=0, ishift=0;
  for (int i=0; i<extras.size(); ++i){
    int inext = extras[i].first;
    xn2(Range(iprevious+ishift,inext+ishift-1)) = xn(Range(iprevious,inext-1));
    xn2(inext+ishift) = extras[i].second;
    //#fn2[(iprevious+ishift):(inext+ishift)] = fn[iprevious:inext]
    //#fn2[inext+ishift] = extra[i][2]
    //cout<<"xn2["<<iprevious+ishift<<":"<<inext+ishift-1<<"]=xn["<<iprevious<<":"<<inext-1<<"]"<<endl;
    //cout<<"xn2["<<inext+ishift<<"]="<<extras[i].second<<endl;
    ishift += 1;
    iprevious=inext;
  }
  xn2(Range(iprevious+ishift,toEnd)) = xn(Range(iprevious,toEnd));
  //#fn2[(iprevious+ishift):] = fn[iprevious:]
  //cout<<"xn2["<<iprevious+ishift<<":]=xn["<<iprevious<<":]"<<endl;
  return xn2;
}

void AddMidpointsToImproveIntegral(Array<double,1>& xn, const Array<double,2>& fx, const Array<double,1>& om, ostream& log, double small=1e-5, int Nwmax=1000, int Nitt=20)
{
  // Here we try to add a point in the midpoint between two successive points on mesh xn. If the integral is changed for more than
  // amount small(1e-5), we accept the point. We then repeat the process up to Nitt(20) times,
  // but we stop when no midpoint point changes the integral for more than small, or, we exceed number of points Nwmax.
  int Npseudo = fx.extent(0);
  for (int itt=0; itt<Nitt; ++itt){
    deque<pair<int,double>> extras;
    for (int idim=0; idim<Npseudo; ++idim){
      auto [extra, fn] = FindExtra(fx(idim,Range::all()),om, xn, small);
      CombineExtra(extra, extras);
    }
    log<<"# "<<itt<<") Adding "<<extras.size()<<" points to xn="<<xn.size()<<endl;
    if (extras.size()==0) break;
    Array<double,1> xn2 = InserExtra(xn, extras);
    xn.resize(xn2.size());
    xn = xn2.copy();
    if (xn.size()>Nwmax) break;
  }
}

vector<Spline1D<double>> SplineSimpleFunctions(const Array<double,2>& App, const Array<double,1>& om)
{
  TinyVector<int,2> shp = App.shape();
  vector<Spline1D<double>> fx(shp[0]);
  if (shp[1]!=om.size()){
    cout<<"ERROR om.size="<<om.size()<<" while App.size[1]="<<shp[1]<<" which does not match"<<endl;
    return fx;
  }
  for (int i=0; i<shp[0]; ++i){
    fx[i].resize(om.size());
    for (int j=0; j<om.size(); ++j) fx[i][j] = App(i,j);
    fx[i].splineIt(om);
  }
  return fx;
}


/**
 * @brief Construct a global frequency-weight function fw2(ω) from a collection of positive “importance” curves App(p,ω).
 *
 * This routine implements the first stage of the adaptive-mesh construction used in our real-frequency solvers.
 * Given a set of curves App(p,ω) ≥ 0 (e.g. spectral weights per block/bath), it builds a single global weight
 * function that indicates where the frequency mesh should be dense.
 *
 * Algorithm:
 *   1) For each curve p, form p1(ω) = App(p,ω) * dh(ω) and normalize it by its integral Z_p = Σ p1.
 *      This removes trivial overall scale differences between different curves.
 *   2) Sum the normalized curves into fw0(ω) = Σ_p p1(ω) / Z_p.
 *   3) Smooth fw0 → fw1 by a simple local 3-point averaging filter.
 *   4) Sharpen the contrast with fw2(ω) = fw1(ω)^alpha.
 *
 * The result fw2(ω) is intended to be fed into NewGrid() (equal-weight grid generation) and optionally
 * into a refinement step (AddMidpointsToImproveIntegral).
 *
 * @param App   2D array (Ncurves × Nom): positive “importance” curves evaluated on the current mesh om.
 * @param dh    1D array (Nom): quadrature weights associated with om (e.g. trapezoid/Newton–Cotes weights).
 * @param alpha Exponent controlling mesh concentration: alpha>1 concentrates more strongly around peaks.
 * @param eps   Small floor value to avoid zero/negative weights (prevents log/zero issues downstream).
 *
 * @return fw2  1D array (Nom): global mesh-weight function.
 *
 * Assumptions:
 *   - App is nonnegative or mostly nonnegative. If a curve integrates to zero (Z_p=0), it is skipped.
 *   - dh corresponds to the same mesh points as App.
 */  
inline bl::Array<double,1>
BuildFw2FromApp(const bl::Array<double,2>& App,    // (Ncurves, Nom)
                const bl::Array<double,1>& dh, double alpha, double eps=1e-10){
  
  const int Ncurves = App.extent(0);
  const int Nom = App.extent(1);

  bl::Array<double,1> fw0(Nom); fw0 = 0.0;

  for (int p=0; p<Ncurves; ++p){
    bl::Array<double,1> p1(App(p, bl::Range::all()) * dh);
    double Z = bl::sum(p1);
    if (Z > 0) fw0 += p1 / Z;
  }
  for (int j=0; j<Nom; ++j) if (fw0(j) <= 0) fw0(j) = eps;

  // smooth fw0 -> fw1
  bl::Array<double,1> fw1(Nom);
  fw1(0) = 0.5*(fw0(0)+fw0(1));
  for (int j=1; j<Nom-1; ++j) fw1(j) = (fw0(j-1)+fw0(j)+fw0(j+1))/3.0;
  fw1(Nom-1) = 0.5*(fw0(Nom-1)+fw0(Nom-2));

  // sharpen: fw2 = fw1^alpha
  bl::Array<double,1> fw2(bl::pow(fw1, alpha));

  return fw2;
}

/**
 * @brief Find index range [ia, ib] on a monotone frequency grid om such that om(ia) ≲ La and om(ib) ≳ Lb.
 *
 * The adaptive mesh is often built in a central frequency window [La, Lb] (where structure is expected),
 * with simpler “tail” meshes outside. This helper locates that window robustly.
 *
 * @param om  1D array (Nom): strictly increasing frequency grid.
 * @param La  Left boundary of the “central” window.
 * @param Lb  Right boundary of the “central” window.
 * @param ia  (output) left index (clamped into [0, Nom-1]).
 * @param ib  (output) right index (clamped into [0, Nom-1]).
 *
 * Notes:
 *   - The routine clamps ia and ib to valid bounds.
 *   - This function assumes om is increasing; behavior is undefined if om is not sorted.
 */
inline std::pair<int,int> FindWindowIndices(const bl::Array<double,1>& om, double La, double Lb){
  const int N = om.extent(0);
  
  int ia = 0;
  while (ia < N && om(ia) < La) ++ia;
  ia = std::max(0, ia - 1);
  
  int ib = 0;
  while (ib < N && om(ib) < Lb) ++ib;
  ib = std::min(N - 1, ib);
  
  return {ia, ib};
}

/**
 * @brief Extract a sub-window of (om, fw2) and remove near-duplicate ω points.
 *
 * Some upstream grid generators (or floating-point issues) can produce nearly coincident points.
 * Such near-duplicates can destabilize interpolation and equal-weight grid generation. This routine
 * extracts the region [ia, ib] and removes points that are too close to the previously kept point.
 *
 * @param om      1D array (Nom): original frequency grid.
 * @param fw2     1D array (Nom): mesh-weight function evaluated on om.
 * @param ia      left index (inclusive).
 * @param ib      right index (inclusive).
 * @param om_out  (output) extracted ω values after deduplication.
 * @param fw_out  (output) extracted fw2 values after deduplication.
 * @param rel     relative spacing threshold used for deduplication.
 *
 * The deduplication criterion is applied to consecutive kept points:
 *   keep ω_k only if (ω_k - ω_prev) > rel * max(1, |ω_k|).
 *
 * Assumptions:
 *   - om is strictly increasing on [ia, ib].
 *   - ia, ib are valid and ia<=ib.
 */  
inline void ExtractWindowDedup(bl::Array<double,1>& om_out, bl::Array<double,1>& fw_out,
			       const bl::Array<double,1>& om, const bl::Array<double,1>& fw2, int ia, int ib, double rel=1e-6){
  const int N0 = ib - ia + 1;
  om_out.resize(N0);
  fw_out.resize(N0);

  int i_insert = 0;
  int last_kept = -1;
  for (int i=0; i<N0; ++i){ // here we remove some points which accidentally collapse on the same point
    int idx = ia + i;
    if (last_kept >= 0){
      double dx = om(idx) - om(last_kept);
      if (dx < rel * om(idx)*om(idx)) continue;
      //if (dx < rel * std::max(1.0, std::abs(om(idx)))) continue;
    }
    om_out(i_insert) = om(idx);
    fw_out(i_insert) = fw2(idx);
    last_kept = idx;
    ++i_insert;
  }
  om_out.resizeAndPreserve(i_insert);
  fw_out.resizeAndPreserve(i_insert);
}
/**
 * @brief Build an adaptive frequency mesh xn from a global weight fw2(ω) on an existing grid om.
 *
 * This routine encapsulates the full “adaptive mesh construction” used in the solver:
 *   - build a dense central mesh in [La, Lb] using equal-weight intervals of fw2(ω),
 *   - attach simpler tail meshes outside that region,
 *   - remove near-duplicates,
 *   - optionally refine by inserting midpoints that significantly change a target integral.
 *
 * Typical workflow:
 *   1) Compute fw2(ω) from App using BuildFw2FromApp().
 *   2) Call BuildAdaptiveMeshFromFw2() to get xn.
 *   3) Interpolate/resample all needed quantities from om to xn.
 *   4) Replace om ← xn and recompute dh on the new mesh.
 *
 * @param om      1D array: current frequency grid (strictly increasing).
 * @param fw2     1D array: weight function on om (positive).
 * @param App     2D array: importance curves used by AddMidpointsToImproveIntegral() to estimate error.
 * @param La,Lb   central window boundaries for the “nonuniform” grid generation.
 * @param Nw_1    number of points reserved for both tails combined (split by tail lengths).
 * @param Nw_2    target number of points in the central region (may change after refinement).
 * @param small   refinement tolerance for AddMidpointsToImproveIntegral().
 * @param Nw_max  maximum number of refinement rounds or maximum points (depending on your implementation).
 * @param Nw_itt  maximum iterations in midpoint refinement.
 * @param clog    output stream for logging decisions (added points, removed duplicates, etc.).
 *
 * @return xn     1D array: new adaptive mesh (strictly increasing).
 *
 * Important notes:
 *   - This routine does NOT broadcast xn; in MPI, you should call it only on the master rank and broadcast.
 *   - It assumes om spans a sufficiently wide interval containing [La, Lb].
 *   - It assumes fw2 is strictly positive (floored at eps upstream).
 */
inline bl::Array<double,1> BuildAdaptiveMeshFromFw2(const bl::Array<double,1>& om, const bl::Array<double,1>& fw2, const bl::Array<double,2>& App, // for AddMidpointsToImproveIntegral
						    double La, double Lb, int Nw_1, int Nw_2, double small, int Nw_max, int Nw_itt, std::ostream& clog, bool debug=false){
  auto [ia,ib] = FindWindowIndices(om, La, Lb);
  
  bl::Array<double,1> omw, fww;
  ExtractWindowDedup(omw, fww,
		     om, fw2, ia, ib, 1e-6);

  // central grid
  bl::Array<double,1> xn0 = NewGrid(fww, omw, Nw_2);
  double La_ = xn0(0), Lb_ = xn0(xn0.extent(0)-1);
  
  int N = om.extent(0);
  double dL = (La_ - om(0)) + (om(N-1) - Lb_);
  if (dL <= 0) {
    // fall back: just use xn0
    bl::Array<double,1> xn = xn0.copy();
    AddMidpointsToImproveIntegral(xn, App, om, clog, small, Nw_max, Nw_itt);
    return xn;
  }

  int Nw1a = int(Nw_1 * (La_ - om(0))/dL + 0.5);
  int Nw1b = int(Nw_1 * (om(N-1) - Lb_)/dL + 0.5);
  Nw1a = std::max(1, Nw1a);
  Nw1b = std::max(1, Nw1b);
  double dh_uniform = dL / (Nw1a + Nw1b);

  // possibly adding one point at the beginning of xn0 and at the end of xn0 mesh when too big
  // differences in points separation occurs at the mesh merging.
  bl::TinyVector<double,2> extra_w(0,0);
  if (xn0.extent(0) >= 2 && (xn0(1) - xn0(0) > 1.5*dh_uniform))
    extra_w[0] = 0.5*(xn0(0)+xn0(1));
  if (xn0.extent(0) >= 2 && (xn0(xn0.extent(0)-1) - xn0(xn0.extent(0)-2) > 1.5*dh_uniform))
    extra_w[1] = 0.5*(xn0(xn0.extent(0)-1)+xn0(xn0.extent(0)-2));

  int size = xn0.extent(0) + Nw1a + Nw1b + (extra_w[0]!=0) + (extra_w[1]!=0);
  bl::Array<double,1> xn(size);
  // left tail
  int offset=0;
  for (int i=0; i<Nw1a; ++i) xn(offset++) = om(0) + (i+0.0)/Nw1a*(La_ - om(0));
  xn(offset++) = xn0(0);
  if (extra_w[0]!=0) xn(offset++) = extra_w[0];
  for (int i=1; i<xn0.extent(0)-1; ++i) xn(offset++) = xn0(i);
  if (extra_w[1]!=0) xn(offset++) = extra_w[1];
  xn(offset++) = xn0(xn0.extent(0)-1);
  // right tail
  for (int i=0; i<Nw1b; ++i) xn(offset++) = Lb_ + (i+1.0)/Nw1b*(om(N-1) - Lb_);
  
  // final dedup pass (simple)
  std::deque<int> remove;
  for (int i=xn.extent(0)-1; i>1; --i)
    if (xn(i) - xn(i-1) < 1e-6 * xn(i)*xn(i)) remove.push_back(i);
  
  if (!remove.empty()){
    for (int ir : remove)
      for (int j=ir+1; j<xn.extent(0); ++j) xn(j-1) = xn(j);
    xn.resizeAndPreserve(xn.extent(0) - (int)remove.size());
  }
  
  if (debug){
    ofstream dbg2("cdebug.2"); dbg2.precision(14);
    dbg2<<"# xn["<<Nw1a-1<<"]="<<xn(Nw1a-1)<<" xn["<<Nw1a<<"]="<<xn(Nw1a)<<" xn["<<Nw1a+1<<"]="<<xn(Nw1a+1)<<" ";
    dbg2<<"xn["<<Nw1a+Nw_2-1<<"]="<<xn(Nw1a+Nw_2-1)<<" xn["<<Nw1a+Nw_2<<"]="<<xn(Nw1a+Nw_2)<<" xn["<<Nw1a+Nw_2+1<<"]="<<xn(Nw1a+Nw_2+1)<<endl;
    for (int i=0; i<xn.size(); ++i)
      dbg2<<i<<" "<<xn(i)<<endl;
    
    ofstream dbg3("cdebug.3"); dbg3.precision(14);
    dbg3<<"# xn0["<<0<<"]="<<xn0(0)<<" xn["<<xn0.size()-1<<"]="<<xn0(xn0.size()-1)<<endl;
    for (int i=0; i<xn0.size(); ++i)
      dbg3<<i<<" "<<xn0(i)<<endl;
  }

  // Here we add points in the midpoint between two successive points on mesh xn. If the integral is changed for more than
  // amount small, we accept the point. We then repeat the process up to Nwmax=20 times,
  // but we stop when no midpoint point changes the integral for more than 1e-5=small, or, we exceed number of points 1000=Nwmax.
  AddMidpointsToImproveIntegral(xn, App, om, clog, small, Nw_max, Nw_itt);
  return xn;
}

/**
 * @brief Resample a 2D complex array A(k,ω) from an old mesh to a new mesh using linear interpolation.
 *
 * This is the “plumbing” routine used after a new mesh xn is built. It takes a set of K complex curves
 * sampled on om_old and interpolates them onto om_new.
 *
 * The interpolation uses InterpLeft(x, pos, om_old) to compute an interval index and fraction,
 * and linInterp(Aold, k, ip) (or equivalent) to evaluate Aold(k, x).
 *
 * @param om_old  1D array: original mesh.
 * @param om_new  1D array: new mesh (monotone increasing).
 * @param Aold    2D complex array (K × Nom_old).
 * @param Anew    (output) 2D complex array resized to (K × Nom_new).
 *
 * Assumptions:
 *   - om_old and om_new are increasing.
 *   - om_new lies within [om_old(0), om_old(end)] or InterpLeft/linInterp must clamp.
 *
 * Performance:
 *   - O(K * Nom_new). Usually dominated by solver work, but can be significant for very large K.
 *   - If you need higher-order interpolation, replace linInterp with a spline evaluator.
 */
inline void ResampleComplex2D(const bl::Array<double,1>& om_old, const bl::Array<double,1>& om_new,
                              const bl::Array<std::complex<double>,2>& Aold, // (K, Nom_old)
                              bl::Array<std::complex<double>,2>& Anew){      // (K, Nom_new)
  const int K = Aold.extent(0);
  Anew.resize(K, om_new.extent(0));

  int pos = 0;
  for (int iw=0; iw<om_new.extent(0); ++iw){
    intpar ip = InterpLeft(om_new(iw), pos, om_old);
    for (int k=0; k<K; ++k)
      Anew(k,iw) = linInterp(Aold, k, ip);
  }
}

/**
 * @brief Build a Fermi-function array few(ω) on a given frequency mesh.
 *
 * In several routines we use a phenomenological factor (1 - few(ω)) to suppress spectral weight
 * above the Fermi level (or to incorporate finite-T occupation effects).
 *
 * @param om   1D array: frequency mesh.
 * @param T    temperature (in the same energy units as om).
 * @param few  (output) array resized to om.size() containing few(i) = fermi(om(i)/T).
 *
 * Note:
 *   - fermi(x) is assumed to be a stable implementation (no overflow for large |x|).
 */  
inline void BuildFewOnMesh(const bl::Array<double,1>& om, double T, bl::Array<double,1>& few){
  few.resize(om.extent(0));
  for (int i=0; i<om.extent(0); ++i) few(i) = fermi(om(i)/T);
}

/*
int main()
{
  cout.precision(15);
  double small = 1e-5;
  double alpha = 0.4;
  int Nwstart = 300;
  int Nwmax = 1000;
  int Nitt = 20;
  
  int Npseudo=2;
  auto [om,App] = ReadSigma("Sigma.000~", Npseudo);
  
  Array<double,1> dh = GiveDh(om);
  Array<double,1> integral0(Npseudo);
  for (int i=0; i<Npseudo; ++i){
    integral0(i) = sum(App(i,Range::all())*dh);
    cout <<"# Integral["<<i<<"]="<< integral0(i) << endl;
  }

  vector<Spline1D<double>> fx = SplineSimpleFunctions(App, om);

  // construct combined function used in constructing new grid. This allows a lot of flexibility.
  Array<double,1> fw1(om.size()); fw1=0;
  for (int i=0; i<Npseudo; ++i){
    Array<double,1> p1(App(i,Range::all())*dh);
    fw1 += p1/sum(p1);
  }
  Array<double,1> fw2( pow(fw1,alpha) );

  // construct new grid, which is optimized for fw2 function on original om grid. The mesh has only Nwstart points.
  Array<double,1> xn = NewGrid(fw2, om, Nwstart);


  // Here we add points in the midpoint between two successive points on mesh xn. If the integral is changed for more than
  // amount small, we accept the point. We then repeat the process up to Nwmax=20 times,
  // but we stop when no midpoint point changes the integral for more than 1e-5=small, or, we exceed number of points 1000=Nwmax.
  AddMidpointsToImproveIntegral(xn, fx, om, small,Nwmax,Nitt);
  
  // Now recomputing fx on new mesh xn
  Array<double,2> fn(Npseudo,xn.size());
  int pos=0;
  for (int i=0; i<xn.size(); ++i){
    intpar ip=InterpLeft(xn(i),pos,om);
    for (int idim=0; idim<Npseudo; ++idim)
      fn(idim,i) = fx[idim].LinInterp(ip);
  }
  
  // Check quality of the integral
  Array<double,1> dhn = GiveDh(xn);
  for (int i=0; i<Npseudo; ++i){
    double intg = sum(fn(i,Range::all())*dhn);
    cout<<"# Final integral["<<i<<"]="<<intg<<" and original "<<integral0(i) << " with diff="<<intg-integral0(i)<<endl;
  }
  // Printing final result
  for (int i=0; i<xn.size(); ++i){
    cout << xn(i)<<"  ";
    for (int idim=0; idim<Npseudo; ++idim) cout<<fn(idim,i)<<" ";
    cout<<endl;
  }
  return 0;
}
*/
#endif
