#ifndef _nca_oca_imag_
#define _nca_oca_imag_

#include <blitz/array.h>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <climits>
#include <functional>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <numeric>
#include <type_traits>
#include "read_cix.h"
#include "oca_dlr_types.h"
using namespace blitz;

inline double intg_dh(int it, int Nlast)
{// returns coefficient for integrating a function with O(d^4) error, which has the following form:
  // \int_0^b dx f(x)
  // where we have Nlast+1 points, distributed at i/Nlast, with i=0....,Nlast
  // The distance between points is b/Nlast, hence dx=b/Nlast, except at the edge, where
  // trapezoid rule would reguest dx=0.5*b/Nlast.
  // For higher order precision, we need to use the coefficient such as those specified here to achive higher order
  // precision
  static const double coefN[]={17/48.,59/48.,43/48.,49./48.,1.0};
  static const double coef6[]={41./140, 216./140, 27./140, 272./140, 27./140, 216./140, 41./140};
  static const double coef5[]={95./288,125./96,125./144,125./144,125./96,95./288};
  static const double coef4[]={14./45,64./45,24./45,64./45,14./45};
  static const double coef3[]={3./8,9./8,9./8,3./8};
  static const double coef2[]={1./3,4./3,1./3};
  if (Nlast>=7){
    if (it<=3 or it >=Nlast-3)
      return (it<=3) ? coefN[it] : coefN[Nlast-it];
    else
      return 1;
  }else{
    switch(Nlast){
    case 0: return 0;
    case 1: return 0.5;
    case 2: return coef2[it];
    case 3: return coef3[it];
    case 4: return coef4[it];
    case 5: return coef5[it];
    case 6: return coef6[it];
    }
  }
  return 0;
}

double I_NCA(Array<double,2>& Sm, Array<double,2>& G00t, Array<double,2>& F00t, Array<double,2>& F01t,
	     bool cmpPhysG, double Q, double beta, const Array<double,1>& tau, const Array<double,1>& tau_dh, const Array<double,2>& Gm, const Array<double,2>& Delta, Cix& cix, MPIcls& mpi)
{ // ********** NCA on imaginary  ***************
  // output:
  //          Sm(Npseudo,Ntau)   -- pseudoparticle self-energy
  //          G00t(Nbath,Ntau)   -- local Green's function
  //          F00t(Nbath,Ntau)   -- the vertex function to be used for local Green's function
  //          F01t(Nbath,Ntau)   -- the second vertex function (destroys cooper pair or down)
  // input :  
  //          Q                 -- <Q> as pseudo charge
  //          tau(Ntau)         -- imaginary time mesh
  //          Gm(Npseudo,Ntau)  -- input pseudoparticle propagators
  //          Delta(Nbath,Ntau) -- input hybridization
  //          cix               -- class with information about the exact diagonalization of the atom from cix file.
  //          beta              -- inverse temperature
  //          cix.Fi(2*Nbath,  Npseudo+1) -- <m|psi^\dagger|m'>, and can be 0,+m,-m. Negative means state |m| with sign -1.
  //          cix.Fi(2*Nbath+1,Npseudo+1) -- <m|psi|m'>, and can be 0,+m,-m.
  //
  //      states[1]
  //      --->---
  //    /  D[s,t] \
  //   |----<------|
  //    \         /
  //      ---<---
  //      states[0]
  //     deg_pp=[1,2,1]
  bool debug=false;
  //const double coef[]={17/48.,59/48.,43/48.,49./48.,1.0};// O(dh^4) integration coefficients.
  //
  //  NCA-type diagrams have the form:
  //   Phi = G[ii](tau) * G[jj](beta-tau) * D[orb,tau] * sign(Fi[2*orb+1,j1])*sign(Fi[2*orb,i2])
  //            t=0
  //           j1  i1   where Fi[2*orb+1,j1]=i1 or psi_{orb}|j1>=|i1>
  //          -->----
  //        /    |    \   block ib
  // G[jj] | D[orb,t]  |  G[ii]=G(i1,i2)=G[ib](_i1_,_i2_)
  //       |     v     | 
  //block jb\    |    /
  //          ---<---
  //           j2  i2   where Fi[2*orb,i2]=j2 or psi^+_{orb}|i2>=|j2>
  //             t
  //  G[jj]=G(j2,j1)=G[jb](_j2_,_j1_)
  //
  //    G[jj]=G(j2,j1)=G[jb](_j2_,_j1_)
  //      --->---
  // j2 /  D[s,t] \ j1
  //   |----<------|
  // i2 \  iorb   / i1
  //      ---<---
  //    G[ii]=G(i1,i2)=G[ib](_i1_,_i2_)
  //
  int Nt = tau.extent(0);
  
  if (Gm.extent(1)!=Nt){ (*mpi.out) << "ERROR in nca_oca_imag_axis:I_NCA Gm.extent(1)="<<Gm.extent(1)<<"!="<<Nt<<endl;}
  if (Delta.extent(1)!=Nt) { (*mpi.out) << "ERROR in nca_oca_imag_axis:I_NCA Delta.extent(1)="<<Delta.extent(1)<<"!="<<Nt<<endl; }
  for (int it=0; it<Nt/2; ++it)
    if (fabs(tau(it)+tau(Nt-1-it)-beta)>1e-5){ (*mpi.out)<<"ERROR in nca_oca_imag_axis:I_NCA tau mesh should be symmetric around beta/2"<<endl;}
  
  map<TinyVector<int,3>,int,Class1Compare> states_unique;
  deque<TinyVector<double,3>> mm_nds;
  // mm_nds[i,0] = mm,
  // mm_nds[i,1] = nd_end/normal_up_or_create_cooper_pair
  // mm_nds[i,2] = nd_end/normal_dn_or_destroys_cooper_pair
  for (int ib=0; ib<cix.Nblock; ++ib){          // over all blocks of PseudoParticles
    int Ni = cix.Block[ib].size();              // dimension of the block
    for (int iorb=0; iorb<cix.Nbath; ++iorb){   // over all baths for back propagating NCA diagram
      int s = abs(cix.bath_sigind0(iorb));            // which unique bath is that
      double s_sgn = cix.bath_sigind0(iorb) < 0 ? -1 : 1;
      // Our assumed convention in normal state is that Operators[2*i] is a creation operator and Operators[2*i+1] is corresponding anhilation operator
      // For bath k, Fi_indx[2*k] and Fi_indx[2*k+1] contains its first (usually creation) and second (usually anhilation) operator.
      // When both Fi_indx[2*k] and Fi_indx[2*k+1] are odd, we have two anhilation operators, hence this destroys Cooper pair.
      // For 2x2 plaquette, this should happen for iorb=9 and iorb=11
      bool destroys_Cooper_pair = cix.Fi_indx(2*iorb)%2==1 and cix.Fi_indx(2*iorb+1)%2==1;
      // normal_down is true when one of the operators is creation and one anhilation (Fi_indx[2*k]+Fi_indx[2*k+1] is odd) and spin of that orbital is negative (-1/2), and the orbital is part of the SC block.
      // normal_down is hence true for 22 block when we have superconductivity. It is false in normal state.
      bool normal_down = cix.SC_bath(iorb) and ((cix.Fi_indx(2*iorb)+cix.Fi_indx(2*iorb+1))%2==1 and cix.bath_spin(iorb)<0);
      // Here jb is block after applying operator psi_{iorb}^+, i.e., |jb> = F^+_{iorb}|ib>
      int jb = cix.Fbc(2*iorb,ib);              // if we apply F^+_{iorb}|ib> = |jb> . Here |jb> is block, not a state. Notice that for iorb=9 or 11 we actually have anhilation operator
      if (jb<0) continue;                       // Pauli principle sometimes makes this vanish
      int Nj = cix.Block[jb].size();            // size of the |jb> block
      //
      //int njmc = cix.Njm_c(ib,iorb); // index to identify whether njm is proportional to F or not.
      NjmResult nx = cix.Njm_ind(ib,2*iorb);    // always replace the end of hybridization (psi^+ in normal state) with the njm operator. In SC state it could be anhilation operator, therefore we need Njm_ind function.
      // Njm_ind function gives Njm(ib,2*iorb) when psi^+ is creation operator, but when psi is anhilation, we need to take different Njm and nx.transposed becomes true.
      Array<double,2> njm;           // matrix when njm is not proportional to F
      double njmr = nx.njm_r;        // proportionality constant for njm
      if (nx.njm_c==2){              // njm is not proportional
	if (not nx.transposed){
	  njm.resize(nx.njm->extent(0),nx.njm->extent(1));
	  njm = *nx.njm;
	}else{
	  njm.resize(nx.njm->extent(1),nx.njm->extent(0));
	  for (int n1=0; n1<njm.extent(0); ++n1)
	    for (int n2=0; n2<njm.extent(1); ++n2)
	      njm(n1,n2) = (*nx.njm)(n2,n1); // here we just transpose, because njm is real. But in general it should be dagger, i.e., also conjugated.
	  // note that this is for the case the first operator is anhilation and not creation operator.
	}
      }
      //  NCA-type diagrams have the form:
      //   Phi = G[ii](tau) * G[jj](beta-tau) * D[orb,tau] * sign(Fi[2*orb+1,j1])*sign(Fi[2*orb,i2])
      //            t=0
      //           j1  i1   where Fi[2*orb+1,j1]=i1 or psi^{start}_{orb}|j1>=|i1>
      //          -->----
      //        /    |    \   block ib
      // G[jj] | D[orb,t]  |  G[ii]=G(i1,i2)=G[ib](_i1_,_i2_)
      //       |     v     | 
      //block jb\    |    /
      //          ---<---
      //           j2  i2   where Fi[2*orb,i2]=j2 or psi^{end}_{orb}|i2>=|j2>
      //             t
      //  G[jj]=G(j2,j1)=G[jb](_j2_,_j1_)
      //
      // |jb> = psi^{end}_{iorb}|ib>  , here {end} is \dagger for normal state, and is anhilation operator for iorb=9,11
      // we first apply psi_{iorb}|j1>
      for (int _j1_=0; _j1_<Nj; ++_j1_){ // over matrix corresponding to block jb
	int j1 = cix.Block[jb](_j1_);               // |j1> between [0,...,Nj]
	int i1 = cix.Fi(2*iorb+1, j1);              // |i1> = psi^{start}_{iorb}|j1>
	if (i1==0) continue;
	int _i1_ = cix.inBlock(abs(i1));
	for (int _i2_=0; _i2_<Ni; ++_i2_){
	  int i2 = cix.Block[ib](_i2_);
	  int j2 = cix.Fi(2*iorb, i2);                // |j2> = psi^{end}_{iorb}|i2>
	  int _ii_ = cix.Sgm_symmetry[ib](_i1_,_i2_); // symmetry of propagator G(i1,i2)==G(_ii_). Here _ii_ starts with 1 (if nonzero) and can be negative (when propagator is equal but opposite to some other propagator).
	  if (_ii_==0) continue;                      // this G(_ii_) vanishes.
	  int ii = abs(_ii_)-1;                       // actual place where G(ii) is found.
	  if (j2!=0){                                 // if j2=0, diagram vanishes.
	    int jb2 = cix.iBlck(abs(j1))-1;         // jb is block in which j1 and j2 are located. Could be obtained by cix.Fbc(2*iorb,ib+1)
	    if (jb2!=jb) cout<<"ERROR jb="<<jb<<" jb2="<<jb2<<endl;
	    int _j2_ = cix.inBlock(abs(j2));          // j2 is also counted from the empty state, while _j2_ is index within the block jb
	    int _jj_ = cix.Sgm_symmetry[jb](_j2_,_j1_); // G(j1,j2)=G(_jj_)
	    if (_jj_!=0){                               // G(_jj_) can sometimes vanish
	      int jj = abs(_jj_)-1;                     // index for G(_jj_) == G[jj]
	      TinyVector<int,3> st(s,jj,ii);            // store ii,jj,s in such a way that it is sorted with respect to ii, and after that jj, and finally s.
	      //           s
	      //      ii   |  jj
	      //    ------ v -----
	      auto [it, inserted] = states_unique.try_emplace(st, mm_nds.size()); // try to insert st as a new key and mm_nds.size() as the value.
	                                                // If key is not yet present, inserted=true and "it" points to its current place. If key is present, "it" points to its existing place.
	      int iind = it->second;                    // The key st has value iind, either just inserted, or found from before
	      if (inserted) mm_nds.push_back(TinyVector<double,3>(0, 0, 0)); // if it was just inserted, deque needs to have this entry set to zero.
	      //              sign of F^+ and F operator and sign of PP propagator can also be negative compared to what is being stored.
	      double mm = sign(static_cast<double>(i1)*static_cast<double>(j2)*static_cast<double>(_ii_)*static_cast<double>(_jj_))*s_sgn; // actual matrix element for NCA diagram.
	      // Here sign(i1)*sign(j2) is matrix element due to F^+ and F operators. The sign of _ii_ and _jj_ comes from the fact that some propagators have opposite signs as the value in G(_ii_) or G(_jj_).
	      if (DEBUG) cout<<"1) ib="<<ib<<" iorb="<<iorb<<" (i1,i2)=("<<i1<<","<<i2<<") (j1,j2)=("<<j1<<","<<j2<<") _ii_="<<_ii_<<" _jj_="<<_jj_<<" st="<<st<<" mm="<<mm<<endl;
	      double nds = 0;
	      if (nx.njm_c==1 and cmpPhysG){  // if njm is proportional to F
		nds=mm*njmr; // mm contains all other terms, we just add njmr, which is the ratio.
	      }
	      if (destroys_Cooper_pair or normal_down)   // this is true for SC state and G_{12} and G_{22} (bath 9, 11, is for G_{12} and 3, 5 for G_{22})
		mm_nds[iind] += TinyVector<double,3>(mm, 0, nds);
	      else                                       // this is always true for normal states, and in SC state for G_{11} and G_{21} (bath 8 and 10)
		mm_nds[iind] += TinyVector<double,3>(mm, nds, 0);
	    }
	  }
	  if (nx.njm_c==2 and cmpPhysG){ // if njmc is not proportional, but completely different matrix, we need to go over all matrix elements.
	    // We are inserting Njm instead of psi^+_{iorb}|i2>=|j2> operator. But now |i2> and |j2> have to be unrestricted, not connected by psi^+_{iorb} but Njm[i2,j2]
	    // Notice that we keep F^+_{iorb}|i1>=|j1> on the other end.
	    //           s
	    //    ii     |      jj
	    // --------- v ----------
	    // i1   i2  njm  j2    j1
	    //   G(ii)*njm*i2,j2)*G(jj)
	    //G(i1,i2)*njm(i2,j2)*G(j2,j1)
	    for (int _j2_=0; _j2_<Nj; ++_j2_){// _j2_ is now arbitrary index in block jb (size Nj), and not given by F^+_{iorb}|i2>!=|j2>
	      bool njm_finite = abs(njm(_i2_,_j2_))>1e-10; // is njm finite?
	      int _jj_ = cix.Sgm_symmetry[jb](_j2_,_j1_);  // for this general j2, what is G(j2,j1)=G[jj] and is it finite?
	      if (_jj_!=0 and njm_finite){
		double nds = njm(_i2_,_j2_)*sign(_ii_*_jj_)*sign(i1)*s_sgn;    // we here use njm(i2,j2) 
		int j2=cix.Block[jb](_j2_);                                    // what is actual value of j2 in index which starts from the empty state.
		int jj = abs(_jj_)-1;
		TinyVector<int,3> st(s,jj,ii);        // F00[s] has contribution G[ii]*G[jj]*nds
		auto [it, inserted] = states_unique.try_emplace(st, mm_nds.size()); // check if st already exist? If yes, "it" contains its place. If not, it is inserted now and "it" contains its place.
		TinyVector<double,3> nds_val = (destroys_Cooper_pair or normal_down) ? TinyVector<double,3>(0, 0, nds) : TinyVector<double,3> (0, nds, 0);
		int iind = it->second;                // this is the value associated with key st. It contains position where nds for (ii,jj,s) is located.
		if (inserted)
		  mm_nds.push_back(nds_val); // was just now inserted, hence need an empty value in mm_nds at the back of deque.
		else
		  mm_nds[iind] += nds_val; // the value is njm(i2,j2) and also sign(i1) for the existing F^+|i1> and possible negative sign of G[ii] and G[jj] propagators.
		//cout<<"1i) _ii_="<<_ii_<<" _jj_="<<_jj_<<" i1="<<i1<<" i2="<<i2<<" j1="<<j1<<" j2="<<j2<<" _j1_="<<_j1_<<" _j2_="<<_j2_<<" st="<<st<<" iorb="<<iorb<<" nds="<<nds<<endl;
	      }
	    }
	  }
	}
      }
    }
  }
  if (debug){
    if (mpi.rank==mpi.master){
      const int Nd=2;
      int ijk=0;
      for (const auto& [st, iind] : states_unique){
	double mm = mm_nds[iind][0];
	double nds = mm_nds[iind][1];
	int s = st[0];
	TinyVector<int,Nd> states(st[2],st[1]);
	
	double maxv0 = max(Gm(states[0],Range::all()));
	double maxv1 = max(Gm(states[1],Range::all()));
      
	fmt::print(cout," {} ii0={:3d},ii1={:3d},s={:d}: mm={} nds={} deg_pp0={} maxv0={} maxv1={}", ijk,states[0]+1,states[1]+1,s,mm,nds,cix.deg_pp(states[0]), maxv0, maxv1);
	ijk+=1;
      
	fmt::print(cout,"[");
	for (int l=0; l<cix.Sgm_indx[states[0]].size(); ++l){
	  auto v = cix.Sgm_indx[states[0]][l];
	  int ib = abs(v[0]);
	  fmt::print(cout, "({},{}),", cix.Block[ib](v[1]), cix.Block[ib](v[2]) );
	}
	fmt::print(cout,"] [");
	for (int l=0; l<cix.Sgm_indx[states[1]].size(); ++l){
	  auto v = cix.Sgm_indx[states[1]][l];
	  int ib = abs(v[0]);
	  fmt::print(cout, "({},{}),", cix.Block[ib](v[1]), cix.Block[ib](v[2]) );
	}
	fmt::print(cout,"]\n");
      }
    }
  }
  
  //auto t0 = high_resolution_clock::now();
  auto t0 = get_cpu_time();
  Sm.resize(cix.Npseudo_unique,Nt); Sm=0;   // output PP self-energy
  Array<double,2> G00(cix.Nbath_unique_unique,Nt); G00=0; // physical Green's function
  Array<double,2> F00(cix.Nbath_unique_unique,Nt); F00=0; // physical vertex function
  Array<double,2> F01(cix.Nbath_unique_unique,Nt); F01=0; // physical vertex function
  //double dh = beta/(Nt-1); // for integration on the uniform mesh
  double Phi_tot=0;        // Phi
  for (const auto& [st, iind] : states_unique){
    int s = st[0];    // bath index as found in bath_unique
    int ip1 = st[1];  // called jj (of G[jj]) above
    int ip0 = st[2];  // called ii (of G[ii]) above
    double mm   = mm_nds[iind][0];// matrix element for Phi or PPSig or G00 
    double nds1 = mm_nds[iind][1];// matrix element for F00, which comes from Njm
    double nds2 = mm_nds[iind][2];// matrix element for F01, which comes from Njm
    for (int it=0; it<Nt; ++it){ // over non-uniform time mesh, which is symmetric around beta/2, such that tau[Nt-i-1]==beta-tau[i].
      double G1 = Gm(ip0,it);     // Gt(tau)
      double G2 = Gm(ip1,Nt-1-it);// G2(beta-tau)
      if (mm!=0){
	double Dl = Delta(s,it);
	Sm(ip0,Nt-1-it) += mm * Dl * G2; // Sm(ii, beta-tau) = d Phi/ d G[ii](tau) = G(jj,beta-tau) * Delta_s(tau) * mm
	Sm(ip1,it)      += mm * Dl * G1; // Sm(jj, tau) =  d Phi / d G2[jj](beta-tau) = G(ii,tau) * Delta_s(tau) * mm
	G00(s,Nt-1-it)  += mm * G1 * G2; // G{physica}(s,beta-tau) = G(ii,tau)*G(jj,beta-tau) * mm
	//double dx=dh;
	//if (it<4 or it >=Nt-4) dx *= (it<4) ? coef[it] : coef[Nt-1-it];
	//Phi_tot += G1*G2*Dl*mm*dx; // integral(Phi)
	Phi_tot += G1*G2*Dl*mm*tau_dh(it); // integral(Phi)
      }
      if (nds1!=0) F00(s,Nt-1-it) += nds1 * G1*G2; // F0 for up or create cooper pair
      if (nds2!=0) F01(s,Nt-1-it) += nds2 * G1*G2; // F01 for down or destroy cooper pair
    }
  }
  for (int ip=0; ip<cix.Npseudo_unique; ++ip)
    Sm(ip,Range::all()) /= cix.deg_pp(ip);

  G00t.resize(cix.Nbath_unique,Nt); G00t=0;
  F00t.resize(cix.Nbath_unique,Nt); F00t=0;
  F01t.resize(cix.Nbath_unique,Nt); F01t=0;
  if (cmpPhysG){
    Array<int,1> sigind(cix.Nbath_unique);
    for (int iorb=0; iorb<cix.Nbath; ++iorb) sigind(cix.bath_sigind(iorb)) = cix.bath_sigind0(iorb);
    Array<double,1> deg_bth(cix.Nbath_unique_unique); deg_bth=0; // degeneracy of each bath
    for (int iorb=0; iorb<cix.Nbath; ++iorb) deg_bth(abs(cix.bath_sigind0(iorb)))+= 1;
    Array<bool,1> bath_SC(cix.Nbath_unique); bath_SC=false;
    for (int iorb=0; iorb<cix.Nbath; ++iorb) bath_SC(cix.bath_sigind(iorb)) = cix.SC_bath(iorb);
    
    for (int ib=0; ib<cix.Nbath_unique; ++ib){
      int s = abs(sigind(ib));
      double s_sgn = sigind(ib) < 0 ? -1 : 1;
      double extra_2 = bath_SC(ib) ? 2.0 : 1.0;
      G00t(ib,Range::all()) = G00(s,Range::all())*(s_sgn/(deg_bth(s)*Q));
      F00t(ib,Range::all()) = F00(s,Range::all())*(extra_2*s_sgn/(deg_bth(s)*Q));
      F01t(ib,Range::all()) = F01(s,Range::all())*(extra_2*s_sgn/(deg_bth(s)*Q));
    }
  }
  Phi_tot /= Q;
  //auto t1 = high_resolution_clock::now();
  auto t1 = get_cpu_time();
  //fsec t_nca = t1-t0;
  double t_nca = t1-t0;
  fmt::print(*mpi.out,"  # t_nca={:.4f}s  Phi_NCA={}\n", t_nca, Phi_tot);
  return Phi_tot;
}

// Dense post-construction representation of one fully reduced OCA path.
// unordered_map is very useful while discovering and accumulating paths, but
// retaining its buckets/nodes wastes memory and gives poor traversal locality.
// This record is also the natural basis for the forthcoming MPI wire format.
struct OCACompactPath {
  OCAKey key{};                       // {s1,s0,p3,p2,p1,p0}
  std::array<double,5> weights{};    // {mm,nds1,nds0,anomalous1,anomalous0}
  std::array<int,4> inner_index{};   // legacy precomputed convolution indices
};

// Everything produced by the CIX/operator walk and needed by either OCA
// integrator.  The enumerator contains no time mesh, splines, AAA, or DLR
// objects: it is a purely discrete description of the allowed diagrams.
struct OCAEnumeration {
  // The enumeration is reusable only while every discrete input which can
  // remove or add a path remains unchanged.  Keeping this metadata beside
  // the large tables makes cache validation explicit and cheap.
  bool includes_physical=false;
  int npseudo=0;
  int nbath_unique=0;
  int cut_sum=0;
  vector<int> importance;
  // Replicated mode stores every reduced path on every rank.  Distributed
  // mode stores only the reduced bath-pair shards owned by this rank; the
  // compact DLR plans and their numerical contractions remain sharded too.
  bool paths_distributed=false;
  // A distributed topology may be replicated inside small MPI worker groups.
  // `topology_shard_count` is the number of disjoint path sets globally;
  // ranks with equal rank%topology_shard_count hold the same set and divide
  // its AAA poles.  Group-rank zero is the unique owner for disk-cache output
  // and global plan-count diagnostics.
  int topology_shard_count=1;
  int topology_shard_id=0;
  int topology_group_rank=0;
  int topology_group_size=1;
  unsigned long long global_nonzero_path_count=0;
  double estimated_replicated_megabytes=0.0;

  bool matches(bool cmpPhysG, const Cix& cix,
               const vector<int>& G_importance,
               int OCA_cut_sum) const
  {
    return includes_physical==cmpPhysG &&
      npseudo==cix.Npseudo_unique &&
      nbath_unique==cix.Nbath_unique_unique &&
      cut_sum==OCA_cut_sum && importance==G_importance;
  }

  // Sorted contiguous paths replace the construction hash table and its two
  // parallel arrays.  Sorting makes diagnostics and future MPI reductions
  // deterministic while improving the legacy integration traversal.
  vector<OCACompactPath> compact_paths;

  // Legacy uniform-mesh integration reuses four kinds of inner pseudo pair.
  vector<vector<TinyVector<int,4>>> unique_inner_states;

  // Lightweight representation consumed by oca_dlr.h.  It is populated only
  // when DLR production or diagnostics request it, avoiding duplicate storage
  // in a normal legacy-only run.
  vector<OCAPathRecord> paths;
};

// Stable signature for the optional on-disk path cache.  The cache stores
// only the mesh-independent OCAPathRecord stream, never Eigen sparse-matrix
// internals, so it remains portable across MPI sizes and plan layouts built
// by the same executable/architecture.
struct OCAPathCacheHeader {
  char magic[8]={'O','C','A','D','L','R','P','1'};
  std::uint32_t version=1;
  std::uint32_t record_size=sizeof(OCAPathRecord);
  std::uint64_t cix_hash=0;
  std::uint64_t importance_hash=0;
  std::uint64_t record_count=0;
  std::int32_t npseudo=0;
  std::int32_t nbath=0;
  std::int32_t cut_sum=0;
  std::int32_t includes_physical=0;
};
static_assert(std::is_trivially_copyable_v<OCAPathCacheHeader>);
static_assert(std::is_trivially_copyable_v<OCAPathRecord>);

inline std::uint64_t OCAFNV1a(const void* data,std::size_t bytes,
                             std::uint64_t hash=1469598103934665603ULL)
{
  const auto* p=static_cast<const unsigned char*>(data);
  for (std::size_t i=0;i<bytes;++i) {
    hash^=p[i];
    hash*=1099511628211ULL;
  }
  return hash;
}

inline std::uint64_t OCAFileHash(const std::string& filename)
{
  std::ifstream in(filename,std::ios::binary);
  if (!in) return 0;
  std::uint64_t hash=1469598103934665603ULL;
  std::array<char,1<<15> buffer{};
  while (in) {
    in.read(buffer.data(),buffer.size());
    const auto count=in.gcount();
    if (count>0) hash=OCAFNV1a(buffer.data(),static_cast<std::size_t>(count),hash);
  }
  return hash;
}

inline std::uint64_t OCAImportanceHash(const vector<int>& importance)
{
  return OCAFNV1a(importance.data(),importance.size()*sizeof(int));
}

// Choose the smallest number of disjoint topology shards whose estimated
// compiled plan fits the configured per-rank budget.  Every additional rank
// is then useful as a pole worker sharing one of those shards.  Explicit
// "replicated" and "distributed" retain their old meanings for controlled
// comparisons; only "auto" selects the hybrid two-dimensional layout.
inline int OCAChooseTopologyShardCount(const std::string& topology_mode,
                                       double estimated_replicated_megabytes,
                                       double max_replicated_megabytes,
                                       int world_size)
{
  if (world_size<=1 || topology_mode=="replicated") return 1;
  if (topology_mode=="distributed") return world_size;
  const double needed=std::ceil(
    estimated_replicated_megabytes/max_replicated_megabytes);
  return std::max(1,std::min(world_size,static_cast<int>(needed)));
}

inline void OCASetTopologyGroupMetadata(OCAEnumeration& result,MPIcls& mpi)
{
  const int shards=std::max(1,result.topology_shard_count);
  result.topology_shard_id=mpi.rank%shards;
  result.topology_group_rank=mpi.rank/shards;
  result.topology_group_size=1+(mpi.size-1-result.topology_shard_id)/shards;
}

#ifdef _MPI
// Convert the unique rank-local path stream produced by enumeration or a
// parallel cache read into the two-dimensional topology layout.  When there
// are at least as many nonempty reduced bath pairs as shards, an entire
// (s0,s1) pair remains on one shard, minimizing duplicate pseudo-pair catalogs.
// Only when more shards are required than bath pairs are individual records
// hash-partitioned, which is necessary to honor the memory bound.
inline void DistributeOCAPathRecordsHybridMPI(
    vector<OCAPathRecord>& paths,
    int nbath,
    int shard_count,
    bool already_pair_partitioned,
    MPIcls& mpi)
{
  if (mpi.size<=1 || shard_count<=1) return;
  shard_count=std::min(shard_count,mpi.size);

  vector<unsigned long long> local_pair_counts(nbath*nbath,0);
  for (const auto& path:paths) {
    if (path.s0<0 || path.s0>=nbath || path.s1<0 || path.s1>=nbath)
      throw std::out_of_range("OCA path bath index exceeds hybrid shard table");
    ++local_pair_counts[path.s0*nbath+path.s1];
  }
  vector<unsigned long long> global_pair_counts(local_pair_counts.size(),0);
  MPI_Allreduce(local_pair_counts.data(),global_pair_counts.data(),
                static_cast<int>(global_pair_counts.size()),
                MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  const int active_pairs=static_cast<int>(std::count_if(
    global_pair_counts.begin(),global_pair_counts.end(),
    [](unsigned long long count){return count!=0;}));

  // Fresh enumeration with one shard per rank and enough bath pairs already
  // has exactly this locality.  Avoid a second all-to-all in that common case.
  const bool keep_existing=already_pair_partitioned &&
    shard_count==mpi.size && shard_count<=active_pairs;
  if (!keep_existing) {
    vector<int> pair_shard(nbath*nbath,-1);
    if (shard_count<=active_pairs) {
      vector<int> order;
      for (int pair=0;pair<nbath*nbath;++pair)
        if (global_pair_counts[pair]!=0) order.push_back(pair);
      std::stable_sort(order.begin(),order.end(),[&](int a,int b) {
        return global_pair_counts[a]>global_pair_counts[b];
      });
      vector<unsigned long long> load(shard_count,0);
      for (int pair:order) {
        const int shard=static_cast<int>(
          std::min_element(load.begin(),load.end())-load.begin());
        pair_shard[pair]=shard;
        load[shard]+=global_pair_counts[pair];
      }
    }

    const auto destination_shard=[&](const OCAPathRecord& path) {
      if (shard_count<=active_pairs)
        return pair_shard[path.s0*nbath+path.s1];
      std::uint64_t hash=1469598103934665603ULL;
      hash=OCAFNV1a(path.pseudo.data(),path.pseudo.size()*sizeof(int),hash);
      hash=OCAFNV1a(&path.s0,sizeof(path.s0),hash);
      hash=OCAFNV1a(&path.s1,sizeof(path.s1),hash);
      return static_cast<int>(hash%shard_count);
    };

    vector<int> send_counts(mpi.size,0),recv_counts(mpi.size,0);
    for (const auto& path:paths) {
      const int target=destination_shard(path); // shard leader world rank
      if (send_counts[target]==INT_MAX)
        throw std::overflow_error("OCA hybrid send count exceeds INT_MAX");
      ++send_counts[target];
    }
    MPI_Alltoall(send_counts.data(),1,MPI_INT,
                 recv_counts.data(),1,MPI_INT,MPI_COMM_WORLD);
    vector<int> send_displacements(mpi.size),recv_displacements(mpi.size);
    long long send_total=0,recv_total=0;
    for (int rank=0;rank<mpi.size;++rank) {
      if (send_total>INT_MAX || recv_total>INT_MAX)
        throw std::overflow_error(
          "OCA hybrid redistribution displacement exceeds INT_MAX");
      send_displacements[rank]=static_cast<int>(send_total);
      recv_displacements[rank]=static_cast<int>(recv_total);
      send_total+=send_counts[rank];
      recv_total+=recv_counts[rank];
    }
    if (send_total>INT_MAX || recv_total>INT_MAX)
      throw std::overflow_error("OCA hybrid redistribution size exceeds INT_MAX");
    vector<OCAPathRecord> send_buffer(paths.size());
    vector<int> cursor=send_displacements;
    for (const auto& path:paths)
      send_buffer[cursor[destination_shard(path)]++]=path;
    vector<OCAPathRecord> receive_buffer(static_cast<size_t>(recv_total));
    MPI_Datatype record_type;
    MPI_Type_contiguous(sizeof(OCAPathRecord),MPI_BYTE,&record_type);
    MPI_Type_commit(&record_type);
    MPI_Alltoallv(send_buffer.data(),send_counts.data(),send_displacements.data(),
                  record_type,receive_buffer.data(),recv_counts.data(),
                  recv_displacements.data(),record_type,MPI_COMM_WORLD);
    MPI_Type_free(&record_type);
    paths.swap(receive_buffer);
  }

  // World ranks 0..shard_count-1 are the leaders of the cyclic groups.  After
  // the redistribution only those ranks own paths; broadcast each leader's
  // shard to its pole workers without ever constructing the global topology.
  MPI_Comm shard_comm=MPI_COMM_NULL;
  const int shard_id=mpi.rank%shard_count;
  MPI_Comm_split(MPI_COMM_WORLD,shard_id,mpi.rank,&shard_comm);
  int shard_rank=0;
  MPI_Comm_rank(shard_comm,&shard_rank);
  unsigned long long shard_paths=shard_rank==0 ? paths.size() : 0;
  MPI_Bcast(&shard_paths,1,MPI_UNSIGNED_LONG_LONG,0,shard_comm);
  if (shard_paths>static_cast<unsigned long long>(INT_MAX)) {
    MPI_Comm_free(&shard_comm);
    throw std::overflow_error(
      "CATASTROPHIC ERROR: one OCA topology shard exceeds INT_MAX paths");
  }
  if (shard_rank!=0) paths.resize(static_cast<size_t>(shard_paths));
  MPI_Datatype record_type;
  MPI_Type_contiguous(sizeof(OCAPathRecord),MPI_BYTE,&record_type);
  MPI_Type_commit(&record_type);
  MPI_Bcast(paths.data(),static_cast<int>(shard_paths),record_type,0,shard_comm);
  MPI_Type_free(&record_type);
  MPI_Comm_free(&shard_comm);
}
#endif

inline bool LoadOCAPathCache(const std::string& filename,
                             std::uint64_t cix_hash,
                             bool cmpPhysG,
                             const Cix& cix,
                             const vector<int>& importance,
                             int cut_sum,
                             const std::string& topology_mode,
                             double max_replicated_megabytes,
                             MPIcls& mpi,
                             OCAEnumeration& result)
{
  if (filename.empty() || cix_hash==0) return false;
  OCAPathCacheHeader header;
  bool valid=false;
  if (mpi.rank==mpi.master) {
    std::ifstream in(filename,std::ios::binary);
    if (in) {
      in.seekg(0,std::ios::end);
      const std::uint64_t file_bytes=static_cast<std::uint64_t>(in.tellg());
      in.seekg(0,std::ios::beg);
      in.read(reinterpret_cast<char*>(&header),sizeof(header));
      const char expected[8]={'O','C','A','D','L','R','P','1'};
      valid=static_cast<bool>(in) &&
        std::memcmp(header.magic,expected,8)==0 && header.version==1 &&
        header.record_size==sizeof(OCAPathRecord) &&
        header.cix_hash==cix_hash &&
        header.importance_hash==OCAImportanceHash(importance) &&
        header.npseudo==cix.Npseudo_unique &&
        header.nbath==cix.Nbath_unique_unique &&
        header.cut_sum==cut_sum &&
        header.includes_physical==static_cast<int>(cmpPhysG) &&
        header.record_count<=static_cast<std::uint64_t>(INT_MAX) &&
        file_bytes==sizeof(OCAPathCacheHeader)+
                    header.record_count*sizeof(OCAPathRecord);
    }
  }
#ifdef _MPI
  int valid_int=valid ? 1 : 0;
  MPI_Bcast(&valid_int,1,MPI_INT,mpi.master,MPI_COMM_WORLD);
  valid=valid_int!=0;
  if (!valid) return false;
  MPI_Bcast(&header,sizeof(header),MPI_BYTE,mpi.master,MPI_COMM_WORLD);
#else
  if (!valid) return false;
#endif

  result.includes_physical=cmpPhysG;
  result.npseudo=cix.Npseudo_unique;
  result.nbath_unique=cix.Nbath_unique_unique;
  result.cut_sum=cut_sum;
  result.importance=importance;
  result.global_nonzero_path_count=header.record_count;
  constexpr double compiled_plan_factor=4.5;
  result.estimated_replicated_megabytes=compiled_plan_factor*
    header.record_count*sizeof(OCACompactPath)/(1024.0*1024.0);
  result.topology_shard_count=OCAChooseTopologyShardCount(
    topology_mode,result.estimated_replicated_megabytes,
    max_replicated_megabytes,mpi.size);
  result.paths_distributed=result.topology_shard_count>1;
  OCASetTopologyGroupMetadata(result,mpi);

  std::uint64_t first=0,count=header.record_count;
  if (result.paths_distributed) {
    // Read one unique flat interval per world rank first.  The inexpensive
    // hybrid redistribution below restores bath-pair locality and then
    // replicates only the selected shard inside each pole-worker group.
    const std::uint64_t base=header.record_count/mpi.size;
    const std::uint64_t extra=header.record_count%mpi.size;
    count=base+(static_cast<std::uint64_t>(mpi.rank)<extra ? 1 : 0);
    first=base*mpi.rank+std::min<std::uint64_t>(mpi.rank,extra);
  }
  result.paths.resize(static_cast<std::size_t>(count));
  const std::uint64_t data_offset=sizeof(OCAPathCacheHeader)+
    first*sizeof(OCAPathRecord);
#ifdef _MPI
  MPI_File file=MPI_FILE_NULL;
  if (MPI_File_open(MPI_COMM_WORLD,filename.c_str(),MPI_MODE_RDONLY,
                    MPI_INFO_NULL,&file)!=MPI_SUCCESS)
    throw std::runtime_error("failed to open OCA DLR path cache with MPI-IO");
  MPI_Datatype record_type;
  MPI_Type_contiguous(sizeof(OCAPathRecord),MPI_BYTE,&record_type);
  MPI_Type_commit(&record_type);
  MPI_File_read_at_all(file,static_cast<MPI_Offset>(data_offset),
                       result.paths.data(),static_cast<int>(count),record_type,
                       MPI_STATUS_IGNORE);
  MPI_Type_free(&record_type);
  MPI_File_close(&file);
#else
  std::ifstream in(filename,std::ios::binary);
  in.seekg(static_cast<std::streamoff>(data_offset));
  in.read(reinterpret_cast<char*>(result.paths.data()),
          static_cast<std::streamsize>(count*sizeof(OCAPathRecord)));
  if (!in) throw std::runtime_error("failed to read OCA DLR path cache");
#endif
#ifdef _MPI
  if (result.paths_distributed)
    DistributeOCAPathRecordsHybridMPI(
      result.paths,result.nbath_unique,result.topology_shard_count,
      false,mpi);
#endif
  if (mpi.rank==mpi.master)
    (*mpi.out) << "OCA disk topology cache hit: file=" << filename
               << " paths=" << header.record_count << " topology="
               << (result.paths_distributed ? "hybrid-distributed" : "replicated")
               << " topology_shards=" << result.topology_shard_count
               << " workers_per_shard_range=["
               << mpi.size/result.topology_shard_count << ","
               << (mpi.size+result.topology_shard_count-1)/
                    result.topology_shard_count << "]"
               << std::endl;
  return true;
}

inline void WriteOCAPathCache(const std::string& filename,
                              std::uint64_t cix_hash,
                              const OCAEnumeration& enumeration,
                              MPIcls& mpi)
{
  if (filename.empty() || cix_hash==0) return;
  OCAPathCacheHeader header;
  header.cix_hash=cix_hash;
  header.importance_hash=OCAImportanceHash(enumeration.importance);
  header.npseudo=enumeration.npseudo;
  header.nbath=enumeration.nbath_unique;
  header.cut_sum=enumeration.cut_sum;
  header.includes_physical=enumeration.includes_physical;

  // Replicated enumeration is written by the master only.  Distributed
  // enumeration contributes one disjoint record interval per rank.  MPI-IO
  // records the intervals in rank order; a later run may repartition this
  // flat stream over any different number of ranks.
  std::uint64_t local_count=enumeration.paths.size();
  if ((!enumeration.paths_distributed && mpi.rank!=mpi.master) ||
      (enumeration.paths_distributed &&
       enumeration.topology_group_rank!=0))
    local_count=0;
  std::uint64_t first=0,global_count=local_count;
#ifdef _MPI
  MPI_Exscan(&local_count,&first,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
  if (mpi.rank==mpi.master) first=0;
  MPI_Allreduce(&local_count,&global_count,1,MPI_UINT64_T,MPI_SUM,
                MPI_COMM_WORLD);
#endif
  header.record_count=global_count;
  const std::string temporary=filename+".tmp";
#ifdef _MPI
  MPI_File file=MPI_FILE_NULL;
  if (MPI_File_open(MPI_COMM_WORLD,temporary.c_str(),
                    MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,
                    &file)!=MPI_SUCCESS)
    throw std::runtime_error("failed to create OCA DLR path cache with MPI-IO");
  MPI_File_set_size(file,static_cast<MPI_Offset>(
    sizeof(header)+global_count*sizeof(OCAPathRecord)));
  if (mpi.rank==mpi.master)
    MPI_File_write_at(file,0,&header,sizeof(header),MPI_BYTE,MPI_STATUS_IGNORE);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Datatype record_type;
  MPI_Type_contiguous(sizeof(OCAPathRecord),MPI_BYTE,&record_type);
  MPI_Type_commit(&record_type);
  const MPI_Offset offset=sizeof(header)+first*sizeof(OCAPathRecord);
  MPI_File_write_at_all(file,offset,
                        const_cast<OCAPathRecord*>(enumeration.paths.data()),
                        static_cast<int>(local_count),record_type,
                        MPI_STATUS_IGNORE);
  MPI_Type_free(&record_type);
  MPI_File_close(&file);
  int rename_error=0;
  if (mpi.rank==mpi.master)
    rename_error=std::rename(temporary.c_str(),filename.c_str())!=0;
  MPI_Bcast(&rename_error,1,MPI_INT,mpi.master,MPI_COMM_WORLD);
  if (rename_error)
    throw std::runtime_error("failed to publish OCA DLR path cache");
  MPI_Barrier(MPI_COMM_WORLD);
#else
  std::ofstream out(temporary,std::ios::binary|std::ios::trunc);
  out.write(reinterpret_cast<const char*>(&header),sizeof(header));
  out.write(reinterpret_cast<const char*>(enumeration.paths.data()),
            enumeration.paths.size()*sizeof(OCAPathRecord));
  out.close();
  if (!out || std::rename(temporary.c_str(),filename.c_str())!=0)
    throw std::runtime_error("failed to publish OCA DLR path cache");
#endif
  if (mpi.rank==mpi.master)
    (*mpi.out) << "OCA disk topology cache written: file=" << filename
               << " paths=" << global_count << std::endl;
}

inline bool OCAPathWeightsAreZero(const std::array<double,5>& weights)
{
  return weights[0]==0.0 && weights[1]==0.0 && weights[2]==0.0 &&
    weights[3]==0.0 && weights[4]==0.0;
}

// Reconstruct every legacy convolution index from the globally reduced path
// records.  Local integer indices cannot be communicated between MPI ranks:
// each rank discovered a different subset and therefore assigned different
// insertion-order numbers.  Rebuilding after the global key reduction is
// deterministic, linear in the final path count, and also discards inner
// pairs belonging only to paths whose accumulated weights cancelled to zero.
inline void RebuildOCAInnerCatalog(OCAEnumeration& result)
{
  result.unique_inner_states.assign(4,{});
  std::array<unordered_map<TinyVector<int,4>,int,TV4Hash,TV4Equal>,4> index;
  for (auto& map:index) map.reserve(result.compact_paths.size());

  const auto add_inner=[&](int cut, TinyVector<int,4> state,
                           bool needed_for_self_energy) {
    // The fourth component is a usage flag rather than part of the identity.
    // Hash only the three physical indices, exactly as the original code did
    // by inserting the key with flag=1 and changing stored values to zero.
    state[3]=1;
    auto [it,inserted]=index[cut].try_emplace(
      state,static_cast<int>(result.unique_inner_states[cut].size()));
    if (inserted) {
      state[3]=needed_for_self_energy ? 1 : 0;
      result.unique_inner_states[cut].push_back(state);
    } else if (needed_for_self_energy) {
      result.unique_inner_states[cut][it->second][3]=1;
    }
    return it->second;
  };

  for (auto& path:result.compact_paths) {
    path.inner_index={-1,-1,-1,-1};
    const TinyVector<int,6> st=unpack6(path.key);
    const int s1=st[0],s0=st[1];
    const int p3=st[2],p2=st[3],p1=st[4],p0=st[5];
    const bool self=path.weights[0]!=0.0;
    if (self) {
      path.inner_index[0]=add_inner(0,TinyVector<int,4>(p0,p1,s0,1),true);
      path.inner_index[1]=add_inner(1,TinyVector<int,4>(p1,p2,s1,1),true);
    }
    // Physical G/F cuts need only the p2-p3 and p3-p0 inner pairs.  They are
    // also required by every ordinary self-energy path.
    path.inner_index[2]=add_inner(2,TinyVector<int,4>(p2,p3,s0,1),self);
    path.inner_index[3]=add_inner(3,TinyVector<int,4>(p3,p0,s1,1),self);
  }
}

#ifdef _MPI
struct OCAWirePath {
  OCAKey key{};
  std::array<double,5> weights{};
};
static_assert(std::is_trivially_copyable_v<OCAWirePath>);

// Gather rank-local path records without a root bottleneck.  The scheduler
// assigns an entire reduced bath pair (s0,s1) to one rank.  Since (s0,s1) is
// part of every packed key, two ranks can never generate the same key: all
// duplicate contributions have already met in the producing rank's local
// hash table.  We can consequently skip the former hash-owner Alltoallv and
// directly Allgatherv the already unique records needed by every evaluator.
inline vector<OCACompactPath> GatherOCAPathsMPI(
  const vector<OCACompactPath>& local_paths,MPIcls& mpi)
{
  if (mpi.size==1) return local_paths;

  if (local_paths.size()>static_cast<size_t>(INT_MAX))
    throw std::overflow_error(
      "CATASTROPHIC ERROR: one MPI rank generated more than INT_MAX OCA paths.");
  vector<OCAWirePath> local_wire;
  local_wire.reserve(local_paths.size());
  for (const auto& path:local_paths) {
    local_wire.push_back({path.key,path.weights});
  }

  const int local_count=static_cast<int>(local_wire.size());
  vector<int> owned_counts(mpi.size),owned_displs(mpi.size);
  MPI_Allgather(&local_count,1,MPI_INT,
                owned_counts.data(),1,MPI_INT,MPI_COMM_WORLD);
  long long global_total=0;
  for (int rank=0;rank<mpi.size;++rank) {
    if (global_total>INT_MAX)
      throw std::overflow_error(
        "CATASTROPHIC ERROR: OCA MPI Allgatherv displacement exceeds INT_MAX. "
        "Chunked communication is required.");
    owned_displs[rank]=static_cast<int>(global_total);
    global_total+=owned_counts[rank];
  }
  if (global_total>INT_MAX)
    throw std::overflow_error(
      "CATASTROPHIC ERROR: global OCA path count exceeds INT_MAX. "
      "Chunked communication is required.");
  vector<OCAWirePath> global_wire(static_cast<size_t>(global_total));
  MPI_Datatype wire_type;
  MPI_Type_contiguous(static_cast<int>(sizeof(OCAWirePath)),MPI_BYTE,&wire_type);
  MPI_Type_commit(&wire_type);
  MPI_Allgatherv(local_wire.data(),local_count,wire_type,
                 global_wire.data(),owned_counts.data(),owned_displs.data(),wire_type,
                 MPI_COMM_WORLD);
  MPI_Type_free(&wire_type);

  vector<OCACompactPath> global_paths;
  global_paths.reserve(global_wire.size());
  for (const auto& wire:global_wire) {
    OCACompactPath path;
    path.key=wire.key;
    path.weights=wire.weights;
    path.inner_index={-1,-1,-1,-1};
    global_paths.push_back(path);
  }
  std::sort(global_paths.begin(),global_paths.end(),
            [](const OCACompactPath& a,const OCACompactPath& b) {
              return a.key<b.key;
            });
  // A duplicate here would prove that the scheduling invariant was broken:
  // equal packed keys necessarily have equal (s0,s1), and an entire pair must
  // be owned by exactly one rank.  Never silently double-count such a bug.
  for (size_t i=1;i<global_paths.size();++i)
    if (global_paths[i-1].key==global_paths[i].key)
      throw std::logic_error(
        "CATASTROPHIC ERROR: duplicate OCA key crossed MPI ranks despite "
        "reduced-bath-pair scheduling.");
  return global_paths;
}

// Rebalance a non-replicated path stream over every MPI rank after the
// bath-pair-owned enumeration.  Enumeration has at most Nbath^2 useful
// owners, but integration may profit from more ranks.  Hash partitioning the
// already compact records is a single O(number_of_paths) Alltoallv, keeps
// exactly one copy globally, and makes all ranks useful during plan building
// and every subsequent AAA-pole contraction.
inline void RedistributeOCAPathRecordsMPI(vector<OCAPathRecord>& paths,
                                          MPIcls& mpi)
{
  if (mpi.size<=1) return;
  vector<int> send_counts(mpi.size,0),recv_counts(mpi.size,0);
  const auto destination=[&](const OCAPathRecord& path) {
    std::uint64_t hash=1469598103934665603ULL;
    hash=OCAFNV1a(path.pseudo.data(),path.pseudo.size()*sizeof(int),hash);
    hash=OCAFNV1a(&path.s0,sizeof(path.s0),hash);
    hash=OCAFNV1a(&path.s1,sizeof(path.s1),hash);
    return static_cast<int>(hash%mpi.size);
  };
  for (const auto& path:paths) {
    const int target=destination(path);
    if (send_counts[target]==INT_MAX)
      throw std::overflow_error("OCA path redistribution count exceeds INT_MAX");
    ++send_counts[target];
  }
  MPI_Alltoall(send_counts.data(),1,MPI_INT,
               recv_counts.data(),1,MPI_INT,MPI_COMM_WORLD);
  vector<int> send_displacements(mpi.size),recv_displacements(mpi.size);
  long long send_total=0,recv_total=0;
  for (int rank=0;rank<mpi.size;++rank) {
    if (send_total>INT_MAX || recv_total>INT_MAX)
      throw std::overflow_error(
        "OCA path redistribution displacement exceeds INT_MAX");
    send_displacements[rank]=static_cast<int>(send_total);
    recv_displacements[rank]=static_cast<int>(recv_total);
    send_total+=send_counts[rank];
    recv_total+=recv_counts[rank];
  }
  if (send_total>INT_MAX || recv_total>INT_MAX)
    throw std::overflow_error("OCA path redistribution size exceeds INT_MAX");
  vector<OCAPathRecord> send_buffer(paths.size());
  vector<int> cursor=send_displacements;
  for (const auto& path:paths)
    send_buffer[cursor[destination(path)]++]=path;
  vector<OCAPathRecord> receive_buffer(static_cast<size_t>(recv_total));
  MPI_Datatype record_type;
  MPI_Type_contiguous(sizeof(OCAPathRecord),MPI_BYTE,&record_type);
  MPI_Type_commit(&record_type);
  MPI_Alltoallv(send_buffer.data(),send_counts.data(),send_displacements.data(),
                record_type,receive_buffer.data(),recv_counts.data(),
                recv_displacements.data(),record_type,MPI_COMM_WORLD);
  MPI_Type_free(&record_type);
  paths.swap(receive_buffer);

  unsigned long long local=paths.size(),minimum=local,maximum=local,total=local;
  MPI_Allreduce(&local,&minimum,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&local,&maximum,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&local,&total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  if (mpi.rank==mpi.master)
    (*mpi.out) << "OCA distributed-plan rebalance: rank_path_range=["
               << minimum << "," << maximum << "] total=" << total
               << std::endl;
}
#endif

// Enumerate the OCA topology once, before selecting a numerical integrator.
// Importance filtering, CIX selection rules, operator signs, and Bulla matrix
// replacements therefore remain exactly common to legacy and DLR paths.
// A bounded consumer receives partial key sums, NOT complete diagrams. It must
// combine equal full keys (including both bath indices) before evaluating them.
// Only DLR production uses this hook; the original two-pass path remains intact.
// The map and weight vector have fixed capacity; no global count pass is needed.
struct OCABoundedEnumeration {
  size_t max_paths=100000;
  bool distributed=true;
  std::function<void(const std::vector<OCACompactPath>&)> consume;
};

inline OCAEnumeration EnumerateOCAPaths(bool cmpPhysG,
				       Cix& cix,
				       const vector<int>& G_importance,
				       int OCA_cut_sum,
				       bool build_legacy_catalog,
				       bool export_path_records,
				       const std::string& mpi_topology_mode,
				       double max_replicated_megabytes,
                                       MPIcls& mpi,
                                       const OCABoundedEnumeration* bounded=nullptr)
{
  if (bounded && (build_legacy_catalog || !bounded->consume ||
      bounded->max_paths==0 || bounded->max_paths>static_cast<size_t>(INT_MAX)))
    throw std::invalid_argument("invalid bounded OCA enumeration options");
  if (mpi_topology_mode!="auto" && mpi_topology_mode!="replicated" &&
      mpi_topology_mode!="distributed")
    throw std::invalid_argument(
      "oca_dlr_mpi_topology must be \"auto\", \"replicated\", or "
      "\"distributed\"");
  ////////////////////////////
  // Diagram convention:
  //                    nds1
  //               Fi[2*is1,p1] 
  //                  --->---
  //   states[1]=p1 /  s1^    \  states[2]=p2
  //               |s0   |   s0|  
  //Fi[2*is0+1,p0] |---->------| Fi[2*is0,p2] ; nds0
  //               |     |     |
  //   states[0]=p0 \  s1|    /  states[3]=p3
  //                  ---<---
  //               Fi[2*is1+1,p3]
  //  states[] for four diagrams, where 0=|0>, 1=|up>, 2=|dn>, 3=|2>
  //  Ds[] 0=up, 1=dn
  const int Nd=4;
  Diag2Matrix diag2matrix(6, cix);
  Array<int,2> diag(2,2);     // OCA diagram
  diag = 0,2,1,3;             // OCA diagram
  Array<int,1> dirspin(2);    // dirspin=2*iorb1,2*iorb0

  using EnumerationClock=std::chrono::steady_clock;
  const auto scheduling_start=EnumerationClock::now();

  // Discover only the inexpensive outer work units before entering either
  // combinatorial pass.  Both passes replay this exact ordering and consult
  // task_owner, so a rank works on the same CIX blocks in both passes.
  //
  // The indivisible scheduling unit is a *reduced bath pair* (s0,s1), not an
  // individual block.  Both bath indices are fields of pack6(), so all terms
  // which can possibly reduce to the same packed path key are thereby kept on
  // one rank.  The local unordered_map performs the complete reduction and
  // the final MPI stage needs only an Allgatherv, with no duplicate-key
  // Alltoallv.  For an f shell, 14 baths provide as many as 14*14=196 such
  // independently schedulable groups.
  struct EnumerationTaskCost {
    size_t sequence=0;
    long double cost=0.0;
  };
  struct EnumerationPairCost {
    int s0=0;
    int s1=0;
    long double cost=0.0;
    vector<size_t> tasks;
  };
  vector<EnumerationTaskCost> task_costs;
  vector<EnumerationPairCost> pair_costs;
  map<pair<int,int>,size_t> pair_lookup;
  for (int iorb0=0;iorb0<cix.Nbath;++iorb0) {
    const int s0=abs(cix.bath_sigind0(iorb0));
    for (int iorb1=0;iorb1<cix.Nbath;++iorb1) {
      const int s1=abs(cix.bath_sigind0(iorb1));
      dirspin=2*iorb1,2*iorb0;
      auto& allowed=diag2matrix.Cmp_matrix(diag,dirspin,cix);
      for (const auto& [ib,alw]:allowed) {
        const long double n0=alw[0].extent(0),n1=alw[1].extent(0);
        const long double n2=alw[2].extent(0),n3=alw[3].extent(0);
        if (n0==0 || n1==0 || n2==0 || n3==0) continue;
        long double cost=n0*n1*n2*n3;
        if (cmpPhysG) {
          TinyVector<int,Nd> ibs;
          ibs[0]=cix.Fbc(diag2matrix.op(0),ib);
          ibs[1]=cix.Fbc(diag2matrix.op(1),ibs[0]);
          ibs[2]=cix.Fbc(diag2matrix.op(2),ibs[1]);
          ibs[3]=ib;
          const NjmResult nx1=cix.Njm_ind(ibs[1],2*iorb1);
          const NjmResult nx2=cix.Njm_ind(ibs[2],2*iorb0);
          if (nx2.njm_c==2)
            cost+=n0*n1*n2*cix.Block[ibs[2]].size()*cix.Block[ibs[3]].size();
          if (nx1.njm_c==2)
            cost+=n0*n1*n3*cix.Block[ibs[1]].size()*cix.Block[ibs[2]].size();
        }
        cost=std::max(cost,1.0L);
        const size_t task_sequence=task_costs.size();
        task_costs.push_back({task_sequence,cost});
        const auto [pair_it,pair_inserted]=pair_lookup.try_emplace(
          make_pair(s0,s1),pair_costs.size());
        if (pair_inserted) pair_costs.push_back({s0,s1,0.0L,{}});
        auto& group=pair_costs[pair_it->second];
        group.cost+=cost;
        group.tasks.push_back(task_sequence);
      }
    }
  }
  vector<int> task_owner(task_costs.size(),0);
  vector<long double> rank_load(std::max(mpi.size,1),0.0L);
  vector<size_t> sorted_pairs(pair_costs.size());
  std::iota(sorted_pairs.begin(),sorted_pairs.end(),0);
  std::stable_sort(sorted_pairs.begin(),sorted_pairs.end(),
                   [&](size_t a,size_t b) {
                     return pair_costs[a].cost>pair_costs[b].cost;
                   });
  for (const size_t pair_index:sorted_pairs) {
    const auto& group=pair_costs[pair_index];
    const int owner=static_cast<int>(std::min_element(
      rank_load.begin(),rank_load.end())-rank_load.begin());
    for (const size_t task_sequence:group.tasks)
      task_owner[task_sequence]=owner;
    rank_load[owner]+=group.cost;
  }
  const long double total_load=
    std::accumulate(rank_load.begin(),rank_load.end(),0.0L);
  const auto [min_load_it,max_load_it]=
    std::minmax_element(rank_load.begin(),rank_load.end());
  const auto scheduling_end=EnumerationClock::now();
  if (mpi.rank==mpi.master) {
    const long double mean_load=total_load/std::max(mpi.size,1);
    const int active_ranks=static_cast<int>(std::count_if(
      rank_load.begin(),rank_load.end(),[](long double load){return load>0;}));
    (*mpi.out) << "OCA enumeration scheduler: tasks=" << task_costs.size()
               << " reduced_bath_pairs=" << pair_costs.size()
               << " ranks=" << mpi.size
               << " active_ranks=" << active_ranks
               << " estimated_load_range=[" << static_cast<double>(*min_load_it)
               << "," << static_cast<double>(*max_load_it) << "]"
               << " max_over_mean="
               << static_cast<double>(mean_load>0 ? *max_load_it/mean_load : 0)
               << " scheduling_seconds="
               << std::chrono::duration<double>(scheduling_end-scheduling_start).count()
               << endl;
  }
  const auto pass1_start=EnumerationClock::now();

  // Data which store the components of the diagram, so that integration looks like a problem of a single orbital.
  int states_unique_size=0;
  if (!bounded) {// This part is carbon-copy of the algorithm below, except that we only calculate the size of the hash-table states_unique.
    // You should not try to understand that, because we comment only the complete algorithm below. This is just stripped/off equivalent,
    // which does not compute anything except total size of the storage.
    unordered_map<OCAKey,int,OCAKeyHash> states_unique;
    states_unique.max_load_factor(0.7f);
    const long double reserve_share=total_load>0 ?
      rank_load[mpi.rank]/total_load : 1.0L/std::max(mpi.size,1);
    const long double global_reserve=static_cast<long double>(cix.Natst)*
      cix.Natst*cix.Nbath*cix.Nbath;
    states_unique.reserve(static_cast<size_t>(global_reserve*reserve_share)+1);
    // for each unique key in states_unique map (which uniquely specifies the diagram contribution) the value of the map gives index to array mm_nds, index0,.., index3
    int mm_nds_size=0; // contains the matrix elements for Phi, and Nds matrix elements for vertex F(omega) calculation. There are two Nds values, because of two hybridizations.
    size_t task_cursor=0;
    // Below we find all contributions to the OCA Feynman diagram.
    for (int iorb0=0; iorb0<cix.Nbath; ++iorb0){ // over all baths for first hybridiation
      int s0 = abs(cix.bath_sigind0(iorb0));
      for (int iorb1=0; iorb1<cix.Nbath; ++iorb1){ // bath for the second hybridization
	int s1 = abs(cix.bath_sigind0(iorb1));
	dirspin = 2*iorb1, 2*iorb0;                // directions are spins choosen to go in positive direction
	auto& allowed  = diag2matrix.Cmp_matrix(diag,dirspin,cix); // computes all allowed Green's functions that survice Pauli principle.
	for (const auto& [ib, alw] : allowed){     // going over all allowed contributions from all blocks.
	  if (alw[0].extent(0)>0 and alw[1].extent(0)>0 and alw[2].extent(0)>0 and alw[3].extent(0)>0){
	    const size_t current_task=task_cursor++;
	    if ((!bounded || bounded->distributed) && task_owner[current_task]!=mpi.rank) continue;
	    TinyVector<int,Nd> ibs;  // the four blocks for pseudo-particles (see figure below). ib is the index to the last block.
	    ibs[0] = cix.Fbc(diag2matrix.op(0),ib);      // op(0) is computed before, and can be applied to ib to get ibs[0] (see figure)
	    ibs[1] = cix.Fbc(diag2matrix.op(1),ibs[0]);  // op(1) is applied to ibs[0] to get ibs[1], etc...
	    ibs[2] = cix.Fbc(diag2matrix.op(2),ibs[1]);
	    ibs[3] = ib;
	    // This part is for the vertex F(omega) calculation, where we need Njm matrices for both hybridizations.
	    NjmResult nx1 = cix.Njm_ind(ibs[1],2*iorb1);
	    Array<double,2> njm1;           // matrix when njm is not proportional to F
	    double njmr1 = nx1.njm_r;       // proportionality constant for njm
	    if (nx1.njm_c==2){               // njm is not proportional
	      if (not nx1.transposed){
		njm1.resize(nx1.njm->extent(0),nx1.njm->extent(1));
		njm1 = *nx1.njm;
	      }else{
		njm1.resize(nx1.njm->extent(1),nx1.njm->extent(0));
		for (int n1=0; n1<njm1.extent(0); ++n1)
		  for (int n2=0; n2<njm1.extent(1); ++n2)
		    njm1(n1,n2) = (*nx1.njm)(n2,n1);
	      }
	    }
	    NjmResult nx2 = cix.Njm_ind(ibs[2],2*iorb0);
	    Array<double,2> njm2;           // matrix when njm is not proportional to F
	    double njmr2 = nx2.njm_r;       // proportionality constant for njm
	    if (nx2.njm_c==2){              // njm is not proportional
	      if (not nx2.transposed){
		njm2.resize(nx2.njm->extent(0),nx2.njm->extent(1));
		njm2 = *nx2.njm;
	      }else{
		njm2.resize(nx2.njm->extent(1),nx2.njm->extent(0));
		for (int n1=0; n1<njm2.extent(0); ++n1)
		  for (int n2=0; n2<njm2.extent(1); ++n2)
		    njm2(n1,n2) = (*nx2.njm)(n2,n1);
	      }
	    }
	    TinyVector<int,Nd> js;
	    for (js[0]=0; js[0]<alw[0].extent(0); ++js[0]){
	      int a0 = cix.inBlock(abs(alw[0](js[0],1))); // alw[0] index is given in terms of atomic state index, but a0 is index within the block. cix.inBlock translates this.
	      int b3 = cix.inBlock(abs(alw[0](js[0],0))); // cix.inBlock translates again to index inside the block ib
	      for (js[1]=0; js[1]<alw[1].extent(0); ++js[1]){
		int b0 = cix.inBlock(abs(alw[1](js[1],0)));
		int _ii0_ = cix.Sgm_symmetry[ibs[0]](a0,b0); // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
		if (_ii0_==0) continue;
		int ii0 = abs(_ii0_)-1;
		int imp_ii0 = G_importance[ii0];
		int ii_sign0 = sign(_ii0_);                 // Green's function value for off-diagonal propagators could be negative/opposite, hence sign
		int a1 = cix.inBlock(abs(alw[1](js[1],1)));
		for (js[2]=0; js[2]<alw[2].extent(0); ++js[2]){
		  int b1 = cix.inBlock(abs(alw[2](js[2],0)));
		  int _ii1_ = cix.Sgm_symmetry[ibs[1]](a1,b1); // Green's function idex G(a1,b1)==G[|_ii1_|-1]
		  if (_ii1_==0) continue;
		  int ii1 = abs(_ii1_)-1;
		  int imp_ii1 = G_importance[ii1];
		  int ii_sign1 = sign(_ii1_) * ii_sign0;
		  int a2 = cix.inBlock(abs(alw[2](js[2],1)));
		  for (js[3]=0; js[3]<alw[3].extent(0); ++js[3]){
		    int b2 = cix.inBlock(abs(alw[3](js[3],0)));
		    int _ii2_ = cix.Sgm_symmetry[ibs[2]](a2,b2); // Green's function index G(a2,b2)==G[|_ii2_|-1]
		    if (_ii2_==0) continue;
		    int ii2 = abs(_ii2_)-1;
		    int imp_ii2 = G_importance[ii2];
		    int ii_sign2 = sign(_ii2_) * ii_sign1;
		    int a3 = cix.inBlock(abs(alw[3](js[3],1)));
		    int _ii3_ = cix.Sgm_symmetry[ibs[3]](a3,b3); // Green's function index G(a3,b3)==G[|_ii3_|-1]
		    if (_ii3_==0) continue;
		    int ii3 = abs(_ii3_)-1;
		    int importance = G_importance[ii3] + imp_ii2 + imp_ii1 + imp_ii0;
		    //cout<<"importance="<<importance<<" i3="<<G_importance[ii3] <<" i2="<<imp_ii2 <<" i1="<<imp_ii1 << " i0="<<imp_ii0 << endl;
		    if (importance > OCA_cut_sum ) continue; // this cuts-off PP-states that are less important
		    int ii_sign3 = sign(_ii3_) * ii_sign2;
		    TinyVector<int,4> states(ii0, ii1, ii2, ii3); // atomic states indices between [0,...cix.Npseudo_unique-1] array.
		    //
		    OCAKey st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // we pack the two hybridizations (s1,s0) and atomic states into unique index-vector.
		    auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // try to insert st-vector as a new key and mm_nds.size() as the value.
		    // If the key is not yet present, inserted=true and "it" points to its current place. If key is present, "it" points to its existing place.
		    int iind = it->second;                    // The key st has value iind, either just inserted, or found from before
		    if (inserted) mm_nds_size += 1;
		  }
		  if (cmpPhysG and nx2.njm_c==2){
		    for (int _b2_=0; _b2_<cix.Block[ibs[2]].size(); ++_b2_){
		      int __ii2__ = cix.Sgm_symmetry[ibs[2]](a2,_b2_);
		      if (__ii2__==0) continue;
		      int ii2 = abs(__ii2__)-1;
		      int imp_ii2 = G_importance[ii2];
		      for (int _a3_=0; _a3_<cix.Block[ibs[3]].size(); ++_a3_){
			if (njm2(_b2_,_a3_)==0) continue;
			int __ii3__ = cix.Sgm_symmetry[ibs[3]](_a3_,b3);
			if (__ii3__==0) continue;
			int ii3 = abs(__ii3__)-1;
			int importance = G_importance[ii3] + imp_ii2 + imp_ii1 + imp_ii0;
			if (importance > OCA_cut_sum ) continue;
			TinyVector<int,4> states(ii0, ii1, ii2, ii3); // atomic states indices between [0,...cix.Npseudo_unique-1] array.
			OCAKey st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // all atributes of diagram contribution.
			auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // Did we see that contribution before?
			int _iind_ = it->second; // in mm_nds and indexi, it is stored in _iind_
			if (inserted) mm_nds_size += 1;
		      }
		    }
		  }
		}
		if (cmpPhysG and nx1.njm_c==2){
		  for (js[3]=0; js[3]<alw[3].extent(0); ++js[3]){
		    int b2 = cix.inBlock(abs(alw[3](js[3],0)));
		    int a3 = cix.inBlock(abs(alw[3](js[3],1)));
		    int _ii3_ = cix.Sgm_symmetry[ibs[3]](a3,b3); // Green's function index G(a3,b3)==G[|_ii3_|-1]
		    if (_ii3_==0) continue;
		    int ii3 = abs(_ii3_)-1;
		    int imp_ii3 = G_importance[ii3];
		    for (int _b1_=0; _b1_<cix.Block[ibs[1]].size(); ++_b1_){ // now going over all other possibilities
		      int __ii1__ = cix.Sgm_symmetry[ibs[1]](a1,_b1_); // a1 is the same, but b1 can be anything within the block
		      if (__ii1__==0) continue;
		      int ii1 = abs(__ii1__)-1;
		      int imp_ii1 = G_importance[ii1];
		      for (int _a2_=0; _a2_<cix.Block[ibs[2]].size(); ++_a2_){
			if (njm1(_b1_,_a2_)==0) continue;
			int __ii2__ = cix.Sgm_symmetry[ibs[2]](_a2_,b2);
			if (__ii2__==0) continue;
			int ii2 = abs(__ii2__)-1;
			int importance = G_importance[ii2] + imp_ii1 + imp_ii3 + imp_ii0;
			if (importance > OCA_cut_sum ) continue;
			TinyVector<int,4> states(ii0, ii1, ii2, ii3); // atomic states indices between [0,...cix.Npseudo_unique-1] array.
			OCAKey st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // all atributes of diagram contribution.
			auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // Did we see that contribution before?
			int _iind_ = it->second; // in mm_nds and indexi, it is stored in _iind_
			if (inserted) mm_nds_size += 1;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if (task_cursor!=task_owner.size())
      throw std::logic_error("OCA pass 1 task ordering disagrees with scheduler discovery");
    // Finally we know the size of the hash-table. Below we will go through the same algorithm once more and actually calculate matrix elements.
    states_unique_size = states_unique.size();
  }
  const auto pass1_end=EnumerationClock::now();
  if (!bounded)
    (*mpi.out)<<"local_states_unique.size="<<states_unique_size
              <<" rank="<<mpi.rank<<endl;
  else
    states_unique_size=static_cast<int>(bounded->max_paths);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  const auto pass2_start=EnumerationClock::now();
  
  //map<TinyVector<int,6>,int,Class1Compare> states_unique; // every unique term is stored in this map with indices {s1,s0,states[3],states[2],states[1], states[0]}
  unordered_map<OCAKey,int,OCAKeyHash> states_unique;
  states_unique.reserve(states_unique_size);
  //states_unique.max_load_factor(0.7f);
  //states_unique.reserve(cix.Natst*cix.Natst*cix.Nbath*cix.Nbath);
  // for each unique key in states_unique map (which uniquely specifies the diagram contribution) the value of the map gives index to array mm_nds, index0,.., index3
  int mm_nds_size=0;
  vector<TinyVector<double,5>> mm_nds(states_unique_size); // contains the matrix elements for Phi, and Nds matrix elements for vertex F(omega) calculation. There are two Nds values, because of two hybridizations.
  // Flush only between complete insert/update operations. A later occurrence
  // of the same key can appear in another flush; the disk reducer reunites it.
  // Do not evaluate a flush directly: cancellation must precede plan building.
  const auto flush_bounded=[&]() {
    if (!bounded || states_unique.empty()) return;
    vector<OCACompactPath> records;
    records.reserve(states_unique.size());
    for (const auto& [key,index]:states_unique) {
      OCACompactPath record;
      record.key=key;
      for (int j=0;j<5;++j) record.weights[j]=mm_nds[index][j];
      if (!OCAPathWeightsAreZero(record.weights)) records.push_back(record);
    }
    bounded->consume(records);
    states_unique.clear(); // retain fixed bucket capacity for the next chunk
    mm_nds_size=0;
  };
  const auto prepare_bounded_insert=[&](const OCAKey& key) {
    if (bounded && states_unique.size()>=bounded->max_paths &&
        states_unique.find(key)==states_unique.end())
      flush_bounded();
  };
  // mm_nds[0] = mm    for diagram
  // mm_nds[1] = nds1  first hybridization Nds matrix element
  // mm_nds[2] = nds2  second hybridization Nds matrix element
  // mm_nds[3] = nds3  first hybridization Nds but for destroying Cooper pair or down in SC state
  // mm_nds[4] = nds4  second hybridization Nds but for destroying Cooper pair or down in SC state
  //
  // For efficiency, we store parts of the above unique contribution, which is needed in the inner loop of integration. Namely, many terms in states_unique share
  // some common parts, which appear in the inner integration loops. These can be enumerated and pre-computed. The table for such pre-computation are given by
  // unique_states_inside0, unique_states_inside1, unique_states_inside2, unique_states_inside3,
  // and to keep track of which contributions are unique, we also keep maps index_of_inside0, index_of_inside1, index_of_inside2, index_of_inside3.
  // In addition, index0, index1, index2 and index3 connect this partial contribution from the inside loop with the index of the entire contribution in states_unique.
  // These tables and their per-path indices are construction data for the
  // legacy spline integrator only.  A pure DLR run compiles OCAPathRecords
  // into its own sparse plans and must not pay their potentially enormous
  // allocation and hash-insertion cost.
  vector<vector<TinyVector<int,4>>> unique_states_inside;
  if (build_legacy_catalog) unique_states_inside.resize(4);
  //vector<int> index0(states_unique_size,-1), index1(states_unique_size,-1), index2(states_unique_size,-1), index3(states_unique_size,-1);
  vector<array<int,4>> indexi;
  if (build_legacy_catalog)
    indexi.assign(states_unique_size,array<int,4>{-1,-1,-1,-1});
  {// index_of_inside can be temporary
    unordered_map<TinyVector<int,4>,int,TV4Hash,TV4Equal> index_of_inside0, index_of_inside1, index_of_inside2, index_of_inside3; // value → index
    size_t task_cursor=0;
    // Below we find all contributions to the OCA Feynman diagram.
    for (int iorb0=0; iorb0<cix.Nbath; ++iorb0){ // over all baths for first hybridiation
      int s0 = abs(cix.bath_sigind0(iorb0));
      double s0_sgn = cix.bath_sigind0(iorb0) < 0 ? -1 : 1;
      for (int iorb1=0; iorb1<cix.Nbath; ++iorb1){ // bath for the second hybridization
	int s1 = abs(cix.bath_sigind0(iorb1));
	double s1_sgn = cix.bath_sigind0(iorb1) < 0 ? -1 : 1;
	dirspin = 2*iorb1, 2*iorb0;                // directions are spins choosen to go in positive direction
	bool destroys_Cooper_pair1 = cix.Fi_indx(2*iorb1)%2==1 and cix.Fi_indx(2*iorb1+1)%2==1;
	bool normal_down1 = cix.SC_bath(iorb1) and ((cix.Fi_indx(2*iorb1)+cix.Fi_indx(2*iorb1+1))%2==1 and cix.bath_spin(iorb1)<0);
	bool destroys_Cooper_pair_or_normal_down1 = destroys_Cooper_pair1 or normal_down1;
	bool destroys_Cooper_pair2 = cix.Fi_indx(2*iorb0)%2==1 and cix.Fi_indx(2*iorb0+1)%2==1;
	bool normal_down2 = cix.SC_bath(iorb0) and ((cix.Fi_indx(2*iorb0)+cix.Fi_indx(2*iorb0+1))%2==1 and cix.bath_spin(iorb0)<0);
	bool destroys_Cooper_pair_or_normal_down2 = destroys_Cooper_pair2 or normal_down2;
	auto& allowed  = diag2matrix.Cmp_matrix(diag,dirspin,cix); // computes all allowed Green's functions that survice Pauli principle.
	for (const auto& [ib, alw] : allowed){     // going over all allowed contributions from all blocks.
	  if (alw[0].extent(0)>0 and alw[1].extent(0)>0 and alw[2].extent(0)>0 and alw[3].extent(0)>0){
	    const size_t current_task=task_cursor++;
	    if ((!bounded || bounded->distributed) && task_owner[current_task]!=mpi.rank) continue;
	    TinyVector<int,Nd> ibs;  // the four blocks for pseudo-particles (see figure below). ib is the index to the last block.
	    ibs[0] = cix.Fbc(diag2matrix.op(0),ib);      // op(0) is computed before, and can be applied to ib to get ibs[0] (see figure)
	    ibs[1] = cix.Fbc(diag2matrix.op(1),ibs[0]);  // op(1) is applied to ibs[0] to get ibs[1], etc...
	    ibs[2] = cix.Fbc(diag2matrix.op(2),ibs[1]);
	    ibs[3] = ib;
	    // This part is for the vertex F(omega) calculation, where we need Njm matrices for both hybridizations.
	    //int njmc1 = cix.Njm_c(ibs[1],iorb1); // index to identify whether njm is proportional to F or not.
	    NjmResult nx1 = cix.Njm_ind(ibs[1],2*iorb1);
	    Array<double,2> njm1;           // matrix when njm is not proportional to F
	    double njmr1 = nx1.njm_r;       // proportionality constant for njm
	    if (nx1.njm_c==2){               // njm is not proportional
	      //njm1.reference(cix.Njm[ibs[1]][iorb1]);
	      if (not nx1.transposed){
		njm1.resize(nx1.njm->extent(0),nx1.njm->extent(1));
		njm1 = *nx1.njm;
	      }else{
		njm1.resize(nx1.njm->extent(1),nx1.njm->extent(0));
		for (int n1=0; n1<njm1.extent(0); ++n1)
		  for (int n2=0; n2<njm1.extent(1); ++n2)
		    njm1(n1,n2) = (*nx1.njm)(n2,n1);
	      }
	    }
	    //int njmc2 = cix.Njm_c(ibs[2],iorb0); // index to identify whether njm is proportional to F or not.
	    NjmResult nx2 = cix.Njm_ind(ibs[2],2*iorb0);
	    Array<double,2> njm2;           // matrix when njm is not proportional to F
	    double njmr2 = nx2.njm_r;       // proportionality constant for njm
	    if (nx2.njm_c==2){              // njm is not proportional
	      //njm2.reference(cix.Njm[ibs[2]][iorb0]);
	      if (not nx2.transposed){
		njm2.resize(nx2.njm->extent(0),nx2.njm->extent(1));
		njm2 = *nx2.njm;
	      }else{
		njm2.resize(nx2.njm->extent(1),nx2.njm->extent(0));
		for (int n1=0; n1<njm2.extent(0); ++n1)
		  for (int n2=0; n2<njm2.extent(1); ++n2)
		    njm2(n1,n2) = (*nx2.njm)(n2,n1);
	      }
	    }
	    //   This should make sense for below implementation of loops over j0,j1,j2,j3.
	    //
	    //             b3=alw[0](j0,0)   a0=alw[0](j0,1)
	    //                           op(0)
	    //                         --->---                
	    //states[3],ibs[3],ii3   /   s1    \  states[0],ibs[0],ii0
	    //a3=alw[3](j3,1)       |     |     | b0=alw[1](j1,0)            
	    //                      |----<--s0--|   op(1)        
	    //b2=alw[3](j3,0); njm2 |     |     | a1=alw[1](j1,1)
	    //states[2],ibs[2],ii2   \    v    /  states[1],ibs[1],ii1   
	    //                         ---<---		   
	    //                            njm1           
	    //            a2=alw[2](j2,1)      b1=alw[2](j2,0)             
	    TinyVector<int,Nd> js;
	    for (js[0]=0; js[0]<alw[0].extent(0); ++js[0]){
	      int a0 = cix.inBlock(abs(alw[0](js[0],1))); // alw[0] index is given in terms of atomic state index, but a0 is index within the block. cix.inBlock translates this.
	      int b3 = cix.inBlock(abs(alw[0](js[0],0))); // cix.inBlock translates again to index inside the block ib
	      double mm0 = sign(alw[0](js[0],1)) * s0_sgn * s1_sgn;  // matrix element of operator op(0) <alw[0](j0,1)|op(0)|alw[0](j0,0)> is stored in the sign of op(0) acting on index, hence sign of alw[0](js[0],1) contains op(0) matrix element.
	      for (js[1]=0; js[1]<alw[1].extent(0); ++js[1]){
		int b0 = cix.inBlock(abs(alw[1](js[1],0)));
		int _ii0_ = cix.Sgm_symmetry[ibs[0]](a0,b0); // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
		if (_ii0_==0) continue;
		int ii0 = abs(_ii0_)-1;
		int imp_ii0 = G_importance[ii0];
		int ii_sign0 = sign(_ii0_);                 // Green's function value for off-diagonal propagators could be negative/opposite, hence sign
		int a1 = cix.inBlock(abs(alw[1](js[1],1)));
		double mm1 = mm0 * sign(alw[1](js[1],1));  // matrix element: <alw[1](j1,1)|op(1)|alw[1](j1,0)>
		for (js[2]=0; js[2]<alw[2].extent(0); ++js[2]){
		  int b1 = cix.inBlock(abs(alw[2](js[2],0)));
		  int _ii1_ = cix.Sgm_symmetry[ibs[1]](a1,b1); // Green's function idex G(a1,b1)==G[|_ii1_|-1]
		  if (_ii1_==0) continue;
		  int ii1 = abs(_ii1_)-1;
		  int imp_ii1 = G_importance[ii1];
		  int ii_sign1 = sign(_ii1_) * ii_sign0;
		  int a2 = cix.inBlock(abs(alw[2](js[2],1)));
		  double mm2 = mm1 * sign(alw[2](js[2],1));
		  for (js[3]=0; js[3]<alw[3].extent(0); ++js[3]){
		    int b2 = cix.inBlock(abs(alw[3](js[3],0)));
		    int _ii2_ = cix.Sgm_symmetry[ibs[2]](a2,b2); // Green's function index G(a2,b2)==G[|_ii2_|-1]
		    if (_ii2_==0) continue;
		    int ii2 = abs(_ii2_)-1;
		    int imp_ii2 = G_importance[ii2];
		    int ii_sign2 = sign(_ii2_) * ii_sign1;
		    int a3 = cix.inBlock(abs(alw[3](js[3],1)));
		    double mm3 = mm2 * sign(alw[3](js[3],1));
		    int _ii3_ = cix.Sgm_symmetry[ibs[3]](a3,b3); // Green's function index G(a3,b3)==G[|_ii3_|-1]
		    if (_ii3_==0) continue;
		    int ii3 = abs(_ii3_)-1;
		    int importance = G_importance[ii3] + imp_ii2 + imp_ii1 + imp_ii0;
		    if (importance > OCA_cut_sum ) continue; // this cuts-off PP-states that are less important
		    int ii_sign3 = sign(_ii3_) * ii_sign2;
		    double mm = mm3*ii_sign3;                    // total sign due to off-diagonal sign (ii_sign3) and matrix elements sign mm3
		    TinyVector<int,4> states(ii0, ii1, ii2, ii3); // atomic states indices between [0,...cix.Npseudo_unique-1] array.
		    //
		    OCAKey st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // we pack the two hybridizations (s1,s0) and atomic states into unique index-vector.
		    prepare_bounded_insert(st);
		    auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // try to insert st-vector as a new key and mm_nds.size() as the value.
		    // If the key is not yet present, inserted=true and "it" points to its current place. If key is present, "it" points to its existing place.
		    int iind = it->second;                    // The key st has value iind, either just inserted, or found from before
		    if (inserted){
		      mm_nds[mm_nds_size] = TinyVector<double,5>(mm,0,0,0,0); // if it was just inserted, deque has matrix element for PPself-energy equal mm, and nothing for two particle vertices.
		      mm_nds_size+=1;
		    }else
		      mm_nds[iind] += TinyVector<double,5>(mm,0,0,0,0); // if we found old entry, which just add to it the current value of mm.
		    // We are done with finding unique contributions. Now we need sub-parts or common parts that appear in the inner-most loop.
		    // 
		    // 
		    // Below we are finding sub-parts of a diagram that is common to many terms just found above. This is done for efficiency reasons,
		    // so that inside loop integration does not need to repeat the same integration many times.
		    // This subpart contains only states[0],states[1], and hybridization s0.
		    // The last number 1 means that it is needed for PP-self-energy calculation.
		    // The plot would look like:
		    //      t=0
		    //      --->---                
		    //              \  states[0],ibs[0],ii0
		    //               | b0=alw[1](j1,0)            
		    // t3-----<-s0---| t1  op(1)        
		    //               | a1=alw[1](j1,1)
		    //              /  states[1],ibs[1],ii1   
		    //      |--<---
		    //      t2
		    if (build_legacy_catalog) {
		    TinyVector<int,4> sti0(states[0],states[1],s0,1);
		    auto [it0, inserted0] = index_of_inside0.try_emplace(sti0, static_cast<int>(unique_states_inside[0].size()));
		    if (inserted0) unique_states_inside[0].push_back(sti0); // first time we've seen this triple
		    else if (unique_states_inside[0][it0->second][3]==0){   // it could happen that it was already found below, but only for two-particle vertex calculation.
		      // Now we know that we need it also for PP self-energy calculation, hence we will change the last flag from 0 to 1.
		      //cout<<"WARN at sti0="<<unique_states_inside0[it0->second]<<" was set to 0"<<endl;
		      unique_states_inside[0][it0->second][3]=1;
		    }
		    // For second self-energy contribution we need different parts of the diagram in the inside loop. The part is sketched below:
		    //                           t3
		    //                            |
		    //                   t2 |    s1     |t=0  op(1)        
		    //b2=alw[3](j3,0); njm2 |     v     | a1=alw[1](j1,1)
		    //states[2],ibs[2],ii2   \    |    /  states[1],ibs[1],ii1   
		    //                         ---<---
		    //                           t1
		    TinyVector<int,4> sti1(states[1],states[2],s1,1);
		    auto [it1, inserted1] = index_of_inside1.try_emplace(sti1, static_cast<int>(unique_states_inside[1].size()));
		    if (inserted1) unique_states_inside[1].push_back(sti1); // first time we've seen this triple
		    else if (unique_states_inside[1][it1->second][3]==0){
		      //cout<<"WARN at sti1="<<unique_states_inside1[it1->second]<<" was set to 0"<<endl;
		      unique_states_inside[1][it1->second][3]=1;
		    }
		    // For the third self-energy contribution we need yet another part of the diagram in the inside loop. The sketch below:
		    //                               
		    //                         --->--t2
		    //states[3],ibs[3],ii3   /     
		    //a3=alw[3](j3,1)       |    
		    //                   t1 |-s0-<-------- t3
		    //b2=alw[3](j3,0); njm2 |     
		    //states[2],ibs[2],ii2   \    
		    //                         --<--t=0
		    TinyVector<int,4> sti2(states[2],states[3],s0,1);
		    auto [it2, inserted2] = index_of_inside2.try_emplace(sti2, static_cast<int>(unique_states_inside[2].size()));
		    if (inserted2) unique_states_inside[2].push_back(sti2); // first time we've seen this triple
		    else if (unique_states_inside[2][it2->second][3]==0){
		      //cout<<"WARN at sti2="<<unique_states_inside2[it2->second]<<" was set to 0"<<endl;
		      unique_states_inside[2][it2->second][3]=1;
		    }
		    // For the last forth self-energy we require this part of the diagram:
		    //                           t1
		    //                         --->---                
		    //states[3],ibs[3],ii3   /    |    \  states[0],ibs[0],ii0
		    //a3=alw[3](j3,1)       |     v     | b0=alw[1](j1,0)            
		    //                      |     s1    |   op(1)
		    //                     t=0    |    t2
		    //                            |
		    //                            t3
		    TinyVector<int,4> sti3(states[3],states[0],s1,1);
		    auto [it3, inserted3] = index_of_inside3.try_emplace(sti3, static_cast<int>(unique_states_inside[3].size()));
		    if (inserted3) unique_states_inside[3].push_back(sti3); // first time we've seen this triple
		    else if (unique_states_inside[3][it3->second][3]==0){
		      //cout<<"WARN at sti3="<<unique_states_inside3[it3->second]<<" was set to 0"<<endl;
		      unique_states_inside[3][it3->second][3]=1;
		    }
		    // Here we remember how to access results from unique_states_inside0, and  unique_states_inside1,..,
		    // when we are in the long loop going through all possible contibutions of states_unique
		    if (inserted){
		      indexi[iind][0] = it0->second;
		      indexi[iind][1] = it1->second;
		      indexi[iind][2] = it2->second;
		      indexi[iind][3] = it3->second;
		    }else{
		      // It could happen that we added st entry first when computing vertex F,
		      // and therefore we set index0[] and index1[] to -1.
		      // We now found the same st, and instead of inserting it, we just updated it and added mm to it.
		      // But therefore we need to check if index0[] or index1[] is set to -1, and correct.
		      if (indexi[iind][0]==-1) indexi[iind][0]=it0->second;
		      if (indexi[iind][1]==-1) indexi[iind][1]=it1->second;
		      if (indexi[iind][2]==-1) indexi[iind][2]=it2->second;
		      if (indexi[iind][3]==-1) indexi[iind][3]=it3->second;
		    }
		    }
		    // Below we are adding information needed for the two-particle vertex F(om) calculation.
		    // The simplest case of njmc being proportional to psi^+ operator is handled here. This is the case for
		    // one-dimensional block and sometimes other small blocks. But generically the two are not proportional, hence it
		    // must be handled separately with njmc1==2 below.
		    if (nx1.njm_c==1 and cmpPhysG){  // if njm is proportional to F, its simple: Add to the same iind entry, but nds1 part, i.e., second component
		      mm_nds[iind] += destroys_Cooper_pair_or_normal_down1 ? TinyVector<double,5>(0,0,0,mm*njmr1,0) : TinyVector<double,5>(0,mm*njmr1,0,0,0);
		    }
		    if (nx2.njm_c==1 and cmpPhysG){  // if njm is proportional to F its simple: Add to the same iind entry, but last component.
		      mm_nds[iind] += destroys_Cooper_pair_or_normal_down2 ? TinyVector<double,5>(0,0,0,0,mm*njmr2) : TinyVector<double,5>(0,0,mm*njmr2,0,0) ;
		    }
		  }
		  if (cmpPhysG and nx2.njm_c==2){
		    // WARNING: This can be improved. You don't need extra inside loop terms. Just rewrite with existing inside terms!!!!
		    // Here we handle two particle vertex F(om) calculation, and add all needed additional terms.
		    // Namely, b2 and a3 are not constrained by Pauli principle anymore, and more terms are possible, i.e., anything
		    // that does not vanish in njm2(b2,a3) matrix. We hence need to go through all possibilities for b2 and a3 indices
		    // in their respective blocks. But this has to be done outside js[3] loop, because js[3] loop also iterates through
		    // the same psi^+_2 operator, which we now replaced by njm2(b2,a3) matrix.
		    // Here is the sketch for this part of calculation:
		    //                         --->
		    //states[3],ibs[3],ii3   /  
		    //                   a3 | a3!=alw[3](j3,1)
		    //       njm2[b2,a3]    |--s0--<-- not op(3)
		    //                   b2 | b2!=alw[3](j3,0)
		    //states[2],ibs[2],ii2   \  
		    //                         ---<
		    for (int _b2_=0; _b2_<cix.Block[ibs[2]].size(); ++_b2_){
		      int __ii2__ = cix.Sgm_symmetry[ibs[2]](a2,_b2_);
		      if (__ii2__==0) continue;
		      int ii2 = abs(__ii2__)-1;
		      int imp_ii2 = G_importance[ii2];
		      for (int _a3_=0; _a3_<cix.Block[ibs[3]].size(); ++_a3_){
			if (njm2(_b2_,_a3_)==0) continue;
			int __ii3__ = cix.Sgm_symmetry[ibs[3]](_a3_,b3);
			if (__ii3__==0) continue;
			int ii3 = abs(__ii3__)-1;
			int importance = G_importance[ii3] + imp_ii2 + imp_ii1 + imp_ii0;
			if (importance > OCA_cut_sum ) continue;
			double mmp = ii_sign1*sign(__ii2__)*sign(__ii3__)*mm2; // sign from the all but op(2) (equal to mm2) as well as off-diagonal propagators.
			// This is the entire matrix element, which inclused njm2, instead of op(3) for kink 3, and for all other kinks
			// it implements op(i) operator, giving proper mm matrix element.
			double nds2 =  mmp*njm2(_b2_,_a3_);
			if (nds2==0) continue;
			TinyVector<int,4> states(ii0, ii1, ii2, ii3); // atomic states indices between [0,...cix.Npseudo_unique-1] array.
			OCAKey st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // all atributes of diagram contribution.
			prepare_bounded_insert(st);
			auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // Did we see that contribution before?
			int _iind_ = it->second; // in mm_nds and indexi, it is stored in _iind_
			if (inserted){// not seen before, hence insert into mm_nds, index0, index1, index2, index3
			  if (destroys_Cooper_pair_or_normal_down2){
			    mm_nds[mm_nds_size]=TinyVector<double,5>(0,0,0,0,nds2);
			  }else
			    mm_nds[mm_nds_size]=TinyVector<double,5>(0,0,nds2,0,0);
			  if (build_legacy_catalog) {
			  //index0[mm_nds_size]=-1; // should not be needed, hence we are not going to compute contribution of the inner loop.
			  //index1[mm_nds_size]=-1; // should not be needed, because inner loop for G&F uses only index2 and index3.
			  TinyVector<int,4> sti2(states[2],states[3],s0,1);// this inner loop is needed for G&F calculation
			  auto [_it2_, _inserted2_] = index_of_inside2.try_emplace(sti2, static_cast<int>(unique_states_inside[2].size()));
			  if (_inserted2_){
			    sti2[3]=0; // this marks that this is not needed for PP self-energies
			    unique_states_inside[2].push_back(sti2); // first time we've seen this triple
			  }
			  indexi[mm_nds_size][2]=_it2_->second;
			  TinyVector<int,4> sti3(states[3],states[0],s1,1);// this inner loop is also needed for G&F calculation.
			  auto [_it3_, _inserted3_] = index_of_inside3.try_emplace(sti3, static_cast<int>(unique_states_inside[3].size()));
			  if (_inserted3_){
			    sti3[3]=0; // this marks that this is not needed for PP self-energies
			    unique_states_inside[3].push_back(sti3); // first time we've seen this triple
			  }
			  indexi[mm_nds_size][3]=_it3_->second;
			  }
			  mm_nds_size += 1;
			}else // we have seen it before, hence just add nds2 contribution.
			  if (destroys_Cooper_pair_or_normal_down2)
			    mm_nds[_iind_] += TinyVector<double,5>(0,0,0,0,nds2);
			  else
			    mm_nds[_iind_] += TinyVector<double,5>(0,0,nds2,0,0);
		      }
		    }
		  }
		}
		if (cmpPhysG and nx1.njm_c==2){
		  // Here we handle two particle vertex F(om) calculation for the second type of orbital.
		  // The plot is
		  //                        b3
		  // states[3],ibs[3],ii3  /
		  //                   a3 |
		  //    op(3)  ~ js[3]    v-    |    -^   op(1)        
		  //    alw[3](j3,0) = b2 |     v     | a1=alw[1](j1,1)
		  // states[2],ibs[2],ii2  \  s1|    /  states[1],ibs[1],ii1   
		  //                      a2 ---<---  b1
		  //                        njm1[b1,a2]
		  //          a2!=alw[2](j2,1)      b1!=alw[2](j2,0)
		  // Namely, b1 and a2 are not constrained by Pauli principle anymore, and more terms are possible, i.e., anything
		  // that does not vanish in njm1(b1,a2) matrix. We hence need to go through all possibilities for b1 and a2 indices
		  // in their respective blocks. But this has to be done outside js[2] loop, because js[2] loop also iterates through
		  // the same psi^+_1==op(2) operator, which we now replaced by njm1(b1,a2) matrix. Because the way the loops are
		  // structured above, we will need to skip both js[2] and js[3], hence only js[1] is left.
		  // As a result, we need to iterate through js[3] loop inside this code.
		  for (js[3]=0; js[3]<alw[3].extent(0); ++js[3]){
		    int b2 = cix.inBlock(abs(alw[3](js[3],0)));
		    int a3 = cix.inBlock(abs(alw[3](js[3],1)));
		    int _ii3_ = cix.Sgm_symmetry[ibs[3]](a3,b3); // Green's function index G(a3,b3)==G[|_ii3_|-1]
		    if (_ii3_==0) continue;
		    int ii3 = abs(_ii3_)-1;
		    int imp_ii3 = G_importance[ii3];
		    int ii_sign_03 = sign(_ii3_) * ii_sign0;  // sign due to propagator states[0], states[3]. Later we will add states[1] and states[2].
		    double mm3 = mm1 * sign(alw[3](js[3],1)); // sign due to op(0),op(1) and op(3). Later we will add njm1
		    for (int _b1_=0; _b1_<cix.Block[ibs[1]].size(); ++_b1_){ // now going over all other possibilities
		      int __ii1__ = cix.Sgm_symmetry[ibs[1]](a1,_b1_); // a1 is the same, but b1 can be anything within the block
		      if (__ii1__==0) continue;
		      int ii1 = abs(__ii1__)-1;
		      int imp_ii1 = G_importance[ii1];
		      for (int _a2_=0; _a2_<cix.Block[ibs[2]].size(); ++_a2_){
			if (njm1(_b1_,_a2_)==0) continue;
			int __ii2__ = cix.Sgm_symmetry[ibs[2]](_a2_,b2);
			if (__ii2__==0) continue;
			int ii2 = abs(__ii2__)-1;
			int importance = G_importance[ii2] + imp_ii1 + imp_ii3 + imp_ii0;
			if (importance > OCA_cut_sum ) continue;
			double mmp = ii_sign_03*sign(__ii1__)*sign(__ii2__)*mm3;
			// nds1 is the entire matrix element, which inclused njm1, instead of op(2) for kink 2, and for all other kinks
			// it implements op(i) operator, giving proper mm matrix element. See commen for mm3 above.
			double nds1 =  mmp*njm1(_b1_,_a2_); // Now adding njm1.
			TinyVector<int,4> states(ii0, ii1, ii2, ii3); // atomic states indices between [0,...cix.Npseudo_unique-1] array.
			OCAKey st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // all atributes of diagram contribution.
			prepare_bounded_insert(st);
			auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // Did we see that contribution before?
			int _iind_ = it->second; // in mm_nds and indexi, it is stored in _iind_
			if (inserted){ // not seen before, hence insert into mm_nds, index0, index1, index2, index3
			  if (destroys_Cooper_pair_or_normal_down1)
			    mm_nds[mm_nds_size]=TinyVector<double,5>(0,0,0,nds1,0);
			  else
			    mm_nds[mm_nds_size]=TinyVector<double,5>(0,nds1,0,0,0);
			  if (build_legacy_catalog) {
			  //index0[mm_nds_size]=-1; // should not be needed, hence we are not going to compute contribution of the inner loop.
			  //index1[mm_nds_size]=-1; // should not be needed, because inner loop for G&F uses only index2 and index3.
			  TinyVector<int,4> sti2(states[2],states[3],s0,1); // this inner loop is needed for G&F calculation
			  auto [_it2_, _inserted2_] = index_of_inside2.try_emplace(sti2, static_cast<int>(unique_states_inside[2].size()));
			  if (_inserted2_){
			    sti2[3]=0; // this marks that this is not needed for PP self-energies
			    unique_states_inside[2].push_back(sti2); // first time we've seen this triple
			  }
			  indexi[mm_nds_size][2] = _it2_->second;
			  TinyVector<int,4> sti3(states[3],states[0],s1,1); // this inner loop is needed for G&F calculation
			  auto [_it3_, _inserted3_] = index_of_inside3.try_emplace(sti3, static_cast<int>(unique_states_inside[3].size()));
			  if (_inserted3_){
			    sti3[3]=0; // this marks that this is not needed for PP self-energies
			    unique_states_inside[3].push_back(sti3); // first time we've seen this triple
			  }
			  indexi[mm_nds_size][3] = _it3_->second;
			  }
			  mm_nds_size+=1;
			}else // we have seen it before, hence just add nds1 contribution.
			  if (destroys_Cooper_pair_or_normal_down1)
			    mm_nds[_iind_] += TinyVector<double,5>(0,0,0,nds1,0);
			  else
			    mm_nds[_iind_] += TinyVector<double,5>(0,nds1,0,0,0);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if (task_cursor!=task_owner.size())
      throw std::logic_error("OCA pass 2 task ordering disagrees with scheduler discovery");
  }

  if (bounded) {
    flush_bounded();
    // No path gather, plan construction, or collective is allowed here.
    // Each rank independently spools its owned reduced bath pairs.
    return OCAEnumeration{};
  }
  const auto pass2_end=EnumerationClock::now();
  (*mpi.out)<<"local mm_nds_size="<<mm_nds_size
            <<" rank="<<mpi.rank<<endl;
  OCAEnumeration result;
  result.includes_physical=cmpPhysG;
  result.npseudo=cix.Npseudo_unique;
  result.nbath_unique=cix.Nbath_unique_unique;
  result.cut_sum=OCA_cut_sum;
  result.importance=G_importance;

  // Convert this rank's node-based construction map into one contiguous,
  // key-sorted representation.  inner_index is deliberately not preserved:
  // those insertion-order indices are local to a rank and are invalid after
  // MPI combines paths found by different ranks.
  const auto compact_start=EnumerationClock::now();
  vector<OCACompactPath> local_paths;
  local_paths.reserve(states_unique.size());
  for (const auto& [_st_, iind] : states_unique) {
    OCACompactPath path;
    path.key=_st_;
    for (int j=0; j<5; ++j) path.weights[j]=mm_nds[iind][j];
    path.inner_index={-1,-1,-1,-1};
    local_paths.push_back(path);
  }
  // Packed records now own the values: release construction copies early.
  decltype(states_unique)().swap(states_unique);
  decltype(mm_nds)().swap(mm_nds);
  std::sort(local_paths.begin(),local_paths.end(),
            [](const OCACompactPath& a,const OCACompactPath& b) {
              return a.key<b.key;
            });
  // Every packed key contains its complete reduced bath pair, and the
  // scheduler assigns that pair to exactly one rank.  Exact cancellations
  // have therefore already met locally; no cross-rank reduction is needed.
  local_paths.erase(
    std::remove_if(local_paths.begin(),local_paths.end(),
                   [](const OCACompactPath& path) {
                     return OCAPathWeightsAreZero(path.weights);
                   }),
    local_paths.end());

  unsigned long long global_nonzero=local_paths.size();
#ifdef _MPI
  if (mpi.size>1) {
    const unsigned long long local_nonzero=local_paths.size();
    MPI_Allreduce(&local_nonzero,&global_nonzero,1,MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM,MPI_COMM_WORLD);
  }
#endif
  result.global_nonzero_path_count=global_nonzero;
  // Measurements on the large SC benchmark show that pair catalogs, sparse
  // triplets, and the two compiled plans together require about 4.5 times the
  // compact path bytes at peak.  This conservative estimate selects the
  // memory-safe path before any global Allgatherv allocation is attempted.
  constexpr double compiled_plan_factor=4.5;
  result.estimated_replicated_megabytes=
    compiled_plan_factor*global_nonzero*sizeof(OCACompactPath)/(1024.0*1024.0);
  // The spline integrator addresses one globally rebuilt inner catalog and
  // consequently always requires replicated paths.  Diagnostics that run it
  // retain the old exact behavior even if distributed DLR was requested.
  result.topology_shard_count=build_legacy_catalog ? 1 :
    OCAChooseTopologyShardCount(
      mpi_topology_mode,result.estimated_replicated_megabytes,
      max_replicated_megabytes,mpi.size);
  result.paths_distributed=result.topology_shard_count>1;
  OCASetTopologyGroupMetadata(result,mpi);

  const auto merge_start=EnumerationClock::now();
#ifdef _MPI
  if (result.paths_distributed)
    result.compact_paths=std::move(local_paths);
  else
    result.compact_paths=GatherOCAPathsMPI(local_paths,mpi);
#else
  result.compact_paths=std::move(local_paths);
#endif

  // All ranks now hold the same sorted, globally accumulated path list.
  // Recreate compact convolution catalogs from physical indices only when
  // the legacy spline integrator will consume them.  DLR builds independent
  // sparse plans from OCAPathRecords and never reads this catalog.
  if (build_legacy_catalog) RebuildOCAInnerCatalog(result);
  const auto merge_end=EnumerationClock::now();

  // A deterministic checksum and five independent weight sums make the
  // forthcoming MPI enumeration easy to compare against this serial
  // reference.  Hash double bit patterns, not formatted decimal values.
  // XOR of per-record FNV hashes is independent of record order and MPI
  // partition.  This makes replicated, bath-pair-sharded, and rebalanced
  // enumerations directly comparable in the log.
  uint64_t checksum=0;
  constexpr uint64_t fnv_prime=1099511628211ULL;
  std::array<long double,5> weight_sums{};
  for (const auto& path : result.compact_paths) {
    uint64_t record_hash=1469598103934665603ULL;
    const auto checksum_word=[&](uint64_t word) {
      record_hash^=word;
      record_hash*=fnv_prime;
    };
    checksum_word(path.key.low);
    checksum_word(path.key.high);
    for (int j=0; j<5; ++j) {
      checksum_word(std::bit_cast<uint64_t>(path.weights[j]));
      weight_sums[j]+=path.weights[j];
    }
    checksum^=record_hash;
  }
#ifdef _MPI
  if (result.paths_distributed && mpi.size>1) {
    uint64_t global_checksum=0;
    MPI_Allreduce(&checksum,&global_checksum,1,MPI_UINT64_T,MPI_BXOR,
                  MPI_COMM_WORLD);
    checksum=global_checksum;
    std::array<long double,5> global_weight_sums{};
    MPI_Allreduce(weight_sums.data(),global_weight_sums.data(),5,
                  MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    weight_sums=global_weight_sums;
  }
#endif

  if (export_path_records){
    // compact_paths is already the fully reduced set after symmetry and CIX
    // selection rules.  Exporting here guarantees that the DLR evaluator and
    // legacy integrator consume identical accumulated weights.
    result.paths.reserve(result.compact_paths.size());
    for (const auto& path : result.compact_paths){
      const TinyVector<int,6> st=unpack6(path.key);
      OCAPathRecord record;
      // unpack6 stores baths first and pseudo states in reverse ring order.
      // Convert once to the p0,p1,p2,p3 convention documented in
      // OCAPathRecord and used by every compact cut table.
      record.pseudo={st[5],st[4],st[3],st[2]};
      record.s0=st[1];
      record.s1=st[0];
      record.weights=path.weights;
      result.paths.push_back(record);
    }
#ifdef _MPI
    if (result.paths_distributed) {
      // Enumeration itself used all world ranks.  Consolidate those unique
      // streams into the memory-selected number of topology shards, preserving
      // whole reduced bath pairs whenever possible, then broadcast each shard
      // only to the ranks that will divide its AAA poles.
      const bool already_pair_partitioned=
        result.topology_shard_count==mpi.size &&
        mpi.size<=static_cast<int>(pair_costs.size());
      DistributeOCAPathRecordsHybridMPI(
        result.paths,result.nbath_unique,result.topology_shard_count,
        already_pair_partitioned,mpi);
    }
#endif
  }

  // Only globally reduced compact paths and their reconstructed convolution
  // catalog survive the function.  Construction hash nodes, buckets, local
  // insertion indices, and local inner catalogs are released on return.
  const auto compact_end=EnumerationClock::now();

  double pass1_seconds=
    std::chrono::duration<double>(pass1_end-pass1_start).count();
  double pass2_seconds=
    std::chrono::duration<double>(pass2_end-pass2_start).count();
  double merge_seconds=
    std::chrono::duration<double>(merge_end-merge_start).count();
  double compact_seconds=
    std::chrono::duration<double>(compact_end-compact_start).count();
  unsigned long long local_count=static_cast<unsigned long long>(mm_nds_size);
  unsigned long long local_count_sum=local_count;
  unsigned long long local_count_min=local_count;
  unsigned long long local_count_max=local_count;
#ifdef _MPI
  if (mpi.size>1) {
    double local_times[4]={pass1_seconds,pass2_seconds,
                           merge_seconds,compact_seconds};
    double maximum_times[4]={};
    MPI_Allreduce(local_times,maximum_times,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    pass1_seconds=maximum_times[0];
    pass2_seconds=maximum_times[1];
    merge_seconds=maximum_times[2];
    compact_seconds=maximum_times[3];
    MPI_Allreduce(&local_count,&local_count_sum,1,MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&local_count,&local_count_min,1,MPI_UNSIGNED_LONG_LONG,
                  MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&local_count,&local_count_max,1,MPI_UNSIGNED_LONG_LONG,
                  MPI_MAX,MPI_COMM_WORLD);
  }
#endif
  if (mpi.rank==mpi.master) {
    (*mpi.out) << "OCA enumeration MPI: local_path_count_range=["
               << local_count_min << "," << local_count_max << "]"
               << " local_path_count_sum=" << local_count_sum
               << " global_nonzero_paths=" << result.global_nonzero_path_count
               << " topology="
               << (result.paths_distributed ? "hybrid-distributed" : "replicated")
               << " topology_shards=" << result.topology_shard_count
               << " workers_per_shard_range=["
               << mpi.size/result.topology_shard_count << ","
               << (mpi.size+result.topology_shard_count-1)/
                    result.topology_shard_count << "]"
               << " estimated_replicated_mb="
               << result.estimated_replicated_megabytes
               << " checksum=0x"
               << std::hex << checksum << std::dec << " weight_sums=[";
    for (int j=0; j<5; ++j)
      (*mpi.out) << (j ? "," : "") << static_cast<double>(weight_sums[j]);
    (*mpi.out) << "] pass1_max_seconds=" << pass1_seconds
               << " pass2_max_seconds=" << pass2_seconds
               << " merge_max_seconds=" << merge_seconds
               << " compact_max_seconds=" << compact_seconds
               << endl;
  }
  return result;
}

// Legacy uniform-mesh OCA integrator.  All discrete diagram work has already
// been done by EnumerateOCAPaths; this routine now contains numerical
// integration only (plus consistency/debug checks on the supplied tables).
inline double I_OCA(const OCAEnumeration& enumeration,
		    Array<double,2>& Sm,
		    Array<double,2>& G00t,
		    Array<double,2>& F00t,
		    Array<double,2>& F01t,
		    bool cmpPhysG,
		    double Q,
		    double beta,
		    const Array<double,1>& tau,
		    const Spline1DMany& Gms,
		    const Spline1DMany& Deltas,
		    Cix& cix,
		    json& pars,
		    MPIcls& mpi,
		    std::array<Array<double,2>,4>* cut_references=nullptr)
{
  bool debug=false;
  const double coef[]={17/48.,59/48.,43/48.,49./48.,1.0}; // O(dh^4)
  const int Nd=4;

  const auto& compact_paths=enumeration.compact_paths;
  const auto& unique_states_inside=enumeration.unique_inner_states;

  if (fabs(beta-tau(tau.size()-1))>1e-7)
    cerr << "ERROR in I_OCA: last point in tau=" << tau(tau.size()-1)
         << " should be beta=" << beta << endl;

  // Everything below this point is legacy-only.  In particular, pure DLR
  // mode never allocates this uniform mesh and never evaluates the input
  // splines on it.
  const int Nt=tau.size()*pars["r_OCA_mesh"].get<double>();
  const double dh=beta/(Nt-1.0);
  Array<double,2> Gm(Gms.size(),Nt),Delta(Deltas.size(),Nt);
  int pos=0;
  for (int it=0; it<Nt; ++it){
    const intpar ips=InterpLeft(dh*it,pos,tau);
    for (int a=0; a<Gms.size(); ++a) Gm(a,it)=Gms(ips,a);
    for (int b=0; b<Deltas.size(); ++b) Delta(b,it)=Deltas(ips,b);
  }

  // The uniform OCA mesh ends at beta.  A tiny discrepancy here is harmless
  // interpolation roundoff; report only a meaningful mismatch and only once.
  constexpr double endpoint_tolerance=1e-6;
  if (mpi.rank==mpi.master){
    for (int a=0; a<Gms.size(); ++a){
      if (fabs(Gm(a,Nt-1)-Gms(tau.size()-1,a))>endpoint_tolerance)
        cerr << "WARNING in I_OCA: last point Gm(" << a << "," << Nt-1
             << ") differs from Gms by more than " << endpoint_tolerance
             << endl;
    }
  }

  {
    int ijk=0;
    for (const auto& path : compact_paths){
      TinyVector<int,6> st = unpack6(path.key);
      double mm = path.weights[0];
      double nds1 = abs(path.weights[1])+abs(path.weights[3]);
      double nds2 = abs(path.weights[2])+abs(path.weights[4]);
      int s1 = st[0];
      int s0 = st[1];
      TinyVector<int,Nd> states(st[5],st[4],st[3],st[2]);
      int index_of_0 = path.inner_index[0];
      int index_of_1 = path.inner_index[1];
      int index_of_2 = path.inner_index[2];
      int index_of_3 = path.inner_index[3];
      if (mm!=0){
	TinyVector<int,4> sti0 = unique_states_inside[0][index_of_0];
	if (sti0[0]!=states[0] or sti0[1]!=states[1] or sti0[2]!=s0){
	  cout<<"WARNING0 sti0="<<sti0<<" while states[0],states[1],s0="<<states[0]<<","<<states[1]<<","<<s0<<" and st="<<st<<endl;
	}
	TinyVector<int,4> sti1 = unique_states_inside[1][index_of_1];
	if (sti1[0]!=states[1] or sti1[1]!=states[2] or sti1[2]!=s1){
	  cout<<"WARNING1 sti1="<<sti1<<" while states[1],states[2],s1="<<states[1]<<","<<states[2]<<","<<s1<<" and st="<<st<<endl;
	}
      }
      TinyVector<int,4> sti2 = unique_states_inside[2][index_of_2];
      if (sti2[0]!=states[2] or sti2[1]!=states[3] or sti2[2]!=s0){
	cout<<"WARNING2 sti2="<<sti2<<" while states[2],states[3],s0="<<states[2]<<","<<states[3]<<","<<s0<<" and st="<<st<<endl;
      }
      TinyVector<int,4> sti3 = unique_states_inside[3][index_of_3];
      if (sti3[0]!=states[3] or sti3[1]!=states[0] or sti3[2]!=s1){
	cout<<"WARNING3 sti3="<<sti3<<" while states[3],states[0],s1="<<states[3]<<","<<states[0]<<","<<s1<<" and st="<<st<<endl;
      }
      if (index_of_2<0){
	cout << "ERROR id[2]="<<index_of_2<<" at st="<<st<<" compact_index="<<ijk<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
      }
      if (index_of_3<0){
	cout << "ERROR id[3]="<<index_of_3<<" at st="<<st<<" compact_index="<<ijk<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
      }
      ++ijk;
    }
  }
  
  if (debug){
    // Below we print out all contributions for OCA, i.e., all possible triple integrations we need to do.
    for (int i=0; i<unique_states_inside[0].size(); ++i){
      TinyVector<int,4> sti0 = unique_states_inside[0][i];
      fmt::print(cout,"#R0: {} ii0={} ii1={} s0={} pp={}\n", i, sti0[0]+1, sti0[1]+1, sti0[2], sti0[3]);
    }
    for (int i=0; i<unique_states_inside[1].size(); ++i){
      TinyVector<int,4> sti1 = unique_states_inside[1][i];
      fmt::print(cout,"#R1: {} ii1={} ii2={} s1={} pp={}\n", i, sti1[0]+1, sti1[1]+1, sti1[2], sti1[3]);
    }
    for (int i=0; i<unique_states_inside[2].size(); ++i){
      TinyVector<int,4> sti2 = unique_states_inside[2][i];
      fmt::print(cout,"#R2: {} ii2={} ii3={} s0={} pp={}\n", i, sti2[0]+1, sti2[1]+1, sti2[2], sti2[3]);
    }
    for (int i=0; i<unique_states_inside[3].size(); ++i){
      TinyVector<int,4> sti3 = unique_states_inside[3][i];
      fmt::print(cout,"#R3: {} ii3={} ii0={} s1={} pp={}\n", i, sti3[0]+1, sti3[1]+1, sti3[2], sti3[3]);
    }
    int ijk=0;
    for (const auto& path : compact_paths){
      TinyVector<int,6> st = unpack6(path.key);
      double mm = path.weights[0];
      double nds1 = path.weights[1];//abs(path.weights[1])+abs(path.weights[3]);
      double nds2 = path.weights[2];//abs(path.weights[2])+abs(path.weights[4]);
      if (fabs(path.weights[3])>1e-7 or fabs(path.weights[4])>1e-7) cout<<"ERROR expecting nds3="<<path.weights[3]<<" and nds4="<<path.weights[4]<<" to vanish"<<endl;
      int s1 = st[0];
      int s0 = st[1];
      TinyVector<int,Nd> states(st[5],st[4],st[3],st[2]);
      int index_of_0 = path.inner_index[0];
      int index_of_1 = path.inner_index[1];
      if (mm!=0){
	TinyVector<int,4> sti0 = unique_states_inside[0][index_of_0];
	if (sti0[0]!=states[0] or sti0[1]!=states[1] or sti0[2]!=s0){
	  cout<<"WARNING0 sti0="<<sti0<<" while states[0],states[1],s0="<<states[0]<<","<<states[1]<<","<<s0<<" and st="<<st<<endl;
	}
	TinyVector<int,4> sti1 = unique_states_inside[1][index_of_1];
	if (sti1[0]!=states[1] or sti1[1]!=states[2] or sti1[2]!=s1){
	  cout<<"WARNING1 sti1="<<sti1<<" while states[1],states[2],s1="<<states[1]<<","<<states[2]<<","<<s1<<" and st="<<st<<endl;
	}
      }
      int index_of_2 = path.inner_index[2];
      TinyVector<int,4> sti2 = unique_states_inside[2][index_of_2];
      if (sti2[0]!=states[2] or sti2[1]!=states[3] or sti2[2]!=s0){
	cout<<"WARNING2 sti2="<<sti2<<" while states[2],states[3],s0="<<states[2]<<","<<states[3]<<","<<s0<<" and st="<<st<<endl;
      }
      int index_of_3 = path.inner_index[3];
      TinyVector<int,4> sti3 = unique_states_inside[3][index_of_3];
      if (sti3[0]!=states[3] or sti3[1]!=states[0] or sti3[2]!=s1){
	cout<<"WARNING3 sti3="<<sti3<<" while states[3],states[0],s1="<<states[3]<<","<<states[0]<<","<<s1<<" and st="<<st<<endl;
      }
      fmt::print(cout," {:6d} ii0={:<3d},ii1={:<3d},ii2={:<3d},ii3={:<3d},s0={:d},s1={:d}: mm={:2g} nds1={:3g},{:3g} nds0={:3g},{:3g}  index0={:4d} index1={:4d} index2={:4d} index3={:4d} ",
		 ijk,st[5]+1,st[4]+1,st[3]+1,st[2]+1,st[1],st[0],
		 mm, path.weights[1],path.weights[3],path.weights[2],path.weights[4],
		 index_of_0,index_of_1,index_of_2,index_of_3);
      ijk+=1;
      /*
      fmt::print(cout,"[");
      for (int l=0; l<cix.Sgm_indx[states[0]].size(); ++l){
	auto v = cix.Sgm_indx[states[0]][l];
	int ib = abs(v[0]);
	fmt::print(cout, "({},{}),", cix.Block[ib](v[1]), cix.Block[ib](v[2]) );
      }
      fmt::print(cout,"] [");
      for (int l=0; l<cix.Sgm_indx[states[1]].size(); ++l){
	auto v = cix.Sgm_indx[states[1]][l];
	int ib = abs(v[0]);
	fmt::print(cout, "({},{}),", cix.Block[ib](v[1]), cix.Block[ib](v[2]) );
      }
      fmt::print(cout,"] [");
      for (int l=0; l<cix.Sgm_indx[states[2]].size(); ++l){
	auto v = cix.Sgm_indx[states[2]][l];
	int ib = abs(v[0]);
	fmt::print(cout, "({},{}),", cix.Block[ib](v[1]), cix.Block[ib](v[2]) );
      }
      fmt::print(cout,"] [");
      for (int l=0; l<cix.Sgm_indx[states[3]].size(); ++l){
	auto v = cix.Sgm_indx[states[3]][l];
	int ib = abs(v[0]);
	fmt::print(cout, "({},{}),", cix.Block[ib](v[1]), cix.Block[ib](v[2]) );
      }
      fmt::print(cout,"]");
      */
      fmt::print(cout,"\n");
    }
  }

  // Below we decide how to redistribute jobs across the mpi run.
  std::vector<BlockJob> jobs;
  jobs.reserve(Nt+Nt); // Parallelization with respect to outside times for the two parts, computing PPSigma, and Gimp.
  for (int it3=0; it3<Nt; ++it3){ // this is for computing PPSigma
    double x = it3/(Nt+0.0); // we can estimate complexity of integrations
    jobs.push_back({it3, 0, x*x+0.022*x+0.002}); // the formula x^2 comes from the integration in the inside loop, and 0.022*x from outside loop
  }
  for (int it2=0; it2<Nt; ++it2){ // this is for computing Gimp&Fimp
    double x = it2/(Nt+0.0); // we can estimate complexity of integrations
    jobs.push_back({it2+Nt, 1, 1.7*(1-x)*(x+0.08)}); // the formula (1-x)*x comes from the integration in the inside loop, and (1-x)*0.08 from outside loop
  }
  std::vector<std::vector<BlockJob>> assign(mpi.size);
  if (mpi.size==1)
    assign[0]=jobs;
  else
    mpi.Assign_jobs(assign, jobs); // Here we actually assign jobs.
  
  // Below the actual time integration is performed.
  
  fsec t_ins=zero_t, t_out=zero_t;
  // For convenience of more efficient MPI communication, we transpose pseudo Self-energy 
  Array<double,2> Smc(Nt,cix.Npseudo_unique);
  std::array<Array<double,2>,4> Smc_cut;
  if (cut_references)
    for (auto& cut:Smc_cut){ cut.resize(Nt,cix.Npseudo_unique); cut=0; }
  double Phi=0;

  {  
    // contribution from the inside loop is contained in dSmwi
    vector<vector<Array<double,1>>> dSmwi(Nd);
    dSmwi[0].resize(unique_states_inside[0].size());
    dSmwi[1].resize(unique_states_inside[1].size());
    dSmwi[2].resize(unique_states_inside[2].size());
    dSmwi[3].resize(unique_states_inside[3].size());
    for (int i=0; i<unique_states_inside[0].size(); ++i){ dSmwi[0][i].resize(Nt); dSmwi[0][i]=0; }
    for (int i=0; i<unique_states_inside[1].size(); ++i){ dSmwi[1][i].resize(Nt); dSmwi[1][i]=0; }
    for (int i=0; i<unique_states_inside[2].size(); ++i){ dSmwi[2][i].resize(Nt); dSmwi[2][i]=0; }
    for (int i=0; i<unique_states_inside[3].size(); ++i){ dSmwi[3][i].resize(Nt); dSmwi[3][i]=0; }
  
    TinyVector<double,4> Smw;
    
    Smc=0;
    {  
      for (int ij=0; ij<assign[mpi.rank].size(); ++ij){
	if (assign[mpi.rank][ij].n!=0) continue;
	int it3 = assign[mpi.rank][ij].id;
	double dz=dh*intg_dh(it3,Nt-1);
	auto t0 = high_resolution_clock::now();      
	//                          op(0)                                                        op(0)				  	    
	//          alw[0](j0,0) =  b3   a0 = alw[0](j0,1)	               alw[0](j0,0) =  b3   a0 = alw[0](j0,1)		  
	//                             t=0			                                   t3				  
	//                            -->---        		                                        --->-          			  
	//Smc[states[3]](b3,a3)(t3)   s1|    \  dSmwi[0]		           states[3],ibs[3]   /  s1|       Smc[states[0]](b0,a0)(t3)  	    	            
	//      alw[3](j3,1) = a3   s0  |  s0 | b0		             alw[3](j3,1) = a3 |s0   |   s0  b0 = alw[1](j1,0)        
	//                      t3|----<------|t1  dSmwi[0]	                                   t2|----<------|t=0 			  	        	    
	//      alw[3](j3,0) = b2 |     v     | a1		                   alw[3](j3,0) = b2 |     v     | a1			  
	//      states[2],ibs[2]   \  s1|    /  dSmwi[0]		                   dSmwi[1]   \  s1|    /  dSmwi[1]		    	    	            	
	//                           ---<---    			                                ---<---				  
	//                             t2	    			                              a2  t1   b1				  
	//           alw[2](j2,1) = a2    b1 = alw[2](j2,0)	                               dSmwi[1]				  
	//							                                                                    
	//							                                                                    
	//                          op(0)                                                             op(0)				  
	//          alw[0](j0,0) =  b3   a0 = alw[0](j0,1)	                                   dSmwi[3]			  
	//                              t2			                                  b3  t1	a0		    	  
	//                           --->---        		                                           --->---        		  
	//              dSmwi[2]   /  s1|    \  states[0],ibs[0]	                      dSmwi[3]   /  s1|    \    dSmwi[3]		  
	//                     a3 |s0   |   s0| b0 = alw[1](j1,0)               alw[3](j3,1) = a3 |s0   |   s0| b0 = alw[1](j1,0)     
	//         dSmwi[2]     t1|----<------|t3   op(1)        	                             t=0|----<------|t2   op(1)        	  
	//                     b2 |     v     | a1 = alw[1](j1,1)	              alw[3](j3,0) = b2 |     v     | a1 = alw[1](j1,1)	  
	//              dSmwi[2]   \  s1|     Smc[states[1]](b1,a1)(t3)         Smc[states[2]](b2,a2)(t3)   s1|    /  states[1],ibs[1]	  
	//                           ---<-  			                                             -<---			  
	//                             t=0			                                             t3			    	  
	//           alw[2](j2,1) = a2    b1 = alw[2](j2,0)	                          alw[2](j2,1) = a2    b1 = alw[2](j2,0)           
	//
	// Precomputing the inside loop using unique_states_inside0
	for (int i=0; i<unique_states_inside[0].size(); ++i){
	  TinyVector<int,4> st = unique_states_inside[0][i];
	  if (st[3]==0) continue; // this means that this term is only needed for cmpPhysG
	  for (int it2=0; it2<=it3; it2++){
	    if (it2==0 && it3==0) continue;
	    int istart = it2!=0 ? 0 : 1;
	    double dsum=0;
	    for (int it1=istart; it1<=it2; it1++){
	      double dx=dh*intg_dh(it1,it2);
	      dsum += dx* Gm(st[0],it1) * Gm(st[1],it2-it1) * Delta(st[2],it3-it1);
	    }
	    dSmwi[0][i](it2) = dsum;
	  }
	}
	// Precomputing the inside loop using unique_states_inside1
	for (int i=0; i<unique_states_inside[1].size(); ++i){
	  TinyVector<int,4> st = unique_states_inside[1][i];
	  if (st[3]==0) continue; // this means that this term is only needed for cmpPhysG
	  for (int it2=0; it2<=it3; it2++){
	    if (it2==0 && it3==0) continue;
	    int istart = it2!=0 ? 0 : 1;
	    double dsum=0;
	    for (int it1=istart; it1<=it2; it1++){
	      double dx=dh*intg_dh(it1,it2);
	      dsum += dx* Gm(st[0],it1) * Gm(st[1],it2-it1) * Delta(st[2],Nt-1-(it3-it1));
	    }
	    dSmwi[1][i](it2) = dsum;
	  }
	}
	// Precomputing the inside loop using unique_states_inside2
	for (int i=0; i<unique_states_inside[2].size(); ++i){
	  TinyVector<int,4> st = unique_states_inside[2][i];
	  if (st[3]==0) continue;
	  for (int it2=0; it2<=it3; it2++){
	    if (it2==0 && it3==0) continue;
	    int istart = it2!=0 ? 0 : 1;
	    double dsum=0;
	    for (int it1=istart; it1<=it2; it1++){
	      double dx=dh*intg_dh(it1,it2);
	      dsum += dx* Gm(st[0],it1) * Gm(st[1],it2-it1) * Delta(st[2],Nt-1-(it3-it1));
	    }
	    dSmwi[2][i](it2) = dsum;
	  }
	}
	// Precomputing the inside loop using unique_states_inside3
	for (int i=0; i<unique_states_inside[3].size(); ++i){
	  TinyVector<int,4> st = unique_states_inside[3][i];
	  if (st[3]==0) continue;
	  for (int it2=0; it2<=it3; it2++){
	    if (it2==0 && it3==0) continue;
	    int istart = it2!=0 ? 0 : 1;
	    double dsum=0;
	    for (int it1=istart; it1<=it2; it1++){
	      double dx=dh*intg_dh(it1,it2);
	      dsum += dx* Gm(st[0],it1) * Gm(st[1],it2-it1) * Delta(st[2],it3-it1);
	    }
	    dSmwi[3][i](it2) = dsum;
	  }
	}
	auto t1 = high_resolution_clock::now();
	t_ins += t1-t0;
      
	for (const auto& path : compact_paths) {
	  TinyVector<int,6> st = unpack6(path.key);
	  double mm = -path.weights[0]; // minus comes from the crossing of hybridization lines
	  if (mm==0) continue;
	  TinyVector<int,4> states(st[5],st[4],st[3],st[2]);
	  int s0=st[1], s1=st[0];
	  int ind0=path.inner_index[0], ind1=path.inner_index[1], ind2=path.inner_index[2], ind3=path.inner_index[3];
	  Smw=0;
	  for (int it2=0; it2<=it3; it2++){
	    if (it2==0 && it3==0) continue;
	    double dy=dh*intg_dh(it2,it3);
	    int m_t2 = Nt-1-it2;
	    Smw[0] += dy* dSmwi[0][ind0](it2) * Gm(states[2],it3-it2) * Delta(s1,it2);
	    Smw[1] += dy* dSmwi[1][ind1](it2) * Gm(states[3],it3-it2) * Delta(s0,it2);
	    Smw[2] += dy* dSmwi[2][ind2](it2) * Gm(states[0],it3-it2) * Delta(s1,m_t2);
	    Smw[3] += dy* dSmwi[3][ind3](it2) * Gm(states[1],it3-it2) * Delta(s0,m_t2);
	  }
	  int mit3 = Nt-1-it3;
	  const double sigma_p3=Smw[0]*mm/cix.deg_pp(states[3]);
	  Smc(it3,states[3]) += sigma_p3;
	  const double sigma_p0=Smw[1]*mm/cix.deg_pp(states[0]);
	  const double sigma_p1=Smw[2]*mm/cix.deg_pp(states[1]);
	  const double sigma_p2=Smw[3]*mm/cix.deg_pp(states[2]);
	  Smc(it3,states[0]) += sigma_p0;
	  Smc(it3,states[1]) += sigma_p1;
	  Smc(it3,states[2]) += sigma_p2;
	  if (cut_references){
	    Smc_cut[3](it3,states[3]) += sigma_p3;
	    Smc_cut[0](it3,states[0]) += sigma_p0;
	    Smc_cut[1](it3,states[1]) += sigma_p1;
	    Smc_cut[2](it3,states[2]) += sigma_p2;
	  }
	  Phi += 0.25*dz*mm*(Smw[1]*Gm(states[0],mit3)+Smw[2]*Gm(states[1],mit3)+Smw[3]*Gm(states[2],mit3)+Smw[0]*Gm(states[3],mit3));
	}
	auto t2 = high_resolution_clock::now();
	fsec dt1 = t1-t0;
	fsec dt2 = t2-t1;
	t_out += dt2;
	if (debug)
	  fmt::print(cout, "{:d} {:.4f} {:.4f} {:.4f} {:.4f} {:d}\n", it3, dt1.count(), dt2.count(), t_ins.count(), t_out.count(), Nt);
      }
      if (mpi.size>1){
	mpi.AllReduceSum_simple_array(Smc);
	if (cut_references)
	  for (auto& cut:Smc_cut) mpi.AllReduceSum_simple_array(cut);
	mpi.AllReduceSum_rawdata(&Phi,1);
      }
    }
    Phi /= Q;
    fsec t_tot=t_ins+t_out;
    fmt::print(*mpi.out, "  # t_OCA_PPSigma={:.4f}s t_inside={:.4f} t_outside={:.4f} \n", t_tot.count(), t_ins.count(), t_out.count());
  
    /////////////////////////????????????????????
    double Phi_tot=0;
    // For faster MPI communication we transpose G00t and F00t
    Array<double,2> G00tc(Nt,cix.Nbath_unique_unique); G00tc=0;
    Array<double,2> F00tc(Nt,cix.Nbath_unique_unique); F00tc=0;
    Array<double,2> F01tc(Nt,cix.Nbath_unique_unique); F01tc=0;
    if (cmpPhysG){
      t_ins = zero_t;
      t_out = zero_t;
      for (int ij=0; ij<assign[mpi.rank].size(); ++ij){
	if (assign[mpi.rank][ij].n!=1) continue;
	// here we use the order : 0 < t1 < t2 < t3 < beta
	int it2 = assign[mpi.rank][ij].id % Nt;
	//for (int it2=0; it2<Nt; ++it2){ // outside loop is not with respect to t2, not t3 as above
	auto t0 = high_resolution_clock::now();    
	// Precomputing the inside loop using unique_states_inside2
	for (int i=0; i<unique_states_inside[2].size(); ++i){
	  TinyVector<int,4> st = unique_states_inside[2][i];
	  for (int it3=it2; it3<Nt; it3++){
	    if (it2==0 && it3==0) continue;
	    int istart = it2!=0 ? 0 : 1;
	    double dsum=0;
	    for (int it1=istart; it1<=it2; it1++){
	      double dx=dh*intg_dh(it1,it2);
	      dsum += dx* Gm(st[0],it1) * Gm(st[1],it2-it1) * Delta(st[2],Nt-1-(it3-it1));
	    }
	    dSmwi[2][i](it3) = dsum;
	  }
	}
	// Precomputing the inside loop using unique_states_inside3
	for (int i=0; i<unique_states_inside[3].size(); ++i){
	  TinyVector<int,4> st = unique_states_inside[3][i];
	  for (int it3=it2; it3<Nt; it3++){
	    if (it2==0 && it3==0) continue;
	    int istart = it2!=0 ? 0 : 1;
	    double dsum=0;
	    for (int it1=istart; it1<=it2; it1++){
	      double dx=dh*intg_dh(it1,it2);
	      dsum += dx* Gm(st[0],it1) * Gm(st[1],it2-it1) * Delta(st[2],it3-it1);
	    }
	    dSmwi[3][i](it3) = dsum;
	  }
	}
	auto t1 = high_resolution_clock::now();    
	t_ins += t1-t0;
	double dz=dh*intg_dh(it2,Nt-1);
	for (const auto& path : compact_paths) {
	  TinyVector<int,6> st = unpack6(path.key);
	  double mm = -path.weights[0];   // minus comes from the crossing of hybridization lines in all three terms below.
	  double nds1  = -path.weights[1]; // this is for computing vertex F(omega) and 
	  double nds1_ = -path.weights[3]; // this is for computing vertex F(omega) and 
	  double nds0  = -path.weights[2]; // minus comes from the crossing of hybridization lines
	  double nds0_ = -path.weights[4]; // minus comes from the crossing of hybridization lines
	  TinyVector<int,4> states(st[5],st[4],st[3],st[2]);
	  int s0=st[1], s1=st[0];
	  int ind2=path.inner_index[2], ind3=path.inner_index[3];
	  TinyVector<double,2> tg00; tg00=0;
	  // Below we have the integration over t3 and over all possible contributions, using precomputed parts
	  // to avoid integration over t1.
	  //							                                                              
	  //          alw[0](j0,0) =  b3   a0 = alw[0](j0,1)	                             
	  //                              t2			                            
	  //                              ^                                                   dSmwi[3] 
	  //                              | s1                                                b3  t1  a0
	  //                           --->---        		                             --->---        		  
	  //              dSmwi[2]   /         \  states[0],ibs[0]                dSmwi[3]   /  s1|    \    dSmwi[3]		  
	  //                     a3 |    s0     | b0 = alw[1](j1,0)       alw[3](j3,1) = a3 |     |     | b0 = alw[1](j1,0)     
	  //         dSmwi[2]     t1|----<------|t3   op(1)                      s0, t=0 ->-|nds0 |     |->- t2   s0
	  //                     b2 |           | a1 = alw[1](j1,1)       alw[3](j3,0) = b2 |     v     | a1 = alw[1](j1,1)	  
	  //              dSmwi[2]   \  nds1   /  states[1],ibs[1]          states[2],ibs[2] \  s1|    /  states[1],ibs[1]	  
	  //                           ---<---			                             ---<---
	  //                              ^                                                        t3
	  //                              | s1                                   alw[2](j2,1) = a2    b1 = alw[2](j2,0)           
	  //                             t=0			                 
	  //           alw[2](j2,1) = a2    b1 = alw[2](j2,0)	             
	  //
	  if (mm!=0){
	    for (int it3=it2; it3<Nt; it3++){
	      if (it2==0 && it3==0) continue;
	      double dy=dh*intg_dh(it3-it2,Nt-1-it2);
	      int m_t3 = Nt-1-it3;
	      tg00[1] += dy* dSmwi[2][ind2](it3) * Gm(states[0],it3-it2) * Gm(states[1],m_t3);
	      tg00[0] += dy* dSmwi[3][ind3](it3) * Gm(states[1],it3-it2) * Gm(states[2],m_t3);
	    }
	    if (mm!=0){
	      TinyVector<double,2> gc; gc=tg00*mm;
	      G00tc(it2,s1) += tg00[1]*mm;
	      G00tc(it2,s0) += tg00[0]*mm;
	      int m_t2 = Nt-1-it2;
	      Phi_tot += gc[1]*0.5*dz*Delta(s1,m_t2); // 1/2 because each diagram has two contributions to g00, hence each term is counted twice.
	      Phi_tot += gc[0]*0.5*dz*Delta(s0,m_t2); // 1/2 because each diagram has two contributions to g00, hence each term is counted twice.
	    }
	    F00tc(it2,s1) += tg00[1]*nds1;
	    F00tc(it2,s0) += tg00[0]*nds0;
	    F01tc(it2,s1) += tg00[1]*nds1_;
	    F01tc(it2,s0) += tg00[0]*nds0_;
	  }else{
	    if (nds1!=0 or nds1_!=0){
	      if (ind2>=0){
		for (int it3=it2; it3<Nt; it3++){
		  if (it2==0 && it3==0) continue;
		  double dy=dh*intg_dh(it3-it2,Nt-1-it2);
		  int m_t3 = Nt-1-it3;
		  tg00[1] += dy* dSmwi[2][ind2](it3) * Gm(states[0],it3-it2) * Gm(states[1],m_t3);
		}
		F00tc(it2,s1) += tg00[1]*nds1;
		F01tc(it2,s1) += tg00[1]*nds1_;
	      }else{
		cout<<"ERROR nds1="<<nds1<<" st="<<st<<" key="<<path.key<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
	      }
	    }
	    if (nds0!=0 or nds0_!=0){
	      /// WARNING!!! HERE WAS THE ISSUE. Is else here allowed? 
	      if (ind3>=0){
		for (int it3=it2; it3<Nt; it3++){
		  if (it2==0 && it3==0) continue;
		  double dy=dh*intg_dh(it3-it2,Nt-1-it2);
		  int m_t3 = Nt-1-it3;
		  tg00[0] += dy* dSmwi[3][ind3](it3) * Gm(states[1],it3-it2) * Gm(states[2],m_t3);
		}
		F00tc(it2,s0) += tg00[0]*nds0;
		F01tc(it2,s0) += tg00[0]*nds0_;
	      }else{
		cout<<"ERROR nds0="<<nds0<<" st="<<st<<" key="<<path.key<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
	      }
	    }
	  }
	}
	auto t2 = high_resolution_clock::now();
	fsec dt1=t1-t0;
	fsec dt2=t2-t1;
	t_out += dt2;
	if (debug)
	  fmt::print(cout, "{:d}  {:.4f} {:.4f} {:.4f} {:.4f} {:d}\n", it2, dt1.count(), dt2.count(), t_ins.count(), t_out.count(), Nt);
      }
      if (mpi.size>1){
	mpi.AllReduceSum_simple_array(G00tc);
	mpi.AllReduceSum_simple_array(F00tc);
	mpi.AllReduceSum_simple_array(F01tc);
	mpi.AllReduceSum_rawdata(&Phi_tot,1);
      }
      Phi_tot /= Q;
      (*mpi.out)<<"  # Phi_OCA="<<Phi<<" ; "<<Phi_tot<<endl;

      // Next we normalize the results by Q and degeneracy
      Array<double,1> deg_bth(cix.Nbath_unique_unique); deg_bth=0; // degeneracy of each bath
      for (int iorb=0; iorb<cix.Nbath; ++iorb) deg_bth(abs(cix.bath_sigind0(iorb)))+= 1;
      for (int s=0; s<cix.Nbath_unique_unique; ++s){
	double nrm = 1.0/(deg_bth(s)*Q);
	G00tc(Range::all(),s) *= nrm;
	F00tc(Range::all(),s) *= nrm;
	F01tc(Range::all(),s) *= nrm;
      }
      // Next we spline the result, which is still on uniform mesh
      Array<double,1> tau_uniform(Nt); // this is the uniform mesh we were calculating on.
      for (int it=0; it<Nt; ++it) tau_uniform(it) = beta*it/(Nt-1.0);
      Spline1DMany fG00(Nt,3*cix.Nbath_unique_unique);
      for (int s=0; s<cix.Nbath_unique_unique; ++s){
	for (int it=0; it<Nt; ++it){
	  fG00(it,3*s  ) = G00tc(it,s);
	  fG00(it,3*s+1) = F00tc(it,s);
	  fG00(it,3*s+2) = F01tc(it,s);
	}
      }
      fG00.splineIt(tau_uniform);
      // Finally we write it back in terms of larger cix.Nbath_unique array, and also in non-uniform tau mesh
      Array<int,1> sigind(cix.Nbath_unique);
      for (int iorb=0; iorb<cix.Nbath; ++iorb) sigind(cix.bath_sigind(iorb)) = cix.bath_sigind0(iorb);
      Array<bool,1> bath_SC(cix.Nbath_unique); bath_SC=false;
      for (int iorb=0; iorb<cix.Nbath; ++iorb) bath_SC(cix.bath_sigind(iorb)) = cix.SC_bath(iorb);
      
      G00t.resize(cix.Nbath_unique,tau.extent(0)); G00t=0;
      F00t.resize(cix.Nbath_unique,tau.extent(0)); F00t=0;
      F01t.resize(cix.Nbath_unique,tau.extent(0)); F01t=0;
      int pos=0;
      for (int it=0; it<tau.extent(0); it++){
	intpar ips = InterpLeft(tau(it), pos,tau_uniform);
	for (int ib=0; ib<cix.Nbath_unique; ib++){
	  int s = abs(sigind(ib));
	  double s_sgn = sigind(ib) < 0 ? -1 : 1;
	  double extra_2 = bath_SC(ib) ? 2.0 : 1.0;
	  G00t(ib,it) = fG00(ips,3*s  ) * s_sgn;
	  F00t(ib,it) = fG00(ips,3*s+1) * (extra_2*s_sgn);
	  F01t(ib,it) = fG00(ips,3*s+2) * (extra_2*s_sgn);
	}
      }
      fsec t_tot = t_out+t_ins;
      fmt::print(*mpi.out, "  # t_OCA_Gimp={:.4f}s t_inside={:.4f} t_outside={:.4f} \n", t_tot.count(), t_ins.count(), t_out.count());
      /*
      ofstream dout1("debug3.dat");
      for (int it=0; it<Nt; ++it){
	dout1 << (dh*it) << " ";
	for (int s=0; s<cix.Nbath_unique_unique; ++s) dout1 << G00tc(it,s) << " " << F00tc(it,s) << " ";
	dout1 << endl;
      }
      */
    }
  }

  {
    // Pseudo-self energies are in Smc(Nt,Npseudo_unique). We interpolate them on origonal mesh now.
    Array<double,1> tau_uniform(Nt); // this is the uniform mesh we were calculating on.
    for (int it=0; it<Nt; ++it) tau_uniform(it) = beta*it/(Nt-1.0);
    
    Spline1DMany fSmc(Nt,cix.Npseudo_unique);
    for (int k=0; k<cix.Npseudo_unique; ++k)
      for (int it=0; it<Nt; ++it) fSmc(it,k) = Smc(it,k);
    fSmc.splineIt(tau_uniform);
    
    Sm.resize(cix.Npseudo_unique,tau.extent(0));
    int pos=0;
    for (int it=0; it<tau.extent(0); it++){
      intpar ips = InterpLeft(tau(it), pos,tau_uniform);
      for (int k=0; k<cix.Npseudo_unique; k++) Sm(k,it) = fSmc(ips,k);
    }

    if (cut_references){
      for (int cut=0; cut<4; ++cut){
        Spline1DMany fcut(Nt,cix.Npseudo_unique);
        for (int k=0; k<cix.Npseudo_unique; ++k)
          for (int it=0; it<Nt; ++it) fcut(it,k)=Smc_cut[cut](it,k);
        fcut.splineIt(tau_uniform);
        (*cut_references)[cut].resize(cix.Npseudo_unique,tau.extent(0));
        pos=0;
        for (int it=0; it<tau.extent(0); ++it){
          const intpar ips=InterpLeft(tau(it),pos,tau_uniform);
          for (int k=0; k<cix.Npseudo_unique; ++k) (*cut_references)[cut](k,it)=fcut(ips,k);
        }
      }
    }
  }
  
  return Phi;
}
#endif
