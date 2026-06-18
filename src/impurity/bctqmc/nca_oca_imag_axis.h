#ifndef _nca_oca_imag_
#define _nca_oca_imag_

#include <blitz/array.h>
#include "read_cix.h"
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

/*
vector<int> ImportanceOfPseudo(const Spline1DMany& Gms, const Cix& cix, std::array<double,2> cut={0.1,0.5}, bool debug=false)
{
  if (Gms.size()!=cix.Npseudo_unique){ cerr << "ERROR in ImportanceOfPseudo: Gms.size="<<Gms.size()<<" should be Npseudo_unique="<<cix.Npseudo_unique<<endl;}
  
  Array<double,1> Gwidth(cix.Npseudo_unique);
  int middle = Gms.x_.extent(0)/2;
  // just taking |G(beta/2)| as the marker of PP importance. It's simple and reasonable, and is simply related to G(omega=0)
  for (int ipp=0; ipp<cix.Npseudo_unique; ++ipp) Gwidth(ipp) = fabs(Gms(middle,ipp));
  
  std::vector<int> idx(cix.Npseudo_unique);
  std::iota(idx.begin(), idx.end(), 0); // 0..N-1
  std::sort(idx.begin(), idx.end(), [&](int a, int b){ return Gwidth(a) > Gwidth(b); });

  vector<int> G_importance(cix.Npseudo_unique,0);
  for (int k = 0; k < (int)idx.size(); ++k) {
    int ipp = idx[k];
    if (k>cut[0]*cix.Npseudo_unique) G_importance[ipp]=1;
    if (k>cut[1]*cix.Npseudo_unique) G_importance[ipp]=2;
  }
  if (debug){
    for (int k = 0; k < (int)idx.size(); ++k) {
      int ipp = idx[k];
      fmt::print(std::cout, "{:<3d} {:<4d} {:3d} {:10.3g}\n", k, ipp+1, G_importance[ipp], Gwidth(ipp));
    }
  }
  return G_importance;
}
*/

double I_OCA(Array<double,2>& Sm, Array<double,2>& G00t, Array<double,2>& F00t, Array<double,2>& F01t,
	     bool cmpPhysG, double Q, double beta, const Array<double,1>& tau,
	     const Spline1DMany& Gms, const Spline1DMany& Deltas,
	     Cix& cix, json& pars, const vector<int>& G_importance, int OCA_cut_sum, MPIcls& mpi, double r_OCA_mesh=1)
{
  ////////////////////////////
  // input :
  //          Gms(Ntau,Npseudo)  -- input pseudoparticle splines propagators
  //          Deltas(Ntau,Ns)    -- input hybridization with splines
  //          tau(Ntau)          -- non-equidistant tau mesh
  //
  // output:
  //          Sm(Npseudo,Ntau)  -- pseudoparticle self-energy
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
  bool debug = false;
  const double coef[]={17/48.,59/48.,43/48.,49./48.,1.0};// O(dh^4)
  const int Nd=4;
  Diag2Matrix diag2matrix(6, cix);
  Array<int,2> diag(2,2);     // OCA diagram
  diag = 0,2,1,3;             // OCA diagram
  Array<int,1> dirspin(2);    // dirspin=2*iorb1,2*iorb0

  //int OCA_cut_sum = pars["OCA_cut_sum"].get<int>();
  //std::array<double,2> OCA_cut({pars["OCA_cut1"].get<double>(),pars["OCA_cut2"].get<double>()});

  if (fabs(beta-tau(tau.size()-1))>1e-7){ cerr<<"ERROR in I_OCA: last point in tau="<< tau(tau.size()-1)<<" should be beta="<<beta<<endl;}
  // we need uniform mesh for OCA
  int Nt = tau.size()*r_OCA_mesh;
  double dh = beta/(Nt-1.0);
  // functions on the uniform mesh, evaluated from splines
  Array<double,2> Gm(Gms.size(),Nt), Delta(Deltas.size(),Nt);
  int pos=0;
  for (int it=0; it<Nt; ++it){
    intpar ips = InterpLeft(dh*it, pos,tau);
    for (int a=0; a<Gms.size(); ++a) Gm(a,it) = Gms(ips,a);
    for (int b=0; b<Deltas.size(); ++b) Delta(b,it) = Deltas(ips,b);
  }
  /*
  ofstream dout1("debug1.dat");
  ofstream dout2("debug2.dat");
  for (int it=0; it<Nt; ++it){
    dout1 << (dh*it) << " ";
    for (int a=0; a<cix.Npseudo_unique; ++a) dout1 << Gm(a,it) << " ";
    for (int b=0; b<cix.Nbath_unique_unique; ++b) dout1 << Delta(b,it) << " ";
    dout1 << endl;
  }
  for (int it=0; it<tau.size(); ++it){
    dout2 << tau(it) << " ";
    for (int a=0; a<cix.Npseudo_unique; ++a) dout2 << Gms(it,a) << " ";
    for (int b=0; b<cix.Nbath_unique_unique; ++b) dout2 << Deltas(it,b) << " ";
    dout2 << endl;
  }
  */
  for (int a=0; a<Gms.size(); ++a)
    if (fabs(Gm(a,Nt-1)-Gms(tau.size()-1,a))>1e-8){ cerr<<"ERROR in I_OCA last point Gm("<<a<<","<<Nt-1<<")="<<Gm(a,Nt-1)<<" should be equal to Gms="<<Gms(tau.size()-1,a)<<endl;}
  //vector<int> G_importance = ImportanceOfPseudo(Gms, cix, OCA_cut, false);//{0.1,0.5});
  
  //cout<<"G_importance="<< G_importance << endl;
  
  // Data which store the components of the diagram, so that integration looks like a problem of a single orbital.
  int states_unique_size=0;
  {// This part is carbon-copy of the algorithm below, except that we only calculate the size of the hash-table states_unique.
    // You should not try to understand that, because we comment only the complete algorithm below. This is just stripped/off equivalent,
    // which does not compute anything except total size of the storage.
    unordered_map<uint64_t,int> states_unique;
    states_unique.max_load_factor(0.7f);
    states_unique.reserve(cix.Natst*cix.Natst*cix.Nbath*cix.Nbath);
    // for each unique key in states_unique map (which uniquely specifies the diagram contribution) the value of the map gives index to array mm_nds, index0,.., index3
    int mm_nds_size=0; // contains the matrix elements for Phi, and Nds matrix elements for vertex F(omega) calculation. There are two Nds values, because of two hybridizations.
    // Below we find all contributions to the OCA Feynman diagram.
    for (int iorb0=0; iorb0<cix.Nbath; ++iorb0){ // over all baths for first hybridiation
      int s0 = abs(cix.bath_sigind0(iorb0));
      for (int iorb1=0; iorb1<cix.Nbath; ++iorb1){ // bath for the second hybridization
	int s1 = abs(cix.bath_sigind0(iorb1));
	dirspin = 2*iorb1, 2*iorb0;                // directions are spins choosen to go in positive direction
	auto& allowed  = diag2matrix.Cmp_matrix(diag,dirspin,cix); // computes all allowed Green's functions that survice Pauli principle.
	for (const auto& [ib, alw] : allowed){     // going over all allowed contributions from all blocks.
	  if (alw[0].extent(0)>0 and alw[1].extent(0)>0 and alw[2].extent(0)>0 and alw[3].extent(0)>0){
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
		    uint64_t st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // we pack the two hybridizations (s1,s0) and atomic states into unique index-vector.
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
			uint64_t st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // all atributes of diagram contribution.
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
			uint64_t st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // all atributes of diagram contribution.
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
    // Finally we know the size of the hash-table. Below we will go through the same algorithm once more and actually calculate matrix elements.
    states_unique_size = states_unique.size();
  }
  (*mpi.out)<<"states_unique.size="<<states_unique_size<<endl;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //map<TinyVector<int,6>,int,Class1Compare> states_unique; // every unique term is stored in this map with indices {s1,s0,states[3],states[2],states[1], states[0]}
  unordered_map<uint64_t,int> states_unique;
  states_unique.reserve(states_unique_size);
  //states_unique.max_load_factor(0.7f);
  //states_unique.reserve(cix.Natst*cix.Natst*cix.Nbath*cix.Nbath);
  // for each unique key in states_unique map (which uniquely specifies the diagram contribution) the value of the map gives index to array mm_nds, index0,.., index3
  int mm_nds_size=0;
  vector<TinyVector<double,5>> mm_nds(states_unique_size); // contains the matrix elements for Phi, and Nds matrix elements for vertex F(omega) calculation. There are two Nds values, because of two hybridizations.
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
  vector<vector<TinyVector<int,4>>>  unique_states_inside(4);//0, unique_states_inside1, unique_states_inside2, unique_states_inside3; 
  //vector<int> index0(states_unique_size,-1), index1(states_unique_size,-1), index2(states_unique_size,-1), index3(states_unique_size,-1);
  vector<array<int,4>> indexi(states_unique_size,array<int,4>{-1,-1,-1,-1});
  {// index_of_inside can be temporary
    unordered_map<TinyVector<int,4>,int,TV4Hash,TV4Equal> index_of_inside0, index_of_inside1, index_of_inside2, index_of_inside3; // value → index
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
		    uint64_t st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // we pack the two hybridizations (s1,s0) and atomic states into unique index-vector.
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
			uint64_t st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // all atributes of diagram contribution.
			auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // Did we see that contribution before?
			int _iind_ = it->second; // in mm_nds and indexi, it is stored in _iind_
			if (inserted){// not seen before, hence insert into mm_nds, index0, index1, index2, index3
			  if (destroys_Cooper_pair_or_normal_down2){
			    mm_nds[mm_nds_size]=TinyVector<double,5>(0,0,0,0,nds2);
			  }else
			    mm_nds[mm_nds_size]=TinyVector<double,5>(0,0,nds2,0,0);
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
			uint64_t st = pack6(s1, s0, states[3], states[2], states[1], states[0]); // all atributes of diagram contribution.
			auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // Did we see that contribution before?
			int _iind_ = it->second; // in mm_nds and indexi, it is stored in _iind_
			if (inserted){ // not seen before, hence insert into mm_nds, index0, index1, index2, index3
			  if (destroys_Cooper_pair_or_normal_down1)
			    mm_nds[mm_nds_size]=TinyVector<double,5>(0,0,0,nds1,0);
			  else
			    mm_nds[mm_nds_size]=TinyVector<double,5>(0,nds1,0,0,0);
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
  }

  (*mpi.out)<<"Final mm_nds_size="<<mm_nds_size<<endl;
  {
    int ijk=0;
    for (const auto& [_st_, iind] : states_unique){
      TinyVector<int,6> st = unpack6(_st_);
      double mm = mm_nds[iind][0];
      double nds1 = abs(mm_nds[iind][1])+abs(mm_nds[iind][3]);
      double nds2 = abs(mm_nds[iind][2])+abs(mm_nds[iind][4]);
      int s1 = st[0];
      int s0 = st[1];
      TinyVector<int,Nd> states(st[5],st[4],st[3],st[2]);
      int index_of_0 = indexi[iind][0];
      int index_of_1 = indexi[iind][1];
      int index_of_2 = indexi[iind][2];
      int index_of_3 = indexi[iind][3];
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
	cout << "ERROR id[2]="<<index_of_2<<" at st="<<st<<" iind="<<iind<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
      }
      if (index_of_3<0){
	cout << "ERROR id[3]="<<index_of_3<<" at st="<<st<<" iind="<<iind<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
      }
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
    for (const auto& [_st_, iind] : states_unique){
      TinyVector<int,6> st = unpack6(_st_);
      double mm = mm_nds[iind][0];
      double nds1 = mm_nds[iind][1];//abs(mm_nds[iind][1])+abs(mm_nds[iind][3]);
      double nds2 = mm_nds[iind][2];//abs(mm_nds[iind][2])+abs(mm_nds[iind][4]);
      if (fabs(mm_nds[iind][3])>1e-7 or fabs(mm_nds[iind][4])>1e-7) cout<<"ERROR expecting nds3="<<mm_nds[iind][3]<<" and nds4="<<mm_nds[iind][4]<<" to vanish"<<endl;
      int s1 = st[0];
      int s0 = st[1];
      TinyVector<int,Nd> states(st[5],st[4],st[3],st[2]);
      int index_of_0 = indexi[iind][0];
      int index_of_1 = indexi[iind][1];
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
      int index_of_2 = indexi[iind][2];
      TinyVector<int,4> sti2 = unique_states_inside[2][index_of_2];
      if (sti2[0]!=states[2] or sti2[1]!=states[3] or sti2[2]!=s0){
	cout<<"WARNING2 sti2="<<sti2<<" while states[2],states[3],s0="<<states[2]<<","<<states[3]<<","<<s0<<" and st="<<st<<endl;
      }
      int index_of_3 = indexi[iind][3];
      TinyVector<int,4> sti3 = unique_states_inside[3][index_of_3];
      if (sti3[0]!=states[3] or sti3[1]!=states[0] or sti3[2]!=s1){
	cout<<"WARNING3 sti3="<<sti3<<" while states[3],states[0],s1="<<states[3]<<","<<states[0]<<","<<s1<<" and st="<<st<<endl;
      }
      fmt::print(cout," {:6d} ii0={:<3d},ii1={:<3d},ii2={:<3d},ii3={:<3d},s0={:d},s1={:d}: mm={:2g} nds1={:3g},{:3g} nds0={:3g},{:3g}  index0={:4d} index1={:4d} index2={:4d} index3={:4d} ",
		 ijk,st[5]+1,st[4]+1,st[3]+1,st[2]+1,st[1],st[0],
		 mm, mm_nds[iind][1],mm_nds[iind][3],mm_nds[iind][2],mm_nds[iind][4],
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
      
	for (auto const& [_st_, iind] : states_unique) {
	  TinyVector<int,6> st = unpack6(_st_);
	  double mm = -mm_nds[iind][0]; // minus comes from the crossing of hybridization lines
	  if (mm==0) continue;
	  TinyVector<int,4> states(st[5],st[4],st[3],st[2]);
	  int s0=st[1], s1=st[0];
	  int ind0=indexi[iind][0], ind1=indexi[iind][1], ind2=indexi[iind][2], ind3=indexi[iind][3];
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
	  Smc(it3,states[3]) += Smw[0]*mm/cix.deg_pp(states[3]);
	  Smc(it3,states[0]) += Smw[1]*mm/cix.deg_pp(states[0]);
	  Smc(it3,states[1]) += Smw[2]*mm/cix.deg_pp(states[1]);
	  Smc(it3,states[2]) += Smw[3]*mm/cix.deg_pp(states[2]);
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
	for (auto const& [_st_, iind] : states_unique) {
	  TinyVector<int,6> st = unpack6(_st_);
	  double mm = -mm_nds[iind][0];   // minus comes from the crossing of hybridization lines in all three terms below.
	  double nds1  = -mm_nds[iind][1]; // this is for computing vertex F(omega) and 
	  double nds1_ = -mm_nds[iind][3]; // this is for computing vertex F(omega) and 
	  double nds0  = -mm_nds[iind][2]; // minus comes from the crossing of hybridization lines
	  double nds0_ = -mm_nds[iind][4]; // minus comes from the crossing of hybridization lines
	  TinyVector<int,4> states(st[5],st[4],st[3],st[2]);
	  int s0=st[1], s1=st[0];
	  int ind2=indexi[iind][2], ind3=indexi[iind][3];
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
		cout<<"ERROR nds1="<<nds1<<" st="<<st<<" iind="<<iind<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
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
		cout<<"ERROR nds0="<<nds0<<" st="<<st<<" iind="<<iind<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
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
  }
  
  return Phi;
}
#endif
