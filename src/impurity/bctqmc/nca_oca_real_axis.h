#include <cmath>
#include "fmt/core.h"

#include <unordered_set>

double R_NCA(Array<double,2>& Smt,  Array<complex<double>,2>& Gc, Array<complex<double>,2>& Fc1, Array<complex<double>,2>& Fc2, 
	     bool cmpPhysG, double Q, const Array<double,1>& xn, const Array<complex<double>,2>& Gm, const Array<double,1>& few,
	     const Array<double,2> PhiR, const Array<double,1>& om, const Array<double,2>& Deltai, const Array<double,2>& mDelta,
	     const Array<double,1>& feo, Cix& cix, double T, MPIcls& mpi)
{
  const bool debug=false;
  const int Nd=2;
  //cout<<"Inside NCA xn.size="<< xn.size() << " " << Gm.extent(1) << endl;
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
      //
      bool destroys_Cooper_pair = cix.Fi_indx(2*iorb)%2==1 and cix.Fi_indx(2*iorb+1)%2==1;
      bool normal_down = cix.SC_bath(iorb) and ((cix.Fi_indx(2*iorb)+cix.Fi_indx(2*iorb+1))%2==1 and cix.bath_spin(iorb)<0);
      //
      int jb = cix.Fbc(2*iorb,ib);              // if we apply F^+_{iorb}|ib> = |jb> . Here |jb> is block, not a state.
      if (jb<0) continue;                       // Pauli principle sometimes makes this vanish
      int Nj = cix.Block[jb].size();            // size of the |jb> block
      NjmResult nx = cix.Njm_ind(ib,2*iorb);    // always replace the end of hybridization (psi^+ in normal case) with the njm operator
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
	      njm(n1,n2) = (*nx.njm)(n2,n1);
	}
      }
      // we first apply psi_{iorb}|j1>
      for (int _j1_=0; _j1_<Nj; ++_j1_){ // over matrix corresponding to block jb
	int j1 = cix.Block[jb](_j1_);
	int i1 = cix.Fi(2*iorb+1, j1);              // |i1> = psi^(1)_{iorb}|j1>
	if (i1==0) continue;
	int _i1_ = cix.inBlock(abs(i1));
	for (int _i2_=0; _i2_<Ni; ++_i2_){          // next psi^+_{iorb}|i2>
	  int i2 = cix.Block[ib](_i2_);
	  int j2 = cix.Fi(2*iorb, i2);                // |j2> = psi^(0)_{iorb}|i2>
	  int _ii_ = cix.Sgm_symmetry[ib](_i1_,_i2_); // symmetry of propagator G(i1,i2)==G(_ii_). Here _ii_ starts with 1 (if nonzero) and can be negative (when propagator is equal but opposite to some other propagator).
	  if (_ii_==0) continue; // this G(_ii_) vanishes.
	  int ii = abs(_ii_)-1;                       // actual place where G(ii) is found.
	  if (j2!=0){
	    // this G(_ii_) is finite and psi^+_{iorb}|i2>=|j2> is finite
	    //cout<<"ib="<<ib<<" iorb="<<iorb<<" (i1,i2)=("<<i1<<","<<i2<<") ii="<<_ii_<<" (j1,j2)=("<<j1<<","<<j2<<")"<<endl;
	    int jb2 = cix.iBlck(abs(j1))-1;           // jb is block in which j1 and j2 are located. Could be obtained by cix.Fbc(2*iorb,ib+1)
	    if (jb2!=jb) cout<<"ERROR jb="<<jb<<" jb2="<<jb2<<endl;
	    int _j2_ = cix.inBlock(abs(j2));          // j2 is also counted from the empty state, while _j2_ is index within the block jb
	    int _jj_ = cix.Sgm_symmetry[jb](_j2_,_j1_); // G(j1,j2)=G(_jj_)
	    if (_jj_!=0){                               // G(_jj_) can sometimes vanish
	      int jj = abs(_jj_)-1;                     // index for G(_jj_) == G[jj]
	      TinyVector<int,3> st(s,jj,ii);            // store ii,jj,s in such a way that it is sorted with respect to ii, and after that jj, and finally s.
	      auto [it, inserted] = states_unique.try_emplace(st, mm_nds.size()); // try to insert st as a new key and mm_nds.size() as the value.
	                                                // If key is not yet present, inserted=true and "it" points to its current place. If key is present, "it" points to its existing place.
	      int iind = it->second;                    // The key st has value iind, either just inserted, or found from before
	      if (inserted) mm_nds.push_back(TinyVector<double,3>(0, 0, 0)); // if it was just inserted, deque needs to have this entry set to zero.
	      double mm = sign(static_cast<double>(i1)*static_cast<double>(j2)*static_cast<double>(_ii_)*static_cast<double>(_jj_))*s_sgn; // actual matrix element for NCA diagram.
	      // Here sign(i1)*sign(j2) is matrix element due to F^+ and F operators. The sign of _ii_ and _jj_ comes from the fact that some propagators have opposite signs as the value in G(_ii_) or G(_jj_).
	      //if (debug){
	      //if (ii >= cix.first_anomalous-1 or jj >= cix.first_anomalous-1)
	      //  cout<<"1) ib="<<ib<<" iorb="<<iorb<<" (i1,i2)=("<<i1<<","<<i2<<") (j1,j2)=("<<j1<<","<<j2<<") _ii_="<<_ii_<<" _jj_="<<_jj_<<" st="<<st<<" mm="<<mm<<endl;
	      //}
	      double nds = 0;
	      if (nx.njm_c==1 and cmpPhysG){  // if njm is proportional to F
		nds=mm*njmr; // mm contains all other terms, we just add njmr, which is the ratio.
	      }
	      if (destroys_Cooper_pair or normal_down)
		mm_nds[iind] += TinyVector<double,3>(mm, 0, nds);
	      else
		mm_nds[iind] += TinyVector<double,3>(mm, nds, 0);
	    }
	  }
	  if (nx.njm_c==2 and cmpPhysG){ // if njmc is not proportional to psi^+, but completely different matrix, we need to go over all matrix elements.
	    // We are inserting Njm instead of psi^+_{iorb}|i2>=|j2> operator. But now |i2> and |j2> have to be unrestricted, not connected by psi^+_{iorb} but Njm[i2,j2]
	    // Notice that we keep F^+_{iorb}|i1>=|j1> on the other end.
	    for (int _j2_=0; _j2_<Nj; ++_j2_){// _j2_ is now arbitrary index in block jb (size Nj), and not given by F^+_{iorb}|i2>!=|j2>
	      bool njm_finite = abs(njm(_i2_,_j2_))>1e-10; // is njm finite?
	      int _jj_ = cix.Sgm_symmetry[jb](_j2_,_j1_);  // for this general j2, what is G(j2,j1)=G[jj] and is it finite?
	      if (_jj_!=0 and njm_finite){
		double nds = njm(_i2_,_j2_)*sign(_ii_*_jj_)*sign(i1)*s_sgn;      // we here use 0.5*njm(i2,j2) because we insert the operator with magnitude 0.5 at the constructore and 0.5 at the anhilator.
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
  int first_bath_anomalous=cix.Nbath;
  for (int ib=0; ib<cix.Nbath; ++ib)
    if (cix.bath_SC(ib) and abs(cix.bath_sigind0(ib))<first_bath_anomalous)
      first_bath_anomalous = abs(cix.bath_sigind0(ib));
  //cout << "first_bath_anomalous="<< first_bath_anomalous << endl;
  if (debug){
    if (mpi.rank==mpi.master){
      int ijk=0;
      for (const auto& [st, iind] : states_unique){
	double mm   = mm_nds[iind][0];
	double nds1 = mm_nds[iind][1];
	double nds2 = mm_nds[iind][2];
	int s = st[0];
	TinyVector<int,Nd> states(st[2],st[1]);
	
	double maxv0 = max(imag(Gm(states[0],Range::all())));
	double maxv1 = max(imag(Gm(states[1],Range::all())));

	if ((states[0] >= cix.first_anomalous-1 or states[1] >= cix.first_anomalous-1) and s>=first_bath_anomalous){
	  int opposite = states_unique[TinyVector<int,3>(st[0],st[2],st[1])];
	  double mm_opp = mm_nds[opposite][0];
	  double nds1_opp = mm_nds[opposite][1];
	  double nds2_opp = mm_nds[opposite][2];

	  if (fabs(mm_opp+mm)>1e-6) {fmt::print(cout,"ERROR: mm={} mm_={}\n",mm,mm_opp); }
	  if (fabs(nds1_opp+nds2_opp+nds1+nds2)>1e-6) {fmt::print(cout,"ERROR: nds1={} nds2={} nds1_opp={} nds2_opp={}\n",nds1,nds2,nds1_opp,nds2_opp); }
	  //if (not ( v0_gt_v1 ^ v1_gt_v0) ) {fmt::print(cout,"ERROR: v0_gt_v1={} v1_gt_v0_={}\n",v0_gt_v1,v1_gt_v0); }
	    
	  fmt::print(cout," {:4d} ii0={:<3d},ii1={:<3d},s={:d}: mm={:3g} mm_={:3g} nds1={:3g} nds2={:3g} nds1_opp={:3g} nds2_opp={:3g}", ijk,states[0]+1,states[1]+1,s,mm, mm_opp,nds1,nds2,nds1_opp,nds2_opp);
	  
	  //fmt::print(cout," {:4d} ii0={:<3d},ii1={:<3d},s={:d}: mm={} nds={} deg_pp0={} ", ijk,states[0]+1,states[1]+1,s,mm,nds,cix.deg_pp(states[0]) );
	  fmt::print(cout,"\n");
	}
	ijk+=1;
	
      }
    }
  }
  
  double phi1=0;
  map<double,int> mesh_eps;
  for (int iep=0; iep<xn.size(); ++iep) mesh_eps.insert({xn(iep),iep}); // all points of PP mesh
  Array<double,2> eps_w(xn.size()+om.size()); // mesh for eps-om, where om is mesh for physical spectra and eps mesh for pseudo spectra,
  Array<double,2> ff_w (xn.size()+om.size()); // contains f(eps-w)*f(-eps)/f(-w) with eps dense around zero and at om
  Array<int,2> eps_i(xn.size()+om.size());
  Array<double,2> TSmt(xn.size(),cix.Npseudo_unique);
  TSmt=0;

  fsec t_tot=zero_t, t_dt=zero_t;
  auto [iqs,iqe] = mpi.mpiSplitArray(xn.size());
  //int iqs=0, iqe=xn.size();
  for (int iom=iqs; iom<iqe; ++iom){
    //auto t0 = high_resolution_clock::now();
    double w=xn(iom);
    map<double,int> m_eps(mesh_eps);
    for (int i=0; i<om.size(); ++i) m_eps.insert({om(i)+w,-i-1}); // to PP mesh (mesh_eps) adding om(i)+w points
    int eps_size=m_eps.size();
    int ii=0;
    for (auto [eps, iep] : m_eps) { // going over all combined points for eps
      double ffx=0, f_eps=0, f_eps_w=0;
      if (iep>=0){  // point is from mesh xn
	f_eps = few(iep);
	f_eps_w = fermi((eps-w)/T);
	if (w>0) // for positive w, we compute f(eps-w)*f(-eps)/f(-w)
	  ffx = f_eps_w*(1-f_eps)/(1-few(iom));
	else     // for negative w, f(w-eps)*f(eps)/f(w)
	  ffx = (1-f_eps_w)*f_eps/few(iom);
      }else{ // point is from mesh om[:]+w
	f_eps = fermi(eps/T);
	f_eps_w = feo(-iep-1); //fermi(beta*(eps-w));
	if (w>0)
	  ffx = f_eps_w*(1-f_eps)/(1-few(iom));
	else
	  ffx = (1-f_eps_w)*f_eps/few(iom);
      }
      eps_w(ii) = eps;
      ff_w (ii) = ffx;
      eps_i(ii) = iep;
      ii+=1;
    }
    Array<double,1> Sg1(Gm.extent(0)), Sg0(Gm.extent(0)); Sg1=0, Sg0=0;
    int pos_e=0, pos_x=0, pos_y=om.size()-1;//, posx=0;
    for (int ii=0; ii<eps_size; ++ii){
      double eps=eps_w(ii);
      if (eps<xn(0) or eps>xn(xn.size()-1)) continue;
      double dh=0;
      if (ii>0 and ii<eps_size-1) 
	dh = 0.5*(eps_w(ii+1)-eps_w(ii-1));
      else{
	int N=eps_size;
	dh = (ii==0) ? 0.5*(eps_w(1)-eps_w(0)) : 0.5*(eps_w(N-1)-eps_w(N-2));
      }
      double ff_dh = ff_w(ii)*dh;  // f(eps-w)*f(-eps)/f(-w) 
      int    iep=eps_i(ii);
      intpar ip_x, ip_y, ip_e;
      if (iep>=0){
	ip_x=InterpLeft(eps-w,pos_x,om);
	ip_y=InterpRight(w-eps,pos_y,om);
      }else{
	ip_e=InterpLeft(eps,pos_e,xn);
      }
      for (const auto& [st, iind] : states_unique){
	//for (auto it=states_unique.begin(); it!=states_unique.end(); ++it){
	int s = st[0];
	int ip0 = st[2];
	int ip1 = st[1];
	double mm  = mm_nds[iind][0];
	//double nds = mm_nds[iind][1];
	//
	if (mm!=0){
	  double G1t, G0t;
	  double D1,D0;
	  if (iep>=0){
	    G1t = Gm(ip1,iep).imag();      // Gt[ip1](eps)
	    G0t = Gm(ip0,iep).imag();      // Gt[ip0](eps)
	    D1 = linInterp(Deltai,s,ip_x); // Dl[s](eps-w)
	    D0 = linInterp(Deltai,s,ip_y); // Dl[s](w-eps)
	    // 
	  }else{
	    G1t = linInterp(Gm,ip1,ip_e).imag();// Gt[ip1](eps)
	    G0t = linInterp(Gm,ip0,ip_e).imag();// G[ip0](eps)
	    D1 =  Deltai(s,-iep-1);             // Dl[s](eps-w);
	    D0 =  mDelta(s,-iep-1);             // Dl[s](-eps+w);
	  }
	  Sg1(ip0) += ff_dh*D1*G1t*mm; // f(eps-w)*f(-eps)/f(-w) Dl[s](eps-w) G[ip1](eps)
	  Sg0(ip1) += ff_dh*D0*G0t*mm; // f(eps-w)*f(-eps)/f(-w) Dl[s](w-eps) G[ip1](eps)
	  //double dd1 = ff_dh*D1*G1t*mm;
	  //double dd2 = ff_dh*D0*G0t*mm;
	  //if (std::isnan(dd1) or std::isnan(dd2)){
	  //  cout<<"dd1="<<dd1<<" dd2="<<dd2<<" ff_dh="<<ff_dh<<" D1="<<D1<<" D0="<<D0<<" G1t="<<G1t<<" G0t="<<G0t<<" mm="<<mm<<" iep="<<iep<<" ip_y="<<ip_y<<endl;
	  //}
	}
      }
    }
    TSmt(iom,Range::all()) += Sg1(Range::all());///deg_pp(Range::all());
    TSmt(iom,Range::all()) += Sg0(Range::all());///deg_pp(Range::all());
    //auto t1 = high_resolution_clock::now();
    //t_dt = t1-t0;
    //t_tot += t_dt;
    //fmt::print(cout,"{} {:.6f}  {:.6f} {}\n", iom, t_dt.count(), t_tot.count(), eps_size);
  }
  //mm is taking care of degenracy, hence we do not need to divide by deg_pp!
  mpi.AllGatherNCA1(cix.Npseudo_unique,iqs,iqe,TSmt);
  
  Array<double,1> dxh(xn.size()); dxh = GiveDh(xn);
  Smt.resize(cix.Npseudo_unique,xn.size());
  for (int ii=0; ii<cix.Npseudo_unique; ++ii){
    Smt(ii,Range::all()) = TSmt(Range::all(),ii)/cix.deg_pp(ii);
    phi1 += sum(Smt(ii,Range::all())*PhiR(ii,Range::all())*dxh(Range::all()))*cix.deg_pp(ii);
  }
  phi1 /= (2*Q); // because there are 2-pp propagators, and Phi = Sigma_p0*G_p0 or Phi=Sigma_p1*G_p1 or sum/2.

  if (cmpPhysG){
    Array<complex<double>,2> TGcm(om.size(),cix.Nbath_unique_unique), TGc(om.size(),cix.Nbath_unique_unique); TGcm=0; TGc=0;
    Array<complex<double>,2> TFcm1(om.size(),cix.Nbath_unique_unique), TFc1(om.size(),cix.Nbath_unique_unique); TFcm1=0; TFc1=0;
    Array<complex<double>,2> TFcm2(om.size(),cix.Nbath_unique_unique), TFc2(om.size(),cix.Nbath_unique_unique); TFcm2=0; TFc2=0;
    Array<bool,2> maxv0_gt_maxv1(cix.Npseudo_unique,cix.Npseudo_unique); maxv0_gt_maxv1=false;
    {
      intpar i0 = Interp(0.0,xn); // which point is at zero frequency
      int di = static_cast<int>(xn.extent(0)/15); // let's take 7% of points below and above zero frequency to check how peaked is PP
      Array<double,1> maxv(cix.Npseudo_unique);
      for (int ip0=0; ip0<cix.Npseudo_unique; ++ip0){
	double val = 0, der = 0, dh; // we will check the integral and the absolute value of the integral of the derivative
	for (int iw=i0.i-di; iw<i0.i+di+1; ++iw){
	  dh = om(iw+1)-om(iw);
	  val += abs(Gm(ip0,iw).imag())*dh;
	  der += abs(Gm(ip0,iw+1).imag()-Gm(ip0,iw).imag());
	}
	maxv(ip0) = val+der; // this should be larger when function is strongly peaked around zero frequency
      }
      for (int ip0=0; ip0<cix.Npseudo_unique; ++ip0)
	for (int ip1=0; ip1<ip0; ++ip1){
	  maxv0_gt_maxv1(ip0,ip1) = (maxv(ip0)>maxv(ip1));
	  maxv0_gt_maxv1(ip1,ip0) = not maxv0_gt_maxv1(ip0,ip1);
	}
    }
    
    Array<double,2> fx_eps(xn.size()+om.size()), fx_eps_w(xn.size()+om.size());
    mesh_eps.clear();
    for (int iep=0; iep<om.size(); ++iep) mesh_eps.insert({om(iep),iep}); //starts with physical mesh
    
    auto [iws,iwe] = mpi.mpiSplitArray(om.size());
    //int iws=0, iwe=om.size();
    for (int iom=iws; iom<iwe; ++iom){
      //auto t0 = high_resolution_clock::now();
      double w=om(iom);
      map<double,int> m_eps(mesh_eps);
      for (int i=0; i<xn.size(); ++i) m_eps.insert({xn(i)+w,-i-1}); // adding PP mesh to physical om mesh -> xn(i)+w points
      int eps_size=m_eps.size();
      int ii=0;
      for (auto it=m_eps.begin(); it!=m_eps.end(); ++it){ // going over all combined points for eps
	double eps=it->first;
	int iep = it->second; // if positive, it is from mesh eps and iep is index in eps==xn(iep). If negative, eps==om(-iep-1)+w
	double f_eps=0, f_eps_w=0;
	if (iep>=0){  // point is from mesh xn
	  f_eps = feo(iep);
	  f_eps_w = fermi((eps-w)/T);
	}else{ // point is from mesh om[:]+w
	  f_eps = fermi(eps/T);
	  f_eps_w = few(-iep-1); //fermi(beta*(eps-w));
	}
	eps_w(ii) = eps;
	eps_i(ii) = iep;
	fx_eps(ii) = f_eps;
	fx_eps_w(ii) = f_eps_w;
	ii+=1;
      }
      int pos_e=0, posx=0;
      for (int ii=0; ii<eps_size; ++ii){
	double eps=eps_w(ii);
	if (eps<xn(0) or eps>xn(xn.size()-1)) continue;
	double dh=0;
	if (ii>0 and ii<eps_size-1) 
	  dh = 0.5*(eps_w(ii+1)-eps_w(ii-1));
	else{
	  int N=eps_size;
	  dh = (ii==0) ? 0.5*(eps_w(1)-eps_w(0)) : 0.5*(eps_w(N-1)-eps_w(N-2));
	}
	double f_eps = fx_eps(ii);
	double f_eps_w = fx_eps_w(ii);
	intpar ipe = InterpLeft(eps,pos_e,xn);
	intpar ipx  = InterpLeft(eps-w,posx,xn);
	int ijk=0;
	for (const auto& [st, iind] : states_unique){
	  int s = st[0];
	  int ip0 = st[2];
	  int ip1 = st[1];
	  double mm   = mm_nds[iind][0];
	  double nds1 = mm_nds[iind][1];
	  double nds2 = mm_nds[iind][2];
	  //
	  if (maxv0_gt_maxv1(ip0,ip1)){
	    // if ip0 is more peaked than ip1, we use this formula, otherwise the other
	    // The ground state should be peaked at (eps-w) because eps mesh is composed of Delta(omega) (omega<->eps) mesh and we add xn mesh peaked at omega, i.e.,
	    // so that it describes well a peak of G(eps-w).
	    complex<double> G1c = linInterp(Gm,ip1,ipe);               // Gc[ip1](eps)
	    double G1t = G1c.imag();                                   // Gt[ip1](eps)
	    G1c = complex<double>(G1c.real(),-G1c.imag()*(1-f_eps));   // Gc[ip1](eps)
	    
	    complex<double> G0c = linInterp(Gm,ip0,ipx);               // G[ip0](eps-w)
	    double G0t = G0c.imag();                                   // Gt[ip0](eps-w)
	    G0c = complex<double>(G0c.real(),-G0c.imag()*(1-f_eps_w)); // G[ip0](eps-w)
	    
	    complex<double> gd = G1c*f_eps_w*G0t - conj(f_eps*G1t*G0c);// Gc[ip1](eps)*f(eps-w)*Gt[ip0](eps-w)-(Gt[ip1](eps)*f(eps)*Gc[ip0](eps-w))^*
	    if (mm!=0)   TGc(iom,s) += gd*mm*dh;
	    if (nds1!=0) TFc1(iom,s) += gd*nds1*dh;
	    if (nds2!=0) TFc2(iom,s) += gd*nds2*dh;
	  }else{
	    // here ip1 is more peaked
	    complex<double> G0c = linInterp(Gm,ip0,ipe);               // G[ip0](eps)
	    double G0t = G0c.imag();                                   // Gt[ip0](eps)
	    G0c = complex<double>(G0c.real(),-G0c.imag()*(1-f_eps));   // G[ip0](eps)
	    
	    complex<double> G1c = linInterp(Gm,ip1,ipx);               // Gc[ip1](eps-w)
	    double G1t = G1c.imag();                                   // Gt[ip1](eps-w)
	    G1c = complex<double>(G1c.real(),-G1c.imag()*(1-f_eps_w)); // Gc[ip1](eps-w)
	    
	    complex<double> gd = -conj(G0c*f_eps_w*G1t) + f_eps*G0t*G1c;// Gc[ip1](eps-w)*f(eps)*Gt[ip0](eps)-(Gt[ip1](eps-w)*f(eps-w)*Gc[ip0](eps))^*
	    // -(Gc[ip0](eps)*f(eps-w)*Gt[ip1](eps-w))^* + Gt[ip0](eps)*f(eps)*Gc[ip1](eps-w)
	    if (mm!=0)   TGcm(iom,s)  += gd*mm*dh;
	    if (nds1!=0) TFcm1(iom,s) += gd*nds1*dh;
	    if (nds2!=0) TFcm2(iom,s) += gd*nds2*dh;
	  }
	  ijk+=1;
	}
      }
      //auto t1 = high_resolution_clock::now();
      //t_dt = t1-t0;
      //t_tot += t_dt;
      //fmt::print(cout,"{} {:.6f}  {:.6f} {}\n", iom, t_dt.count(), t_tot.count(), eps_size);
    }
    mpi.AllGatherNCA3(cix.Nbath_unique_unique,iws,iwe,TGc,TGcm,TFc1,TFcm1,TFc2,TFcm2);
    {
      if (debug){
	ofstream out("AllG.debug.dat");
	for (int iw=0; iw<om.size(); ++iw){
	  out <<om(iw)<<" ";
	  for (int ib=0; ib<TGc.extent(1); ++ib)
	    out << TGc(iw,ib).real()/(M_PI*Q)<<" "<<TGc(iw,ib).imag()/(M_PI*Q)<<" "<<TGcm(iw,ib).real()/(M_PI*Q)<<" "<<TGcm(iw,ib).imag()/(M_PI*Q)<<" ";
	  out <<endl;
	}
	ofstream out1("Anomal.F1.dat");
	for (int iw=0; iw<om.size(); ++iw){
	  out1 <<om(iw)<<" ";
	  for (int ib=first_bath_anomalous; ib<TFc1.extent(1); ++ib)
	    out1 << TFc1(iw,ib).real()/(M_PI*Q)<<" "<<TFc1(iw,ib).imag()/(M_PI*Q)<<" "<<TFcm1(iw,ib).real()/(M_PI*Q)<<" "<<TFcm1(iw,ib).imag()/(M_PI*Q)<<" ";
	  out1 <<endl;
	}
	ofstream out2("Anomal.F2.dat");
	for (int iw=0; iw<om.size(); ++iw){
	  out2 <<om(iw)<<" ";
	  for (int ib=first_bath_anomalous; ib<TFc2.extent(1); ++ib)
	    out2 << TFc2(iw,ib).real()/(M_PI*Q)<<" "<<TFc2(iw,ib).imag()/(M_PI*Q)<<" "<<TFcm2(iw,ib).real()/(M_PI*Q)<<" "<<TFcm2(iw,ib).imag()/(M_PI*Q)<<" ";
	  out2 <<endl;
	}
	ofstream out3("Normal.F1.dat");
	for (int iw=0; iw<om.size(); ++iw){
	  out3 <<om(iw)<<" ";
	  for (int ib=0; ib<first_bath_anomalous; ++ib)
	    out3 << TFc1(iw,ib).real()/(M_PI*Q)<<" "<<TFc1(iw,ib).imag()/(M_PI*Q)<<" "<<TFcm1(iw,ib).real()/(M_PI*Q)<<" "<<TFcm1(iw,ib).imag()/(M_PI*Q)<<" ";
	  out3 <<endl;
	}
	ofstream out4("Normal.F1.dat");
	for (int iw=0; iw<om.size(); ++iw){
	  out4 <<om(iw)<<" ";
	  for (int ib=0; ib<first_bath_anomalous; ++ib)
	    out4 << TFc2(iw,ib).real()/(M_PI*Q)<<" "<<TFc2(iw,ib).imag()/(M_PI*Q)<<" "<<TFcm1(iw,ib).real()/(M_PI*Q)<<" "<<TFcm2(iw,ib).imag()/(M_PI*Q)<<" ";
	  out4 <<endl;
	}
      }
      int mpos=om.size()-1;
      for (int iw=0; iw<om.size(); ++iw){
	intpar px=InterpRight(-om(iw),mpos,om);
	for (int ib=0; ib<cix.Nbath_unique_unique; ++ib){
	  TGc(iw,ib) += linInterpT(TGcm,ib,px);
	  TFc1(iw,ib) += linInterpT(TFcm1,ib,px);
	  TFc2(iw,ib) += linInterpT(TFcm2,ib,px);
	}
      }
    }

    Array<int,1> sigind(cix.Nbath_unique);
    for (int iorb=0; iorb<cix.Nbath; ++iorb) sigind(cix.bath_sigind(iorb)) = cix.bath_sigind0(iorb);
    Array<double,1> deg_bth(cix.Nbath_unique_unique); deg_bth=0; // degeneracy of each bath
    for (int iorb=0; iorb<cix.Nbath; ++iorb) deg_bth(abs(cix.bath_sigind0(iorb)))+= 1;
    Array<bool,1> bath_SC(cix.Nbath_unique); bath_SC=false;
    for (int iorb=0; iorb<cix.Nbath; ++iorb) bath_SC(cix.bath_sigind(iorb)) = cix.SC_bath(iorb);
    
    Gc.resize(cix.Nbath_unique,om.size());
    Fc1.resize(cix.Nbath_unique,om.size());
    Fc2.resize(cix.Nbath_unique,om.size());
    for (int ib=0; ib<cix.Nbath_unique; ++ib){
      int s = abs(sigind(ib));
      double s_sgn = sigind(ib) < 0 ? -1 : 1;
      Gc(ib,Range::all())  = TGc(Range::all(),s)*(s_sgn/deg_bth(s));
      double extra_2 = bath_SC(ib) ? 2.0 : 1.0;
      Fc1(ib,Range::all()) = TFc1(Range::all(),s)*(extra_2*s_sgn/deg_bth(s)); // factor of two because deg_bth combines spin-up and down. If not, we need to be more careful????
      Fc2(ib,Range::all()) = TFc2(Range::all(),s)*(extra_2*s_sgn/deg_bth(s));// factor of two because deg_bth combines spin-up and down. If not, we need to be more careful????
      //cout<<"Gc["<<ib<<"]= TGc["<<s<<"]*"<<s_sgn<<"/"<<deg_bth(s)<<endl;
    }
    Gc  *= 1./(M_PI*Q);
    Fc1 *= 1./(M_PI*Q);
    Fc2 *= 1./(M_PI*Q);

    if (debug){
      ofstream out1("G00_debug.dat");
      ofstream out2("F001_debug.dat");
      ofstream out3("F002_debug.dat");
      for (int iw=0; iw<om.size(); ++iw){
	out1 << om(iw) << " ";
	for (int ib=0; ib<Gc.extent(0); ++ib)
	  out1 << Gc(ib,iw).real() << " " << Gc(ib,iw).imag() << " ";
	out1<<endl;
	out2 << om(iw) << " ";
	for (int ib=0; ib<Fc1.extent(0); ++ib)
	  out2 << Fc1(ib,iw).real() <<" "<< Fc1(ib,iw).imag() << " ";
	out2 << endl;
	out3 << om(iw) << " ";
	for (int ib=0; ib<Fc2.extent(0); ++ib)
	  out3 << Fc2(ib,iw).real() <<" "<< Fc2(ib,iw).imag() << " ";
	out3 << endl;
      }
    }
  }
  return phi1;
}

inline double edge_width(const Array<double,1>& a, int i){
  const int N = a.extent(0);
  if (N <= 1) return 0.0;                 // single point → zero measure
  if (i > 0 && i < N-1) return 0.5*(a(i+1)-a(i-1));
  if (i == 0)           return 0.5*(a(1)-a(0));
  /* i == N-1 */        return 0.5*(a(N-1)-a(N-2));
}
//double dh  = edge_width(xc,  jj);     // for the x-loop
//double dhe = edge_width(eps_w, iie);  // for the eps-loop

inline bool too_close(double a, double b, double min_sep_abs){
  const double sep = std::abs(a-b);
  return sep <= min_sep_abs;
}
void enforce_min_spacing(std::map<double,int>& m, double min_sep_abs)
{
  if (m.empty()) return;
  for (auto it = m.begin(); it != m.end(); ) {
    auto nxt = std::next(it);
    if (nxt == m.end()) break;
    if (too_close(it->first, nxt->first, min_sep_abs)) {
      // (policy) keep 'it', erase 'nxt'
      nxt = m.erase(nxt);
    } else {
      it = nxt;
    }
  }
}

/*
vector<int> ImportanceOfPseudo(const Array<double,1>& xn, const Array<complex<double>,2>& Gm, double beta, const Cix& cix, std::array<double,2> cut={0.1,0.5})
{
  bool debug=false;
  Array<double,1> Gwidth(cix.Npseudo_unique);
  Array<double,1> dxh(xn.size()); dxh = GiveDh(xn);
  double x0=0.2;
  Array<double,1> chb(xn.size());
  for (int ix=0; ix<xn.size(); ++ix) chb(ix) = 1.0/cosh(0.5*beta*xn(ix));
  for (int ipp=0; ipp<cix.Npseudo_unique; ++ipp){
    double gi=0;
    for (int ix=0; ix<xn.size(); ++ix){
      double g = fabs(Gm(ipp,ix).imag())*dxh(ix);
      gi += g*chb(ix);
    }
    double sgm=1/gi;
    Gwidth(ipp) = sgm;
    //cout<<"width["<<ipp+1<<"]=" << sgm << endl;
  }
  std::vector<int> idx(cix.Npseudo_unique);
  std::iota(idx.begin(), idx.end(), 0); // 0..N-1
  std::sort(idx.begin(), idx.end(), [&](int a, int b){ return Gwidth(a) < Gwidth(b); });

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
  

double R_OCA(Array<double,2>& Smt, Array<complex<double>,2>& Gc, Array<complex<double>,2>& Fc1, Array<complex<double>,2>& Fc2,
	     bool cmpPhysG, double Q, const Array<double,1>& xn, const Array<complex<double>,2>& Gm, Cix& cix, double T,
	     const Array<double,1>& few, const Array<double,2>& PhiR, const Array<double,1>& om, const Array<double,2>& Deltai, const Array<double,2>& mDelta,
	     const Array<double,1>& feo, json& pars, const vector<int>& G_importance, int OCA_cut_sum, MPIcls& mpi){
  ////////////////////////////
  // input :
  //          Gc(Npseudo,Nw)                                             -- input pseudoparticle propagators
  //          Deltai(Ns,Nom)=Im(Delta(om)), mDelta(Ns,Nom)=Im(Delta(-om))-- input hybridization
  //          xn(Nw)         -- frequency mesh for Gc
  //          om(Nom)        -- frequency for hybridization
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
  //   states[0]=p0 \  s1|    /  states[3]=p3  ibs[3]==ib
  //                  ---<---
  //               Fi[2*is1+1,p3]
  //  states[] for four diagrams, where 0=|0>, 1=|up>, 2=|dn>, 3=|2>
  //  Ds[] 0=up, 1=dn
  bool debug = false;
  const int Nd=4;
  Diag2Matrix diag2matrix(6, cix);
  Array<int,2> diag(2,2);
  diag = 0,2,1,3;
  Array<int,1> dirspin(2);

  //int OCA_cut_sum = pars["OCA_cut_sum"].get<int>();
  /*
  std::array<double,2> OCA_cut({pars["OCA_cut1"].get<double>(),pars["OCA_cut2"].get<double>()});
  vector<int> G_importance = ImportanceOfPseudo(xn, Gm, 1.0/T, cix, OCA_cut);//{0.1,0.5});
  */
  //cout<<"Inside OCA xn.size="<< xn.size() << " " << Gm.extent(1) << endl;
  
  int states_unique_size=0;
  { // This part is carbon-copy of the algorithm below, except that we only calculate the size of the hash-table states_unique.
    // You should not try to understand that, because we comment only the complete algorithm below. This is just stripped/off equivalent,
    // which does not compute anything except total size of the storage.
    unordered_map<OCAKey,int,OCAKeyHash> states_unique;
    states_unique.max_load_factor(0.7f);
    states_unique.reserve(cix.Natst*cix.Natst*cix.Nbath*cix.Nbath);
    int mm_nds_size=0;
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
	    //int njmc1 = cix.Njm_c(ibs[1],iorb1); // index to identify whether njm is proportional to F or not.
	    NjmResult nx1 = cix.Njm_ind(ibs[1],2*iorb1);
	    Array<double,2> njm1;           // matrix when njm is not proportional to F
	    double njmr1 = nx1.njm_r;       // proportionality constant for njm
	    if (nx1.njm_c==2){              // njm is not proportional
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
	    TinyVector<int,Nd> js;
	    for (js[0]=0; js[0]<alw[0].extent(0); ++js[0]){
	      int a0 = cix.inBlock(abs(alw[0](js[0],1))); // alw[0] index is given in terms of atomic state index, but a0 is index within the block. cix.inBlock translates this.
	      int b3 = cix.inBlock(abs(alw[0](js[0],0))); // cix.inBlock translates again to index inside the block ib
	      //double mm0 = sign(alw[0](js[0],1)) * s0_sgn * s1_sgn; // matrix element of operator op(0) <alw[0](j0,1)|op(0)|alw[0](j0,0)> is stored in the sign of op(0) acting on index, hence sign of alw[0](js[0],1) contains op(0) matrix element.
	      for (js[1]=0; js[1]<alw[1].extent(0); ++js[1]){
		int b0 = cix.inBlock(abs(alw[1](js[1],0)));
		int _ii0_ = cix.Sgm_symmetry[ibs[0]](a0,b0); // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
		if (_ii0_==0) continue;
		int ii0 = abs(_ii0_)-1;
		int imp_ii0 = G_importance[ii0];
		int ii_sign0 = sign(_ii0_);                 // Green's function value for off-diagonal propagators could be negative/opposite, hence sign
		int a1 = cix.inBlock(abs(alw[1](js[1],1)));
		//double mm1 = mm0 * sign(alw[1](js[1],1));  // matrix element: <alw[1](j1,1)|op(1)|alw[1](j1,0)>
		for (js[2]=0; js[2]<alw[2].extent(0); ++js[2]){
		  int b1 = cix.inBlock(abs(alw[2](js[2],0)));
		  int _ii1_ = cix.Sgm_symmetry[ibs[1]](a1,b1); // Green's function idex G(a1,b1)==G[|_ii1_|-1]
		  if (_ii1_==0) continue;
		  int ii1 = abs(_ii1_)-1;
		  int imp_ii1 = G_importance[ii1];
		  int ii_sign1 = sign(_ii1_) * ii_sign0;
		  int a2 = cix.inBlock(abs(alw[2](js[2],1)));
		  //double mm2 = mm1 * sign(alw[2](js[2],1));
		  for (js[3]=0; js[3]<alw[3].extent(0); ++js[3]){
		    int b2 = cix.inBlock(abs(alw[3](js[3],0)));
		    int _ii2_ = cix.Sgm_symmetry[ibs[2]](a2,b2); // Green's function index G(a2,b2)==G[|_ii2_|-1]
		    if (_ii2_==0) continue;
		    int ii2 = abs(_ii2_)-1;
		    int imp_ii2 = G_importance[ii2];
		    int ii_sign2 = sign(_ii2_) * ii_sign1;
		    int a3 = cix.inBlock(abs(alw[3](js[3],1)));
		    //double mm3 = mm2 * sign(alw[3](js[3],1));
		    int _ii3_ = cix.Sgm_symmetry[ibs[3]](a3,b3); // Green's function index G(a3,b3)==G[|_ii3_|-1]
		    if (_ii3_==0) continue;
		    int ii3 = abs(_ii3_)-1;
		    int importance = G_importance[ii3] + imp_ii2 + imp_ii1 + imp_ii0;
		    if (importance > OCA_cut_sum ) continue; // this cuts-off PP-states that are less important
		    int ii_sign3 = sign(_ii3_) * ii_sign2;
		    //double mm = mm3*ii_sign3;                    // total sign due to off-diagonal sign (ii_sign3) and matrix elements sign mm3
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
			OCAKey st = pack6(s1, s0, states[3], states[2], states[1], states[0] ); // all atributes of diagram contribution.
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
			OCAKey st = pack6(s1, s0, states[3], states[2], states[1], states[0] ); // all atributes of diagram contribution.
			auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // Did we see that contribution before?
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
  Array<complex<double>,2> _Gc_, _Fc1_, _Fc2_;
  double phi2=0;
  {
    // Data which store the components of the diagram, so that integration looks like a problem of a single orbital.
    //map<TinyVector<int,6>,int,Class1Compare> states_unique; // every unique term is stored in this map with indices {s1,s0,states[3],states[2],states[1], states[0]}
    //// or uncomment below for faster unordered_map, which does not sort the contributions.
    unordered_map<OCAKey,int,OCAKeyHash> states_unique;
    states_unique.reserve(states_unique_size);
    //
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
    vector<vector<TinyVector<int,4>>>  unique_states_inside(4);// unique_states_inside0, unique_states_inside1, unique_states_inside2, unique_states_inside3;
    //vector<vector<int>> indexi(4, vector<int>(states_unique_size, -1));
    //array<vector<int>,4> indexi;
    //for (auto& v : indexi) v.assign(states_unique_size, -1);
    vector<array<int,4>> indexi(states_unique_size,array<int,4>{-1,-1,-1,-1});
    //vector<deque<int>> indexi(4); // index0, index1, index2, index3;
    {
      // index_of_inside can be temporary
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
	      if (nx1.njm_c==2){              // njm is not proportional
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
	      //b2=alw[3](j3,0); njm2 |     v     | a1=alw[1](j1,1)
	      //states[2],ibs[2],ii2   \    |    /  states[1],ibs[1],ii1   
	      //                         ---<---		   
	      //                            njm1           
	      //            a2=alw[2](j2,1)      b1=alw[2](j2,0)             
	      TinyVector<int,Nd> js;
	      for (js[0]=0; js[0]<alw[0].extent(0); ++js[0]){
		int a0 = cix.inBlock(abs(alw[0](js[0],1))); // alw[0] index is given in terms of atomic state index, but a0 is index within the block. cix.inBlock translates this.
		int b3 = cix.inBlock(abs(alw[0](js[0],0))); // cix.inBlock translates again to index inside the block ib
		double mm0 = sign(alw[0](js[0],1)) * s0_sgn * s1_sgn; // matrix element of operator op(0) <alw[0](j0,1)|op(0)|alw[0](j0,0)> is stored in the sign of op(0) acting on index, hence sign of alw[0](js[0],1) contains op(0) matrix element.
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
		      auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // try to insert st-vector as a new key and mm_nds.size() as the value.
		      // If the key is not yet present, inserted=true and "it" points to its current place. If key is present, "it" points to its existing place.
		      int iind = it->second;                    // The key st has value iind, either just inserted, or found from before
		      if (inserted){
			mm_nds[mm_nds_size]=TinyVector<double,5>(mm,0,0,0,0); // if it was just inserted, deque has matrix element for PPself-energy equal mm, and nothing for two particle vertices.
			mm_nds_size+=1;
		      }else
			mm_nds[iind] += TinyVector<double,5>(mm,0,0,0,0); // if we found old entry, which just add to it the current value of mm.
		      // We are done with finding unique contributions. Now we need sub-parts or common parts that appear in the inner-most loop.
		      // 
		      // Below we are finding sub-parts of a diagram that is common to many terms just found above. This is done for efficiency reasons,
		      // so that inside loop integration does not need to repeat the same integration many times.
		      // This subpart contains only states[0](x),states[1](x-e), and hybridization s1(w-x).
		      // The last number 1 means that it is needed for PP-self-energy calculation.
		      // The plot would look like:
		      // 
		      //               op(0)		  		 
		      //              --->---                            |-->---                            		 
		      //states[3]   /w  s1  x \  states[0]	      s1     x \  states[0],ibs[0],ii0    
		      //           |    w-x    | 	             w-x       | b0=alw[1](j1,0)         
		      //           |--e-<--s0--|   op(1)                 |    e <-|   op(1)        	     		 
		      //           |     v     | 		               |        | a1=alw[1](j1,1)	    
		      //states[2]   \w-e | x-e/  states[1]	       v   x-e /  states[1],ibs[1],ii1    
		      //              ---<---		               |--<---			    		 
		      //                                                                              			 
		      //
		      TinyVector<int,4> sti0(states[0],states[1],s1,1);
		      auto [it0, inserted0] = index_of_inside3.try_emplace(sti0, static_cast<int>(unique_states_inside[3].size()));
		      if (inserted0) unique_states_inside[3].push_back(sti0); // first time we've seen this triple
		      else if (unique_states_inside[3][it0->second][3]==0){   // it could happen that it was already found below, but only for two-particle vertex calculation.
			// Now we know that we need it also for PP self-energy calculation, hence we will change the last flag from 0 to 1.
			//cout<<"WARN at sti0="<<unique_states_inside[3][it0->second]<<" was set to 0"<<endl;
			unique_states_inside[3][it0->second][3]=1;
		      }
		      // For second self-energy contribution we need different parts of the diagram in the inside loop. The part is sketched below:
		      //                      ---<--s0----
		      //                      |    w-x    |    op(1)        
		      //b2=alw[3](j3,0)       |x-e      x | a1=alw[1](j1,1)
		      //states[2],ibs[2],ii2   \    v    /  states[1],ibs[1],ii1   
		      //                         ---<---		   
		      TinyVector<int,4> sti1(states[1],states[2],s0,1);
		      auto [it1, inserted1] = index_of_inside0.try_emplace(sti1, static_cast<int>(unique_states_inside[0].size()));
		      if (inserted1) unique_states_inside[0].push_back(sti1); // first time we've seen this triple
		      else if (unique_states_inside[0][it1->second][3]==0){
			//cout<<"WARN at sti1="<<unique_states_inside[0][it1->second]<<" was set to 0"<<endl;
			unique_states_inside[0][it1->second][3]=1;
		      }
		      // For the third self-energy contribution we need yet another part of the diagram in the inside loop. The sketch below:
		      //                         --->---|
		      //states[3],ibs[3],ii3   /x-e    s1
		      //a3=alw[3](j3,1)       |         |
		      //                      |< -e    x-w
		      //b2=alw[3](j3,0); njm2 |         |
		      //states[2],ibs[2],ii2   \ x      v
		      //                         ---<---|
		      TinyVector<int,4> sti2(states[2],states[3],s1,1);
		      auto [it2, inserted2] = index_of_inside1.try_emplace(sti2, static_cast<int>(unique_states_inside[1].size()));
		      if (inserted2) unique_states_inside[1].push_back(sti2); // first time we've seen this triple
		      else if (unique_states_inside[1][it2->second][3]==0){
			//cout<<"WARN at sti2="<<unique_states_inside[1][it2->second]<<" was set to 0"<<endl;
			unique_states_inside[1][it2->second][3]=1;
		      }
		      // For the last forth self-energy we require this part of the diagram:
		      //                            
		      //                         --->---                
		      //states[3],ibs[3],ii3   /    v    \  states[0],ibs[0],ii0
		      //a3=alw[3](j3,1)       |x    e  x-e| b0=alw[1](j1,0)            
		      //                      |           |   op(1)
		      //                      --<-x-w---s0-
		      TinyVector<int,4> sti3(states[3],states[0],s0,1);
		      auto [it3, inserted3] = index_of_inside2.try_emplace(sti3, static_cast<int>(unique_states_inside[2].size()));
		      if (inserted3) unique_states_inside[2].push_back(sti3); // first time we've seen this triple
		      else if (unique_states_inside[2][it3->second][3]==0){
			//cout<<"WARN at sti3="<<unique_states_inside[2][it3->second]<<" was set to 0"<<endl;
			unique_states_inside[2][it3->second][3]=1;
		      }
		      // Here we remember how to access results from unique_states_inside0, and  unique_states_inside1,..,
		      // when we are in the long loop going through all possible contibutions of states_unique
		      if (inserted){
			//indexi[3].push_back(it0->second); // index0
			//indexi[0].push_back(it1->second); // index1
			//indexi[1].push_back(it2->second); // index2
			//indexi[2].push_back(it3->second); // index3
			indexi[iind][3] = it0->second; // index0
			indexi[iind][0] = it1->second; // index1
			indexi[iind][1] = it2->second; // index2
			indexi[iind][2] = it3->second; // index3
		      }else{
			// It could happen that we added st entry first when computing vertex F,
			// and therefore we set index0[] and index1[] to -1.
			// We now found the same st, and instead of inserting it, we just updated it and added mm to it.
			// But therefore we need to check if index0[] or index1[] is set to -1, and correct.
			if (indexi[iind][3]==-1) indexi[iind][3]=it0->second;
			if (indexi[iind][0]==-1) indexi[iind][0]=it1->second;
			if (indexi[iind][1]==-1) indexi[iind][1]=it2->second;
			if (indexi[iind][2]==-1) indexi[iind][2]=it3->second;
		      }
		      // Below we are adding information needed for the two-particle vertex F(om) calculation.
		      // The simplest case of njmc being proportional to psi^+ operator is handled here. This is the case for
		      // one-dimensional block and sometimes other small blocks. But generically the two are not proportional, hence it
		      // must be handled separately with njmc1==2 below.
		      if (nx1.njm_c==1 and cmpPhysG)  // if njm is proportional to F, its simple: Add to the same iind entry, but nds1 part, i.e., second component
			mm_nds[iind] += destroys_Cooper_pair_or_normal_down1 ? TinyVector<double,5>(0,0,0,mm*njmr1,0): TinyVector<double,5>(0,mm*njmr1,0,0,0);
		      if (nx2.njm_c==1 and cmpPhysG)  // if njm is proportional to F its simple: Add to the same iind entry, but last component.
			mm_nds[iind] += destroys_Cooper_pair_or_normal_down2 ? TinyVector<double,5>(0,0,0,0,mm*njmr2): TinyVector<double,5>(0,0,mm*njmr2,0,0);
		    }
		    if (cmpPhysG and nx2.njm_c==2){
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
			  OCAKey st = pack6(s1, s0, states[3], states[2], states[1], states[0] ); // all atributes of diagram contribution.
			  auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // Did we see that contribution before?
			  int _iind_ = it->second; // in mm_nds and indexi, it is stored in _iind_
			  if (inserted){// not seen before, hence insert into mm_nds, index0, index1, index2, index3
			    if (destroys_Cooper_pair_or_normal_down2)
			      mm_nds[mm_nds_size]=TinyVector<double,5>(0,0,0,0,nds2);
			    else
			      mm_nds[mm_nds_size]=TinyVector<double,5>(0,0,nds2,0,0);
			    indexi[mm_nds_size][0]=-1; // should not be needed, hence we are not going to compute contribution of the inner loop.
			    indexi[mm_nds_size][1]=-1; // should not be needed, because inner loop for G&F uses only index2 and index3.
			    indexi[mm_nds_size][2]=-1; // will be added below, but only if its part of PP-self energy contribution
			    indexi[mm_nds_size][3]=-1; // otherwise we choose the configuration such that extra terms vanish
			    mm_nds_size+=1;
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
			  // it implements op(i) operator, giving proper mm matrix element. See comment for mm3 above.
			  double nds1 =  mmp*njm1(_b1_,_a2_); // Now adding njm1.
			  TinyVector<int,4> states(ii0, ii1, ii2, ii3); // atomic states indices between [0,...cix.Npseudo_unique-1] array.
			  OCAKey st = pack6(s1, s0, states[3], states[2], states[1], states[0] ); // all atributes of diagram contribution.
			  auto [it, inserted] = states_unique.try_emplace(st, mm_nds_size); // Did we see that contribution before?
			  int _iind_ = it->second; // in mm_nds and indexi, it is stored in _iind_
			  if (inserted){ // not seen before, hence insert into mm_nds, index0, index1, index2, index3
			    if (destroys_Cooper_pair_or_normal_down1)
			      mm_nds[mm_nds_size]=TinyVector<double,5>(0,0,0,nds1,0);
			    else
			      mm_nds[mm_nds_size]=TinyVector<double,5>(0,nds1,0,0,0);
			    indexi[mm_nds_size][0]=-1; // should not be needed, hence we are not going to compute contribution of the inner loop.
			    indexi[mm_nds_size][1]=-1; // should not be needed, because inner loop for G&F uses only index2 and index3.
			    indexi[mm_nds_size][2]=-1; // will be added below, but only if its part of PP self-energy contribution. We will arrange 
			    indexi[mm_nds_size][3]=-1; // integration such that extra contributions vanish.
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
      (*mpi.out)<<"Final mm_nds_size="<<mm_nds_size<<endl;
      for (const auto& [_st_, iind] : states_unique){
	TinyVector<int,6> st = unpack6(_st_);
	double mm = mm_nds[iind][0];
	double nds1 = abs(mm_nds[iind][1])+abs(mm_nds[iind][3]);
	double nds0 = abs(mm_nds[iind][2])+abs(mm_nds[iind][4]);
	int s1 = st[0];
	int s0 = st[1];
	TinyVector<int,Nd> states(st[5],st[4],st[3],st[2]);
	TinyVector<int,4> id(indexi[iind][0],indexi[iind][1],indexi[iind][2],indexi[iind][3]);
	//
	// unique_states_inside[y][indexi[iind][y]] connected to st.
	//
	if (mm!=0){
	  TinyVector<int,4> sti0 = unique_states_inside[0][id[0]];
	  if (sti0[0]!=states[1] or sti0[1]!=states[2] or sti0[2]!=s0){
	    cout<<"WARNING0 sti0="<<sti0<<" while states[1],states[2],s0="<<states[1]<<","<<states[2]<<","<<s0<<" and st="<<st<<endl;
	  }
	  TinyVector<int,4> sti1 = unique_states_inside[1][id[1]];
	  if (sti1[0]!=states[2] or sti1[1]!=states[3] or sti1[2]!=s1){
	    cout<<"WARNING1 sti1="<<sti1<<" while states[2],states[3],s1="<<states[2]<<","<<states[3]<<","<<s1<<" and st="<<st<<endl;
	  }
	  TinyVector<int,4> sti2 = unique_states_inside[2][id[2]];
	  if (sti2[0]!=states[3] or sti2[1]!=states[0] or sti2[2]!=s0){
	    cout<<"WARNING2 sti2="<<sti2<<" while states[3],states[0],s0="<<states[3]<<","<<states[0]<<","<<s0<<" and st="<<st<<endl;
	  }
	  TinyVector<int,4> sti3 = unique_states_inside[3][id[3]];
	  if (sti3[0]!=states[0] or sti3[1]!=states[1] or sti3[2]!=s1){
	    cout<<"WARNING3 sti3="<<sti3<<" while states[0],states[1],s1="<<states[0]<<","<<states[1]<<","<<s1<<" and st="<<st<<endl;
	  }
	}
	if (cmpPhysG){
	  // This is the crucial part for computing F(omega). In this case we have many additional contributions to states_unique,
	  // for which mm=0 but nds1 or nsd0 is finite. For these cases, we need to find inside integral contributions, such that
	  // we can evaluate inside contributions first. It turns out that we can arrange integrations such that the number of
	  // inside contributions does not increase. The outside contributions are much more numerous, but this is not computationally
	  // very expensive.
	  // In this block below, we just find indexi[2] or indexi[3] (for inside integral) for each term in states_unique, in particular for those
	  // that don't contribute to PP-self-energy.
	  if (nds1!=0){
	    if (id[2]<0){
	      auto it2 = index_of_inside2.find(TinyVector<int,4>(states[3],states[0],s0,1));
	      if (it2!=index_of_inside2.end()){
		indexi[iind][2] = it2->second;
	      }else{
		auto it2 = index_of_inside2.find(TinyVector<int,4>(states[3],states[0],s0,0));
		if (it2!=index_of_inside2.end()){
		  indexi[iind][2] = it2->second;
		  //cout<<" iind="<<iind<<" also gets corrected with indexi["<<iind<<"][2]="<<indexi[iind][2]<<" sti2="<<it2->first<<endl;
		}else{
		  // Added code to correct for the case diagrams are prunned early and therefore this inside contribution newer appears due to prunning.
		  TinyVector<int,4> sti2(states[3],states[0],s0,0); // Here last 0 means that it is needed only for CmpPhysG
		  int isize = static_cast<int>(unique_states_inside[2].size());
		  //cout<<"Correcting iind="<<iind<<" index_of_inside2["<<sti2<<"]="<<isize<<" indexi["<<iind<<"][2]="<<isize<<endl;
		  index_of_inside2[sti2] = isize;
		  unique_states_inside[2].push_back(sti2);
		  indexi[iind][2] = isize;
		}
	      }
	    }
	  }
	  if (nds0!=0){
	    if (id[3]<0){
	      auto it3 = index_of_inside3.find(TinyVector<int,4>(states[0],states[1],s1,1));
	      if (it3!=index_of_inside3.end()){
		indexi[iind][3] = it3->second;
	      }else{
		auto it3 = index_of_inside3.find(TinyVector<int,4>(states[0],states[1],s1,0));
		if (it3!=index_of_inside3.end()){
		  indexi[iind][3] = it3->second;
		  //cout<<" iind="<<iind<<" also gets corrected with indexi["<<iind<<"][3]="<<indexi[iind][3]<<" sti3="<<it3->first<<endl;
		}else{
		  // Added code to correct for the case diagrams are prunned early and therefore this inside contribution newer appears due to prunning.
		  TinyVector<int,4> sti3(states[0],states[1],s1,0); // Here last 0 means that it is needed only for CmpPhysG
		  int isize = static_cast<int>(unique_states_inside[3].size());
		  //cout<<"Correcting iind="<<iind<<" index_of_inside3["<<sti3<<"]="<<isize<<" indexi["<<iind<<"][3]="<<isize<<endl;
		  index_of_inside3[sti3] = isize;
		  unique_states_inside[3].push_back(sti3);
		  indexi[iind][3] = isize;
		}
	      }
	    }
	  }
	}
      }
    
      {// Checking for possible ERROR in assigning indexi[iind][?]
	int ijk=0;
	for (const auto& [_st_, iind] : states_unique){
	  TinyVector<int,6> st = unpack6(_st_);
	  double mm = mm_nds[iind][0];
	  double nds1 = abs(mm_nds[iind][1])+abs(mm_nds[iind][3]);
	  double nds0 = abs(mm_nds[iind][2])+abs(mm_nds[iind][4]);
	  int s1 = st[0];
	  int s0 = st[1];
	  TinyVector<int,Nd> states(st[5],st[4],st[3],st[2]);
	  TinyVector<int,4> id(indexi[iind][0],indexi[iind][1],indexi[iind][2],indexi[iind][3]);
	  if (mm!=0){
	    TinyVector<int,4> sti0 = unique_states_inside[0][id[0]];
	    if (sti0[0]!=states[1] or sti0[1]!=states[2] or sti0[2]!=s0){
	      cout<<"WARNING0 sti0="<<sti0<<" while states[1],states[2],s0="<<states[1]<<","<<states[2]<<","<<s0<<" and st="<<st<<endl;
	    }
	    TinyVector<int,4> sti1 = unique_states_inside[1][id[1]];
	    if (sti1[0]!=states[2] or sti1[1]!=states[3] or sti1[2]!=s1){
	      cout<<"WARNING1 sti1="<<sti1<<" while states[2],states[3],s1="<<states[2]<<","<<states[3]<<","<<s1<<" and st="<<st<<endl;
	    }
	    TinyVector<int,4> sti2 = unique_states_inside[2][id[2]];
	    if (sti2[0]!=states[3] or sti2[1]!=states[0] or sti2[2]!=s0){
	      cout<<"WARNING2 sti2="<<sti2<<" while states[3],states[0],s0="<<states[3]<<","<<states[0]<<","<<s0<<" and st="<<st<<endl;
	    }
	    TinyVector<int,4> sti3 = unique_states_inside[3][id[3]];
	    if (sti3[0]!=states[0] or sti3[1]!=states[1] or sti3[2]!=s1){
	      cout<<"WARNING3 sti3="<<sti3<<" while states[0],states[1],s1="<<states[0]<<","<<states[1]<<","<<s1<<" and st="<<st<<endl;
	    }
	  }
	  if (nds1!=0){
	    if (id[2]<0){
	      cout<<"ERROR id[2]="<<id[2]<<" at st="<<st<<" iind="<<iind<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
	      TinyVector<int,4> sti2(states[3],states[0],s0,1); // this inner loop is needed for G&F calculation
	      auto it2 = index_of_inside2.find(sti2);
	      if (it2!=index_of_inside2.end()){
		indexi[iind][2] = it2->second;
		cout<<"  Found ind="<<it2->second<<" where st="<<unique_states_inside[2][it2->second]<<endl;
	      }
	    }
	  }
	  if (nds0!=0){
	    if (id[3]<0){
	      cout<<"ERROR id[3]="<<id[3]<<" at st="<<st<<" iind="<<iind<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
	      TinyVector<int,4> sti3(states[0],states[1],s1,1); // this inner loop is needed for G&F calculation
	      auto it3 = index_of_inside3.find(sti3);
	      if (it3!=index_of_inside3.end()){
		indexi[iind][3] = it3->second;
		cout<<"  Found ind="<<it3->second<<" where st="<<unique_states_inside[3][it3->second]<<endl;
	      }
	    }
	  }
	  ijk+=1;
	}
      }
      if (debug){
	int first_bath_anomalous=cix.Nbath;
	for (int ib=0; ib<cix.Nbath; ++ib)
	  if (cix.bath_SC(ib) and abs(cix.bath_sigind0(ib))<first_bath_anomalous)
	    first_bath_anomalous = abs(cix.bath_sigind0(ib));
    
	// Below we print out all contributions for OCA, i.e., all possible triple integrations we need to do.
	for (int i=0; i<unique_states_inside[0].size(); ++i){
	  TinyVector<int,4> sti0 = unique_states_inside[0][i];
	  if ((sti0[0]>=cix.first_anomalous-1 or sti0[1]>=cix.first_anomalous-1) and sti0[2]>=first_bath_anomalous)
	    fmt::print(cout,"#R0: {} ii1={} ii2={} s0={} pp={}\n", i, sti0[0]+1, sti0[1]+1, sti0[2], sti0[3]);
	}
	for (int i=0; i<unique_states_inside[1].size(); ++i){
	  TinyVector<int,4> sti1 = unique_states_inside[1][i];
	  if ((sti1[0]>=cix.first_anomalous-1 or sti1[1]>=cix.first_anomalous-1) and sti1[2]>=first_bath_anomalous)
	    fmt::print(cout,"#R1: {} ii2={} ii3={} s1={} pp={}\n", i, sti1[0]+1, sti1[1]+1, sti1[2], sti1[3]);
	}
	for (int i=0; i<unique_states_inside[2].size(); ++i){
	  TinyVector<int,4> sti2 = unique_states_inside[2][i];
	  if ((sti2[0]>=cix.first_anomalous-1 or sti2[1]>=cix.first_anomalous-1) and sti2[2]>=first_bath_anomalous)
	    fmt::print(cout,"#R2: {} ii3={} ii0={} s0={} pp={}\n", i, sti2[0]+1, sti2[1]+1, sti2[2], sti2[3]);
	}
	for (int i=0; i<unique_states_inside[3].size(); ++i){
	  TinyVector<int,4> sti3 = unique_states_inside[3][i];
	  if ((sti3[0]>=cix.first_anomalous-1 or sti3[1]>=cix.first_anomalous-1) and sti3[2]>=first_bath_anomalous)
	    fmt::print(cout,"#R3: {} ii0={} ii1={} s1={} pp={}\n", i, sti3[0]+1, sti3[1]+1, sti3[2], sti3[3]);
	}
	
	int ijk=0;
	for (const auto& [_st_, iind] : states_unique){
	  TinyVector<int,6> st = unpack6(_st_);
	  double mm = mm_nds[iind][0];
	  double nds1 = abs(mm_nds[iind][1])+abs(mm_nds[iind][3]);
	  double nds0 = abs(mm_nds[iind][2])+abs(mm_nds[iind][4]);
	  int s1 = st[0];
	  int s0 = st[1];
	  TinyVector<int,Nd> states(st[5],st[4],st[3],st[2]);
	  TinyVector<int,4> id(indexi[iind][0],indexi[iind][1],indexi[iind][2],indexi[iind][3]);
	  if (mm!=0){
	    TinyVector<int,4> sti0 = unique_states_inside[0][id[0]];
	    if (sti0[0]!=states[1] or sti0[1]!=states[2] or sti0[2]!=s0){
	      cout<<"WARNING0 sti0="<<sti0<<" while states[1],states[2],s0="<<states[1]<<","<<states[2]<<","<<s0<<" and st="<<st<<endl;
	    }
	    TinyVector<int,4> sti1 = unique_states_inside[1][id[1]];
	    if (sti1[0]!=states[2] or sti1[1]!=states[3] or sti1[2]!=s1){
	      cout<<"WARNING1 sti1="<<sti1<<" while states[2],states[3],s1="<<states[2]<<","<<states[3]<<","<<s1<<" and st="<<st<<endl;
	    }
	    TinyVector<int,4> sti2 = unique_states_inside[2][id[2]];
	    if (sti2[0]!=states[3] or sti2[1]!=states[0] or sti2[2]!=s0){
	      cout<<"WARNING2 sti2="<<sti2<<" while states[3],states[0],s0="<<states[3]<<","<<states[0]<<","<<s0<<" and st="<<st<<endl;
	    }
	    TinyVector<int,4> sti3 = unique_states_inside[3][id[3]];
	    if (sti3[0]!=states[0] or sti3[1]!=states[1] or sti3[2]!=s1){
	      cout<<"WARNING3 sti3="<<sti3<<" while states[0],states[1],s1="<<states[0]<<","<<states[1]<<","<<s1<<" and st="<<st<<endl;
	    }
	  }
	  fmt::print(cout," {:6d} ii0={:<3d},ii1={:<3d},ii2={:<3d},ii3={:<3d},s0={:d},s1={:d}: mm={:2g} nds1={:2g},{:2g} nds0={:2g},{:2g} index0={:4d} index1={:4d} index2={:4d} index3={:4d} ",
		     ijk,st[5]+1,st[4]+1,st[3]+1,st[2]+1,st[1],st[0],
		     mm,mm_nds[iind][1],mm_nds[iind][3],mm_nds[iind][2],mm_nds[iind][4],id[0],id[1],id[2],id[3]);
	  fmt::print(cout,"\n");

	  ijk+=1;
	}
      }
    }
    
    Spline1D<double> fex(xn.size()), fep(om.size());
    for (int i=0; i<xn.size(); ++i) fex[i]=few(i);
    for (int i=0; i<om.size(); ++i) fep[i]=feo(i);
    fex.splineIt(xn);
    fep.splineIt(om);
  
    Array<complex<double>,1> Gx(cix.Npseudo_unique), Gw(cix.Npseudo_unique);
    Array<double,1> Gxt(cix.Npseudo_unique), Gwt(cix.Npseudo_unique);;
    Array<double,1> Dip(cix.Nbath_unique), Dim(cix.Nbath_unique);
  
    phi2=0;
    Array<double,1> xc(2*xn.size()+om.size());
    Array<int,1> xi;
    {
      Array<double,2> TSmt(xn.size(),cix.Npseudo_unique); TSmt=0;
      vector<Array<double,2>> Ri(4);
      //vector<Array<double,3>> Ris(4);
      for (int ir=0; ir<4; ++ir) Ri[ir].resize(unique_states_inside[ir].size(),3);

      // Below we decide how to redistribute jobs across the mpi run.
    
      std::vector<BlockJob> jobs;
      jobs.reserve(xn.size());
      for (int iw=0; iw<xn.size(); ++iw){ // this is for computing PPSigma
	double x = iw/(xn.size()+0.0); // we can estimate complexity of integrations
	jobs.push_back({iw, 0, 0.6+tanh(15*x*(1.0-x)) }); // the formula comes from the fact that the frequencies at the edge have smaller integration range.
      }
      std::vector<std::vector<BlockJob>> assign(mpi.size);
      if (mpi.size==1)
	assign[0]=jobs;
      else
	mpi.Assign_jobs(assign, jobs); // Here we actually assign jobs.

      for (int ij=0; ij<assign[mpi.rank].size(); ++ij){
	int iom = assign[mpi.rank][ij].id;
	//cout<<"starting iom="<<iom<<endl;
	auto t0 = high_resolution_clock::now();
      
	double w=xn(iom);
	double fe_w = few(iom);
      
	// requirements for m_eps mesh : eps>=om(0)+w and eps<=om(-1)+w and eps>=xn(0) and eps<=xn(-1)
	double om_m1 = om(om.size()-1);
	double xn_m1 = xn(xn.size()-1);
      
	map<double,int> m_x0;  // m_x0 is dense around 0 as pp mesh, and has also mesh around om as om(i) mesh describing Delta.
	for (int i=0; i<xn.size(); ++i){
	  double eps = xn(i); 
	  if (eps-w>=om(0) and eps-w<=om_m1) m_x0.insert({eps,i});
	}
	for (int i=0; i<om.size(); ++i){
	  double eps = om(i)+w;
	  if(eps>=xn(0) and eps<=xn_m1) m_x0.insert({eps,-i-1}); // adding to m_eps -> om(i)+w points
	}
	enforce_min_spacing(m_x0, 1e-10);
	map<double,int> m_eps;  // m_eps is dense around w-eps as pp mesh, and has also mesh around 0 as om(i) mesh describing Delta.
	for (int i=0; i<om.size(); ++i){
	  double eps = om(i); 
	  if (w-eps>=xn(0) and w-eps<=xn_m1) m_eps.insert({eps,-i-1});
	}
	for (int i=0; i<xn.size(); ++i){
	  double eps = w-xn(i);
	  if(eps>=om(0) and eps<=om_m1) m_eps.insert({eps,i}); 
	}
	enforce_min_spacing(m_eps, 1e-10);
	Array<double,1> eps_w(m_eps.size()); // mesh for eps-om, where om is mesh for physical spectra and eps mesh for pseudo spectra,
	std::size_t iie=0;
	for (const auto &[eps,iep] : m_eps){ // going over all combined points for eps
	  eps_w(iie++) = eps;
	}
	//Array<double,1> ff_w (m_eps.size()); // contains f(eps-w)*f(-eps)/f(-w) with eps dense around zero and at om
	int posb=xn.size()-1, posc=0, mposc=om.size()-1;
	iie=-1;
	for (const auto &[eps,iep] : m_eps){ // going over all combined points for eps
	  iie+=1;
	  if (eps-w<om(0) or eps-w>om(om.size()-1)) continue;//this should not happen anymore
	  double f_eps_w, f_eps; // f(eps-w), f(eps)
	
	  double ffx=0;
	  if (iep>=0){  // point is from mesh xn but as xn=w-eps
	    // I think f_eps_w should be 1-few(iep). Check below. correct.
	    f_eps_w=fermi((eps-w)/T); // 1-few(iep);
	    if (abs(1-few(iep)-f_eps_w)>1e-10){
	      cout<<"ERROR 1-few(iep)=1-few(iep) while f_eps_w="<<f_eps_w<<endl;
	    }
	    f_eps=fermi(eps/T);
	    if (w>0) // for positive w, we compute f(eps-w)*f(-eps)/f(-w)
	      ffx = f_eps_w*(1-f_eps)/(1-fe_w);
	    else     // for negative w, f(w-eps)*f(eps)/f(w)
	      ffx = (1-f_eps_w)*f_eps/fe_w;
	  }else{ // point is from mesh om[:]
	    f_eps_w = fermi((eps-w)/T);
	    f_eps = fermi(eps/T); // should be feo(-iep-1) //fermi(eps/T);
	    if (fabs(f_eps-feo(-iep-1))>1e-6){
	      cout<<"ERROR f_eps="<<f_eps<<" while feo(-iep-1)="<<feo(-iep-1)<<endl;
	    }
	    if (w>0){
	      ffx = f_eps_w*(1-f_eps)/(1-fe_w);
	    }else{
	      ffx = (1-f_eps_w)*f_eps/fe_w;
	    }
	  }
	  //ff_w (iie) = ffx; // f(eps-w)*f(-eps)/f(-w)
	  int Nx=xn.size();
	  int No=om.size();
	  map<double,int> m_x(m_x0); // copy m_x0 mesh into m_x, hence it is dense around 0 (as pp), and around w (as Delta). We will next add mesh dens around eps (as pp).
	  deque<double> keys_to_erase;
	  for (const auto &[x,ix] : m_x){
	    // requirements for m_eps mesh which is also obeyed by m_x now: eps>=om(0)+w and eps<=om(-1)+w and eps>=xn(0) and eps<=xn(-1)
	    if(x-eps<xn(0) or x-eps>xn(Nx-1)) keys_to_erase.push_back(x);
	  }
	  for (int l=0; l<keys_to_erase.size(); ++l) m_x.erase(keys_to_erase[l]);
	  for (int i=0; i<xn.size(); ++i){
	    double x = xn(i)+eps;
	    if(x>=xn(0) and x<=xn(Nx-1) and x-w>=om(0) and x-w<=om(No-1)) m_x.insert({x,10000 + i});
	  }
	  enforce_min_spacing(m_x, 1e-10);
	  xc.resize(m_x.size());
	  int jj=0;
	  for (const auto &[x,ix] : m_x){
	    xc(jj)=x;
	    jj+=1;
	  }
	  for (int ir=0; ir<4; ++ir) Ri[ir]=0;
	  int pos_x=0, pos_o=0, pos_y=0, pos_mo=om.size()-1, pos_z=0;
	  jj=-1;
	  for (const auto &[x,ix] : m_x){ //x is dense around 0 (as pp), and around w (as Delta) and dense around eps (as pp)
	    jj+=1;
	    if(x<xn(0) or x>xn(Nx-1) or x-eps<xn(0) or x-eps>xn(Nx-1) or x-w<om(0) or x-w>om(om.size()-1)) continue;
	    double dh=0;
	    if (jj>0 and jj<xc.size()-1){
	      dh = 0.5*(xc(jj+1)-xc(jj-1));
	    }else{
	      int N=xc.size();
	      dh = (jj==0) ? 0.5*(xc(1)-xc(0)) : 0.5*(xc(N-1)-xc(N-2));
	    }
	    intpar px=InterpLeft(x,pos_x,xn);
	    intpar py=InterpLeft(x-eps,pos_y,xn);
	    intpar po=InterpLeft(x-w,pos_o,om);
	    intpar mpo=InterpRight(w-x,pos_mo,om);
	
	    double f_x     = fex(px);  // f(x)
	    // Note that Gx[0] and Gxt[0] is not computed and should not be used.
	    for (int ip=0; ip<cix.Npseudo_unique; ++ip){
	      Gx(ip) = linInterp(Gm,ip,px);  // Gx[ip]  = Gm[ip](x)
	      // We will not compute the complex green's function (Re(G),Im(G)), but we will keep our previous choice
	      // that contains (Re(G), -Im(G)/f(-w) ). From that we will below compute Im(G), only where needed.
	      //Gxt(ip) = Gx(ip).imag();       // Gxt[ip] = Gt[ip](x)
	      //Gx(ip) = complex<double>(Gx(ip).real(), -Gx(ip).imag()*(1-f_x)); // Gx[ip] = Gc[states[ip]](x)
	    }
	    double f_x_eps = fex(py);  // f(x-eps)
	    for (int ip=0; ip<cix.Npseudo_unique; ++ip){
	      Gw(ip) = linInterp(Gm,ip,py);  // Gw[ip]  = Gm[ip](x-eps)
	      //Gwt(ip) = Gw(ip).imag();       // Gwt[ip] = Gt[ip](x-eps)
	      //Gw(ip) = complex<double>(Gw(ip).real(), -Gw(ip).imag()*(1-f_x_eps)); // Gw[ip] = Gc[states[ip]](x-eps)
	    }
	    intpar pz=InterpLeft(x-w,pos_z,xn);
	    double f_x_w   = fex(pz);  // f(x-w)
	    double ffw = (w>0) ? f_x_w*(1-f_x)/(1-fe_w) : (1-f_x_w)*f_x/fe_w; // f(x-om)*f(-x)/f(-om)
	    double ffy = (eps>w) ? (1-f_x_w)*f_x_eps/(1-f_eps_w) : f_x_w*(1-f_x_eps)/f_eps_w; // f(x-w)*f(eps-x)/f(eps-w)
	    for (int is=0; is<cix.Nbath_unique; ++is){
	      Dip(is) = linInterp(Deltai,is,po)*dh;   // D(x-w)*dh
	      Dim(is) = linInterp(Deltai,is,mpo)*dh;  // D(w-x)*dh
	    }
	    for (int ir=0; ir<4; ir++){
	      //  R[0](:,:)=
	      //           ---<--s0----
	      //           |    w-x    |
	      //           |x-e      x | 
	      //states[2]   \    v    /  states[1]
	      //              ---<---
	      //
	      //  R[1](:,:)=
	      //            --->---|
	      //states[3]  /x-e    s1
	      //          |        |
	      //          |< -e    x-w
	      //          |        |
	      //states[2]  \ x     v
	      //            ---<---|
	      //
	      //  R[2](:,:)=                          
	      //            --->---                
	      //states[3] /    v    \  states[0]
	      //         |x    e  x-e| 
	      //         |           |
	      //         --<-x-w---s0-
	      //
	      // R[3](:,:)= 		 
	      //       |-->---                            		 
	      //      s1     x \  states[0]
	      //      w-x       |
	      //       |    e <-|
	      //       |        |
	      //       v   x-e /  states[1]
	      //       |--<---			    		 
	      for (int i=0; i<unique_states_inside[ir].size(); ++i){
		TinyVector<int,4> st = unique_states_inside[ir][i];
		if (st[3]==0) continue; // now it can happen due to prunning. These contibutions are needed only in CmpPhysG.
		double Ddh = (ir==1 or ir==2) ? Dip(st[2]) : Dim(st[2]); // either D(x-w) or D(w-x) depending on direction of hybridization. Needs to look at the diagram.
		complex<double> G0 = Gx(st[0]), G1 = Gw(st[1]); // this is (G.real, Gt) not (G.real,G.imag)
		// This should be Ddh*real(G0*G1), but because G0.imag==-G0t*f(-x) and G1.imag==-G1t*f(-x+eps), we have to add these fermi functions.
		Ri[ir](i,0) += Ddh*f_x_w * (G0.real()*G1.real()-G0.imag()*(1-f_x)*G1.imag()*(1-f_x_eps)); 
		Ri[ir](i,1) += Ddh*ffw   * G0.imag()*G1.real();
		Ri[ir](i,2) += Ddh*ffy   * G0.real()*G1.imag();
	      }
	    }
	  } // end ix
	  if (fabs(eps-eps_w(iie))>1e-10){ cout<<"ERROR eps="<<eps<<" eps_w(iie)="<<eps_w(iie)<<endl;}
	  int Ne=eps_w.extent(0);
	  double dh=0;
	  if (iie>0 and iie<Ne-1)
	    dh = 0.5*(eps_w(iie+1)-eps_w(iie-1));
	  else
	    dh = (iie==0) ? 0.5*(eps_w(1)-eps_w(0)) : 0.5*(eps_w(Ne-1)-eps_w(Ne-2));
	  //
	  intpar pb=InterpRight(w-eps,posb,xn); // could be iep
	  intpar pc=InterpLeft(eps,posc,om); // could be iep
	  intpar mpc=InterpRight(-eps,mposc,om);
	  for (int ip=0; ip<cix.Npseudo_unique; ++ip){
	    Gx(ip) = (iep>=0) ? Gm(ip,iep) : linInterp(Gm,ip,pb);  // Gx[ip]  = Gm[ip](w-eps); mesh is such that iep>0, than eps=w-xn(iep)
	    Gxt(ip) = Gx(ip).imag();                               // Gxt[ip] = Gt[ip](w-eps)
	    Gx(ip) = complex<double>(Gx(ip).real(), -Gx(ip).imag()*(1-f_eps)); // Gx[ip] = Gc[ip](w-eps)
	  }
	  for (int is=0; is<cix.Nbath_unique; ++is){
	    Dip(is) = (iep<0) ? Deltai(is,-iep-1)*dh : linInterp(Deltai,is,pc)*dh;// D(eps)*dh
	    Dim(is) = linInterp(Deltai,is,mpc)*dh;                                // D(-eps)*dh
	  }

	  //if (mpi.rank==0) cout<<"eps="<<eps<<" done"<<endl;
	
	  double f_meps = 1-f_eps; // f(-eps)
	  Array<double,2> Rix(4,3);
	  for (auto const& [_st_, iind] : states_unique){
	    TinyVector<int,6> st = unpack6(_st_);
	    // Smt[states[0]]=
	    //               --->---
	    //states[3] w-e/ -e |   w\  states[0]
	    //            |     |     |
	    //            --<-s0--w-x--
	    //  R[0](e,w) |x-e  |   x | 
	    // states[2]   \    v    /  states[1]
	    //               ---<---
	    //
	    // Smt[states[1]]=
	    //            --->--|---
	    //states[3]  /x-e   s1   \ w-e states[0]
	    //          |       |     |
	    //R[1](e,w) |< -e --s0----|
	    //          |      x-w    |
	    //states[2]  \ x    v    / w  states[1]
	    //             --<--|---
	    //
	    // Smt[states[2]]=
	    //           R[2](e,w)                          
	    //            --->----                
	    //states[3] /    v     \  states[0]
	    //         |x    |  x-e | 
	    //         --<-x-w---s0-
	    //         |     e      |
	    //states[2] \w   v w-e /  states[1]
	    //            --------
	    //
	    // Smt[states[3]]=
	    //             ----|-->---                            		 
	    //states[3]  /w   s1    x \  states[0]
	    //          |     w-x      |
	    //          |--e-- |--<----|  R[3](e,w)
	    //          |      |       |
	    //states[2]  \w-e  v  x-e /  states[1]
	    //             ----|--<--			    		 
	    double mm = -mm_nds[iind][0]; // minus comes from the crossing of hybridization lines
	    if (mm==0) continue;
	    int s0=st[1], s1=st[0];
	    TinyVector<int,4> states(st[5],st[4],st[3],st[2]);
	    TinyVector<int,4> id(indexi[iind][0],indexi[iind][1],indexi[iind][2],indexi[iind][3]);
	    TinyVector<double,4> Gt( Gxt(states[0]), Gxt(states[1]), Gxt(states[2]), Gxt(states[3]) );
	    TinyVector<double,4> Gr( Gx(states[0]).real(), Gx(states[1]).real(), Gx(states[2]).real(), Gx(states[3]).real() );
	    for (int ir=0; ir<4; ++ir)
	      for (int ip=0; ip<=2; ++ip)
		Rix(ir,ip) =  Ri[ir](id[ir],ip);                                                                // ffx = f(-eps)*f(eps-w)/f(-w)
	    TSmt(iom,states[0]) += Dim(s1)*mm*(ffx*(Gt[3]*Rix(0,0) + Gr[3]*Rix(0,2)) + f_meps*Gr[3]*Rix(0,1));  // Sig0 += D1(-eps)*(ffx*(Gt3(w-eps)*R0_0(eps,w)+Re(G3(w-eps))*R0_2(eps,w))+f(-eps)*Re(G3(w-eps))*R0_1(eps,w))
	    TSmt(iom,states[1]) += Dim(s0)*mm*(ffx*(Gt[0]*Rix(1,0) + Gr[0]*Rix(1,2)) + f_meps*Gr[0]*Rix(1,1));  // Sig1 += D0(-eps)*(ffx*(Gt0(w-eps)*R1_0(eps,w)+Re(G0(w-eps))*R1_2(eps,w))+f(-eps)*Re(G0(w-eps))*R1_1(eps,w))
	    TSmt(iom,states[2]) += Dip(s1)*mm*(ffx*(Gt[1]*Rix(2,0) + Gr[1]*Rix(2,2)) + f_meps*Gr[1]*Rix(2,1));  // Sig2 += D1( eps)*(ffx*(Gt1(w-eps)*R2_0(eps,w)+Re(G1(w-eps))*R2_2(eps,w))+f(-eps)*Re(G1(w-eps))*R2_1(eps,w))
	    TSmt(iom,states[3]) += Dip(s0)*mm*(ffx*(Gt[2]*Rix(3,0) + Gr[2]*Rix(3,2)) + f_meps*Gr[2]*Rix(3,1));  // Sig3 += D0( eps)*(ffx*(Gt2(w-eps)*R3_0(eps,w)+Re(G2(w-eps))*R3_2(eps,w))+f(-eps)*Re(G2(w-eps))*R3_1(eps,w))
	  }
	}
	if (mpi.rank==mpi.master){
	  auto t1 = high_resolution_clock::now();
	  fsec dt = t1-t0;
	  fmt::print(*mpi.out,"{}({:.4f}) ", iom, dt.count());
	  mpi.out->flush();
	}
      }
      if (mpi.rank==mpi.master) *mpi.out<<endl;
    
      if (mpi.size>1){
	mpi.AllReduceSum_simple_array(TSmt);
      }
      Array<double,1> dxh(xn.size()); dxh = GiveDh(xn);
      Smt.resize(cix.Npseudo_unique,xn.size());
      for (int ip=0; ip<cix.Npseudo_unique; ++ip){
	Smt(ip,Range::all()) = TSmt(Range::all(),ip)/cix.deg_pp(ip);
	phi2 += sum(Smt(ip,Range::all())*PhiR(ip,Range::all())*dxh(Range::all()))*cix.deg_pp(ip);
      }
    }
    phi2 /= (4*Q); // because there are 2-pp propagators, hence Phi=Tr(Sigma*G)/(2*n)

    if (cmpPhysG){
      std::vector<BlockJob> jobs;
      jobs.reserve(om.size());
      for (int iw=0; iw<om.size(); ++iw){ // this is for computing PPSigma
	double x = iw/(om.size()+0.0); // we can estimate complexity of integrations
	jobs.push_back({iw, 0, 0.6+tanh(15*x*(1.0-x)) }); // the formula comes from the fact that the frequencies at the edge have smaller integration range.
      }
      std::vector<std::vector<BlockJob>> assign(mpi.size);
      if (mpi.size==1)
	assign[0]=jobs;
      else
	mpi.Assign_jobs(assign, jobs); // Here we actually assign jobs.
    
      _Gc_.resize(cix.Nbath_unique_unique,om.size()); _Gc_=0;
      _Fc1_.resize(cix.Nbath_unique_unique,om.size()); _Fc1_=0;
      _Fc2_.resize(cix.Nbath_unique_unique,om.size()); _Fc2_=0;
    
      vector<Array<complex<double>,2>> Ri(2);
      Ri[0].resize(unique_states_inside[2].size(),4);
      Ri[1].resize(unique_states_inside[3].size(),4);
      Array<complex<double>,2> Rix(2,4);

      //auto [iws,iwe] = mpi.mpiSplitArray(om.size());
      //int iws=0, iwe=om.size();
      //for (int iom=iws; iom<iwe; ++iom){
      for (int ij=0; ij<assign[mpi.rank].size(); ++ij){
	int iom = assign[mpi.rank][ij].id;
	//cout<<"cmpPhysG iom="<<iom<<endl;
	auto t0 = high_resolution_clock::now();
	double w=om(iom);
	double fe_w = feo(iom);
	map<double,int> m_eps;  // m_eps is dense around 0 as pp mesh is, and has also mesh around w with pp-mesh density
	// requirements for m_eps mesh : eps>=om(0)+w and eps<=om(-1)+w and eps>=xn(0) and eps<=xn(-1)
	int Nx = xn.size();
	//requirement for m_eps : xn(0)<=eps+w<=xn(-1)  and xn(0)<=eps<=xn(-1)
	for (int iep=0; iep<xn.size(); ++iep){
	  double eps = xn(iep); 
	  if (eps+w>=xn(0) and eps+w<=xn(Nx-1)) m_eps.insert({eps,iep});
	}
	for (int i=0; i<xn.size(); ++i){
	  double eps = xn(i)+w;
	  if(eps>=xn(0) and eps<=xn(Nx-1))
	    m_eps.insert({eps,-i-1}); // adding to m_eps -> xn(i)+w points
	}
	enforce_min_spacing(m_eps, 1e-10);
	Array<double,1> eps_w(m_eps.size()); // mesh for eps-om, where om is mesh for physical spectra and eps mesh for pseudo spectra,
	int iie=0;
	for (const auto &[eps,iep] : m_eps){
	  eps_w(iie) = eps;
	  iie+=1;
	}
	int pos_w=0, pos_ex=0;
	for (int iie=0; iie<eps_w.size(); ++iie){
	  double eps = eps_w(iie);
	  int Nx=xn.size();
	  int No=om.size();
	  map<double,int> m_x; // m_x requires : xn(0)<=x-w<=xn(-1) and xn(0)<=x<=xn(-1) and om(0)<=x-eps<=om(-1)
	  for (int i=0; i<xn.size(); ++i){
	    double x = xn(i); 
	    if (x-w>=xn(0) and x-w<=xn(Nx-1) and x-eps>=om(0) and x-eps<=om(No-1)) m_x.insert({x,i});
	  }
	  for (int i=0; i<xn.size(); ++i){
	    double x = xn(i)+w;
	    if (x>=xn(0) and x<=xn(Nx-1) and x-eps>=om(0) and x-eps<=om(No-1)) m_x.insert({x,-i-1}); // adding to m_eps -> xn(i)+w points
	  }
	  for (int i=0; i<om.size(); ++i){
	    double x = om(i)+eps;
	    if(x>=xn(0) and x<=xn(Nx-1) and x-w>=xn(0) and x-w<=xn(Nx-1)) m_x.insert({x,10000 + i});
	  }
	  enforce_min_spacing(m_x, 1e-10);
	  xc.resize(m_x.size());
	  xi.resize(m_x.size());
	  int jj=0;
	  for (const auto &[x,ix] : m_x){
	    xc(jj)=x;
	    xi(jj)=ix;
	    jj+=1;
	  }
	  int pos_x=0, pos_y=0, pos_z=0, pos_q=0, pos_o=om.size()-1;
	  Ri[0]=0; Ri[1]=0;
	  for (int jj=0; jj<xc.size(); ++jj){
	    double x = xc(jj);
	    int ix = xi(jj);
	    double dh=0; int N=xc.size();
	    if (jj>0 and jj<N-1){
	      dh = 0.5*(xc(jj+1)-xc(jj-1));
	    }else{
	      dh = (jj==0) ? 0.5*(xc(1)-xc(0)) : 0.5*(xc(N-1)-xc(N-2));
	    }
	    double f_x, f_x_w, f_x_eps;// f(x), f(x-w), f(x-eps)
	    if (ix>=10000){
	      int iix = ix-10000;
	      if (fabs(x-eps-om(iix))>1e-10){ cout<<"ERROR x-eps="<<x-eps<<" om["<<iix<<"]="<<om(iix)<<endl;}
	      f_x_eps = feo(iix);  // f(x-eps)
	    }else{
	      intpar py=InterpLeft(x-eps,pos_y,xn);
	      f_x_eps = fex(py);  // f(x-eps)
	    }
	    if (ix>=0 and ix<10000){
	      if (fabs(x-xn(ix))>1e-10){ cout<<"ERROR x="<<x<<" xn["<<ix<<"]="<<xn(ix)<<endl;}
	      f_x     = fex[ix];  // f(x)
	      for (int ip=0; ip<cix.Npseudo_unique; ++ip){
		Gx(ip) = Gm(ip,ix);                    // Gx(ip)  = Gm[ip](x)
		Gxt(ip) = Gx(ip).imag();               // Gxt(ip) = Gt[ip](x)
		Gx(ip) = complex<double>(Gx(ip).real(), -Gx(ip).imag()*(1-f_x)); // Gx(ip) = Gc[ip](x)
	      }
	    }else{
	      intpar px=InterpLeft(x,pos_x,xn);
	      f_x     = fex(px);  // f(x)
	      for (int ip=0; ip<cix.Npseudo_unique; ++ip){
		Gx(ip) = linInterp(Gm,ip,px);          // Gx(ip)  = Gm[ip](x)
		Gxt(ip) = Gx(ip).imag();               // Gxt(ip) = Gt[ip](x)
		Gx(ip) = complex<double>(Gx(ip).real(), -Gx(ip).imag()*(1-f_x)); // Gx(ip) = Gc[ip](x)
	      }
	    }
	    if (ix<0){
	      int iix=-ix-1; // xn(iix) is actually x-w!
	      if (fabs(x-w-xn(iix))>1e-10){ cout<<"ERROR x-w="<<x-w<<" while xn["<<iix<<"]="<<xn(iix)<<endl;}
	      f_x_w = fex[iix]; // f(x-w)
	      for (int ip=0; ip<cix.Npseudo_unique; ++ip){
		Gw(ip) = Gm(ip,iix);                  // Gw(ip)  = Gm[ip](x-w)
		Gwt(ip) = Gw(ip).imag();              // Gwt(ip) = Gt[ip](x-w)
		Gw(ip) = complex<double>(Gw(ip).real(), -Gw(ip).imag()*(1-f_x_w)); // Gw(ip) = Gc[ip](x-w)
	      }
	    }else{
	      intpar pz=InterpLeft(x-w,pos_z,xn);
	      f_x_w   = fex(pz);  // f(x-w)
	      for (int ip=0; ip<cix.Npseudo_unique; ++ip){
		Gw(ip) = linInterp(Gm,ip,pz);         // Gw(ip)  = Gm[ip](x-w)
		Gwt(ip) = Gw(ip).imag();              // Gwt(ip) = Gt[ip](x-w)
		Gw(ip) = complex<double>(Gw(ip).real(), -Gw(ip).imag()*(1-f_x_w)); // Gw(ip) = Gc[ip](x-w)
	      }
	    }
	    intpar pq=InterpLeft(x-eps,pos_q,om);
	    intpar po=InterpRight(eps-x,pos_o,om);
	    for (int is=0; is<cix.Nbath_unique; ++is){
	      Dip(is) = linInterp(Deltai,is,pq)*dh;          // D_1(x-eps)*dh
	      Dim(is) = linInterp(Deltai,is,po)*dh;          // D_1(eps-x)*dh
	    }
	    {
	      //               ^
	      //R[0](e,w)=     | w
	      //            --->----                
	      //states[3] /x      x-w\  states[0]
	      //         |            | 
	      //         |<--s0-- x-e |
	      //         |            |
	      //states[2] \e     e-w /  states[1]
	      //            --------
	      //               ^ w
	      //               |
	      //               s1
	      // 
	      //             ------>-                            		 
	      //states[3]  /e   |   x \  states[0]
	      //          |     s1     |
	      //    s0--> |     e-x    |->  R[1](e,w)
	      //          |     |      |
	      //states[2]  \e-w v     / x-w states[1]
	      //             ---|--<-			    		 
	      int ir=2;
	      for (int i=0; i<unique_states_inside[ir].size(); ++i){
		TinyVector<int,4> st = unique_states_inside[ir][i];  // st ~ (states[3],states[0],s0,1);
		double Ddh = Dip(st[2]); // either D(x-w) or D(w-x) depending on direction of hybridization. Needs to look at the diagram.
		Ri[0](i,0) += Ddh*f_x_eps          * Gx(st[0])       * Gw(st[1]).real();// D0(x-eps)*dh*f(x-eps)       * G[states[3]](x) * G'[states[0]](x-w)
		Ri[0](i,1) += Ddh*(1-f_x_eps)*f_x_w* Gx(st[0])       * Gwt(st[1]);      // D0(x-eps)*dh*f(eps-x)*f(x-w)* G[states[3]](x) * Gt[states[0]](x-w)
		Ri[0](i,2) += Ddh*f_x_eps          * Gx(st[0]).real()* Gw(st[1]);       // D0(x-eps)*dh*f(x-eps)       * G'[states[3]](x)* G[states[0]](x-w) 
		Ri[0](i,3) += Ddh*(1-f_x_eps)*f_x  * Gxt(st[0])      * Gw(st[1]);       // D0(x-eps)*dh*f(eps-x)*f(x)  * Gt[states[3]](x)* G[states[0]](x-w) 
	      }
	      ir=3;
	      for (int i=0; i<unique_states_inside[ir].size(); ++i){
		TinyVector<int,4> st = unique_states_inside[ir][i]; // st ~ (states[0],states[1],s1,1);
		double Ddh = Dim(st[2]); // either D(x-w) or D(w-x) depending on direction of hybridization. Needs to look at the diagram.
		Ri[1](i,0) += Ddh*f_x_eps          * Gx(st[0]).real() * Gw(st[1]);       // D1(eps-x)*dh*f(x-eps)       * G[st[0]](x) * G'[st[1]](x-w)
		Ri[1](i,1) += Ddh*(1-f_x_eps)*f_x  * Gxt(st[0])       * Gw(st[1]);       // D1(eps-x)*dh*f(eps-x)*f(x-w)* G[st[0]](x) * Gt[st[1]](x-w)
		Ri[1](i,2) += Ddh*f_x_eps          * Gx(st[0])        * Gw(st[1]).real();// D1(eps-x)*dh*f(x-eps)       * G'[st[0]](x)* G[st[1]](x-w)
		Ri[1](i,3) += Ddh*(1-f_x_eps)*f_x_w* Gx(st[0])        * Gwt(st[1]);      // D1(eps-x)*dh*f(eps-x)*f(x)  * Gt[st[0]](x)* G[st[1]](x-w)
	      }
	    }
	  }
	  //
	  double dhe=0; int Ne=eps_w.size();
	  if (iie>0 and iie<Ne-1)
	    dhe = 0.5*(eps_w(iie+1)-eps_w(iie-1));
	  else
	    dhe = (iie==0) ? 0.5*(eps_w(1)-eps_w(0)) : 0.5*(eps_w(Ne-1)-eps_w(Ne-2));
	  double f_eps = fermi(eps/T);     //f(eps)
	  double f_eps_w=fermi((eps-w)/T); //f(eps-w)
	  intpar pew = InterpLeft(eps-w,pos_w,xn);
	  intpar pee = InterpLeft(eps,pos_ex,xn);
	  for (int ip=0; ip<cix.Npseudo_unique; ++ip){
	    Gx(ip) = linInterp(Gm,ip,pee);         // Gx(ip)  = Gm[ip](eps)
	    Gxt(ip) = Gx(ip).imag();               // Gxt(ip) = Gt[ip](eps)
	    Gx(ip) = complex<double>(Gx(ip).real(), -Gx(ip).imag()*(1-f_eps)); // Gx(ip) = Gc[ip](eps)
	  }
	  for (int ip=0; ip<cix.Npseudo_unique; ++ip){
	    Gw(ip) = linInterp(Gm,ip,pew);        // Gw(ip)  = Gm[ip](eps-w)
	    Gwt(ip) = Gw(ip).imag();              // Gwt(ip) = Gt[ip](eps-w)
	    Gw(ip) = complex<double>(Gw(ip).real(), -Gw(ip).imag()*(1-f_eps_w)); // Gw(ip) = Gc[ip](eps-w)
	  }
	  complex<double> Gc0=0, Gc1=0;
	  for (auto const& [_st_, iind] : states_unique){
	    TinyVector<int,6> st = unpack6(_st_);
	    //           
	    //               ^ s1
	    //R[0](e,w)=     | w
	    //            --->----                
	    //states[3] /x      x-w\  states[0]
	    //         |            | 
	    //         |<--s0-- x-e |
	    //         |            |
	    //states[2] \e     e-w /  states[1]
	    //            ---@----
	    //              nds1
	    //               ^ w
	    //               |
	    //
	    //             ------>-                            		 
	    //states[3]  /e   |   x \  states[0]
	    //          |     s1     |
	    //s0-->nds0 @     e-x    |->  R[1](e,w)
	    //          |     |      |
	    //states[2]  \e-w v     / x-w states[1]
	    //             ---|--<-			    		 
	    double mm = -mm_nds[iind][0]; // minus comes from the crossing of hybridization lines
	    double nds1  = -mm_nds[iind][1];
	    double nds1_ = -mm_nds[iind][3];
	    double nds0  = -mm_nds[iind][2];
	    double nds0_ = -mm_nds[iind][4];
	    TinyVector<int,4> states(st[5],st[4],st[3],st[2]);
	    int s0=st[1], s1=st[0];
	    TinyVector<int,4> id(indexi[iind][2],indexi[iind][3]);
	    if (mm!=0 or nds1!=0 or nds1_!=0){
	      if (id[0]<0)
		cout<<"ERROR id[0]="<<id[0]<<" at st="<<st<<" iind="<<iind<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
	      complex<double> Gx2 = Gx(states[2]), Gw1 = Gw(states[1]);
	      double          Gt2 = Gxt(states[2]),Gwt1= Gwt(states[1]);
	      for (int ip=0; ip<4; ++ip) Rix(0,ip) =  Ri[0](id[0],ip);
	      Gc1  = dhe     *Gx2*(f_eps_w*Gwt1*Rix(0,0)+Gw1.real()*Rix(0,1)); // dh*  G[states[2]](eps) * ( f(eps-w)*Gt[states[1]](eps-w) * Ri[2,0] + Re(G[states[1]](eps-w)) * Ri[2,1] )
	      Gc1 -= dhe*conj(Gw1*(f_eps  *Gt2 *Rix(0,2)+Gx2.real()*Rix(0,3)));//-dh*( G[states[1]](eps-w) *( f(eps) *Gt[states[2]](eps)   * Ri[2,2] + Re(G[states[2]](eps))   * Ri[2,3] ) )^*
	      if (mm!=0) _Gc_(s1,iom) += Gc1*mm;
	      if (nds1 !=0) _Fc1_(s1,iom) += Gc1*nds1;
	      if (nds1_!=0) _Fc2_(s1,iom) += Gc1*nds1_;
	    }
	    if (mm!=0 or nds0!=0 or nds0_!=0){
	      if (id[1]<0)
		cout<<"ERROR id[1]="<<id[1]<<" at st="<<st<<" iind="<<iind<<" states="<<states<<" s0="<<s0<<" s1="<<s1<<endl;
	      complex<double> Gx3 = Gx(states[3]), Gw2 = Gw(states[2]);
	      double          Gt3 = Gxt(states[3]),Gwt2= Gwt(states[2]);
	      for (int ip=0; ip<4; ++ip) Rix(1,ip) =  Ri[1](id[1],ip);
	      Gc0  = dhe     *Gx3*(f_eps_w*Gwt2*Rix(1,2)+Gw2.real()*Rix(1,3)); // dh* G[states[3]](eps)  * ( f(eps-w)*Gt[states[2]](eps-w) * Ri[3,2] + Re(G[states[2]](eps-w)) * Ri[3,3] )
	      Gc0 -= dhe*conj(Gw2*(f_eps  *Gt3 *Rix(1,0)+Gx3.real()*Rix(1,1)));//-dh*(G[states[2]](eps-w)*(  f(eps)  *Gt[states[3]](eps)   * Ri[3,0] + Re(G[states[3]](eps))   * Ri[3,1] ) )^*
	      if (mm!=0) _Gc_(s0,iom) += Gc0*mm;
	      if (nds0 !=0) _Fc1_(s0,iom) += Gc0*nds0;
	      if (nds0_!=0) _Fc2_(s0,iom) += Gc0*nds0_;
	    }
	  }
	}
	if (mpi.rank==mpi.master){
	  auto t1 = high_resolution_clock::now();
	  fsec dt=t1-t0;
	  fmt::print(*mpi.out,"{}({:.4f}) ", iom, dt.count());
	  mpi.out->flush();
	}
      }
      if (mpi.rank==mpi.master) *mpi.out<<endl;
    
      if (mpi.size>1){
	mpi.AllReduceSum_simple_array(_Gc_);
	mpi.AllReduceSum_simple_array(_Fc1_);
	mpi.AllReduceSum_simple_array(_Fc2_);
      }
    }
  }
  if (cmpPhysG){
    Gc.resize(cix.Nbath_unique,om.size()); Gc=0;
    Fc1.resize(cix.Nbath_unique,om.size()); Fc1=0;
    Fc2.resize(cix.Nbath_unique,om.size()); Fc2=0;
    Array<int,1> sigind(cix.Nbath_unique);
    for (int iorb=0; iorb<cix.Nbath; ++iorb) sigind(cix.bath_sigind(iorb)) = cix.bath_sigind0(iorb);
    Array<double,1> deg_bth(cix.Nbath_unique_unique); deg_bth=0; // degeneracy of each bath
    for (int iorb=0; iorb<cix.Nbath; ++iorb) deg_bth(abs(cix.bath_sigind0(iorb)))+= 1;
    Array<bool,1> bath_SC(cix.Nbath_unique); bath_SC=false;
    for (int iorb=0; iorb<cix.Nbath; ++iorb) bath_SC(cix.bath_sigind(iorb)) = cix.SC_bath(iorb);
    for (int ib=0; ib<cix.Nbath_unique; ++ib){
      int s = abs(sigind(ib));
      double s_sgn = sigind(ib) < 0 ? -1 : 1;
      Gc(ib,Range::all()) = _Gc_(s,Range::all()) * s_sgn/(deg_bth(s)*Q*M_PI);
      double extra_2 = bath_SC(ib) ? 2.0 : 1.0;
      Fc1(ib,Range::all()) = _Fc1_(s,Range::all()) * (extra_2*s_sgn)/(deg_bth(s)*Q*M_PI);
      Fc2(ib,Range::all()) = _Fc2_(s,Range::all()) * (extra_2*s_sgn)/(deg_bth(s)*Q*M_PI);
    }
    
    if (debug){
      ofstream out1("GOCA_debug.dat");
      ofstream out2("FOCA1_debug.dat");
      ofstream out3("FOCA2_debug.dat");
      for (int iw=0; iw<om.size(); ++iw){
	out1 << om(iw) << " ";
	for (int ib=0; ib<Gc.extent(0); ++ib)
	  out1 << Gc(ib,iw).real() << " " << Gc(ib,iw).imag() << " ";
	out1<<endl;
	out2 << om(iw) << " ";
	for (int ib=0; ib<Fc1.extent(0); ++ib)
	  out2 << Fc1(ib,iw).real() <<" "<< Fc1(ib,iw).imag() << " ";
	out2 << endl;
	out3 << om(iw) << " ";
	for (int ib=0; ib<Fc2.extent(0); ++ib)
	  out3 << Fc2(ib,iw).real() <<" "<< Fc2(ib,iw).imag() << " ";
	out3 << endl;
      }
    }
  }
  return phi2;
}
