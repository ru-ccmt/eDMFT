#ifndef _frequency_loop_real_axis
#define _frequency_loop_real_axis

template <typename T>
int fsign(T x) {
    return (T(0) < x) - (x < T(0));
}
void DetermineLoops(Array<int,2>& wts, const Array<int,2>& diag, const Array<int,1>& dirspin)
{
  // wts(:, iloop) contains the argument for each propagator in loop `iloop`.
  //  - Values: +1, -1, or 0
  //    +1 → loop passes through propagator in the positive direction
  //    -1 → loop passes in the negative direction
  //     0 → not in the loop
  //
  //  Propagator ranges:
  //    - wts( [0, n),  iloop)     : hybridization propagators
  //    - wts( [n, 3n), iloop)     : PP propagators
  //
  //  Loop indices:
  //    - iloop = 0   : full impurity ring (frequency loop), goes through all PP propagators positively
  //    - iloop = i+1 : loop starting from hybridization i and traversing shortest path around impurity ring
  //  We choose to go through the impurity ring always in positive direction. Through hybridization we go either in positive or negative direction.
  //
  // wts[:,iloop] is either +1, -1, 0. If loop goes through propagator in positive (negative) direction it is +1 (-1), otherwise 0
  // wts[ range(0,n),   iloop] is reserved for hybridization propagators, since there are n hybridizations.
  // wts[ range(n,3*n), iloop] stands for PP propagators.
  // iloop=0 is frequency that does not involve any hybridizations, and goes around the impurity ring. It goes in positive direction
  //  hence wts[ range(n,3n), 0] = 1
  // The rest of loops wts[:,ihyb+1] go through hybridization first, and than go around the ring in positive direction, but such that the loop is the shortest.
  // We expect wts[ n+..,ihyb+1] = 1 or 0, for those PP that are (or are not) in the ihyb loop.
  //           wts[ ihyb,ihyb+1] = 1 or -1, depending on whether loop goes through hybridization in positive or negative direction.
  int n = diag.extent(0);
  wts.resize(3*n,n+1); // (3*n propagators, n+1 loops)
  wts=0;
  wts(Range(n,3*n-1),0)=1; // om loop
  for (int ip=0; ip<n; ++ip){ // hybridization loop ip
    int p0 = diag(ip,0);
    int p1 = diag(ip,1);
    int idir = fsign(p1-p0);
    if (idir<0) std::swap(p0,p1); // make sure p0<p1
    int loop_idx=ip+1;
    if (p1-p0<=n){ // shortest path is forward from p0 to p1
      for (int i=0; i<p1-p0; ++i) wts(p0+i +n,ip+1) = 1;
      wts(ip,loop_idx) = -idir * (1-2*(dirspin(ip)%2));
    }else{
      for (int i=0; i<p0+2*n-p1; ++i) wts( (p1+i)%(2*n) +n, ip+1)=1;
      wts(ip,loop_idx) = idir * (1-2*(dirspin(ip)%2));
    }
  }
}


void DetermineLoops(Array<int,2>& wts, const Array<int,2>& diag, int iremove=-1)
{
  // wts[:,iloop] is either +1, -1, 0. If loop goes through propagator in positive (negative) direction it is +1 (-1), otherwise 0
  // wts[ range(0,n),   iloop] is reserved for hybridization propagators, since there are n hybridizations.
  // wts[ range(n,3*n), iloop] stands for PP propagators.
  // iloop=0 is frequency that does not involve any hybridizations, and goes around the impurity ring. It goes in positive direction
  //  hence wts[ range(n,3n), 0] = 1
  // The rest of loops wts[:,ihyb+1] go through hybridization first, and than go around the ring, in the direction where the loops is the smallest.
  // 
  int n = diag.extent(0);
  int N=n;
  if (iremove>=0) N=n+1; // to have compatible wtb when removing one hybridization, we want to have wts larger.
  //if (wts.extent(0)!=3*N)
  wts.resize(N+2*n,N+1);
  wts=0;
  wts(Range(N,N+2*n-1),0)=1; // om loop
  for (int ip=0; ip<n; ++ip){
    int p0, p1;
    if (diag(ip,0)<diag(ip,1)){
      p0 = diag(ip,0);
      p1 = diag(ip,1);
    }else{
      p0 = diag(ip,1);
      p1 = diag(ip,0);
    }
    if (p1-p0<=n){
      int ipi = (iremove>=0 and ip>=iremove) ? ip+1 : ip;
      for (int i=0; i<p1-p0; ++i) wts(p0+i +N,ipi+1) = 1;
      wts(ipi,ipi+1) = diag(ip,1)>diag(ip,0) ? -1 : 1;
    }else{
      int ipi = (iremove>=0 and ip>=iremove) ? ip+1 : ip;
      for (int i=0; i<p0+2*n-p1; ++i) wts( (p1+i)%(2*n) +N, ipi+1)=1;
      wts(ipi,ipi+1) = diag(ip,1)>diag(ip,0) ? 1 : -1;
    }
  }
  //cout<<endl;
  //Print(cout,wts);
  //cout<<endl;
  //if (iremove>=0){
  //  wts(Range::all(),Range(iremove+2,n+1)) = wts(Range::all(),Range(iremove+1,n));
  //  wts(Range::all(),iremove+1)=0;
  //}
}
void CorrectHeader(Array<int,2>& wts, const Array<int,2>& diag){
  int n = diag.extent(0);
  for (int ip=0; ip<n; ++ip){
    int p0, p1;
    if (diag(ip,0)<diag(ip,1)){
      p0 = diag(ip,0);
      p1 = diag(ip,1);
    }else{
      p0 = diag(ip,1);
      p1 = diag(ip,0);
    }
    if (p1-p0<=n){
      wts(ip,ip+1) = diag(ip,1)>diag(ip,0) ? -1 : 1;
    }else{
      wts(ip,ip+1) = diag(ip,1)>diag(ip,0) ? 1 : -1;
    }
  }
}

class FreqKey{
public:
  int N_max_order, Npseudo_unique, Nbath;
  Array<int,1> to_p_3; // This is array for 3^n
  //
  FreqKey(int _N_max_order_, int _Npseudo_unique_, int _Nbath_) :
    N_max_order(_N_max_order_),  Npseudo_unique(_Npseudo_unique_), Nbath(_Nbath_), to_p_3(N_max_order+2)
  {
    //Array<int,1> to_p_3(N_max_order+2); // This is array for 3^n
    to_p_3(0)=1;
    for (int i=1; i<to_p_3.size(); ++i) to_p_3(i) = to_p_3(i-1)*3;
    //cout<<"3^i="<<to_p_3<<endl;
  }
  int ToKey(Array<int,2>& wtb, int i) const{
    int key=0;
    for (int j=0; j<wtb.extent(1); ++j){
      int ik = wtb(i,j);
      if (ik<0) ik=2;
      key += ik*to_p_3(j);
    }
    return key;
  }
  int ToKeyWithRemove(Array<int,2>& wtb, int i, int iremove) const{
    int key=0;
    for (int j=0; j<wtb.extent(1); ++j){
      int ik = wtb(i,j);
      if (ik<0) ik=2;
      int iy = j<=iremove ? j : j+1;
      key += ik*to_p_3(iy);
    }
    return key;
  }
  int AdjustKeyRemovingFreq(int iremove, int xkey) const{
    int iloop = iremove;
    if ( (xkey % to_p_3(iloop+1)) >= to_p_3(iloop) ){
      return -1; // We return negative key indicating that such key can not be properly modified as it contains this frequency. 
      // meaning that key corresponds to x = w + y_i + ... and we want to remove y_i from xkey.
    }else{
      // y_i can be safely removed from xkey, as y_i does not appear in frequency x
      int new_xkey = static_cast<int>(xkey / (to_p_3(iremove)*3)) * to_p_3(iremove);
      new_xkey += xkey % to_p_3(iremove);
      return new_xkey;
    }
  }
  int PPKey(int ipp, int key) const{
    return Npseudo_unique*key + ipp;
  }
  int DKey(int pm, int s, int iy) const{
    return pm + 2*s + (2*Nbath)*iy;
  }
  TinyVector<int,2> FromPPKey(int ppkey) const{
    int ipp = ppkey % Npseudo_unique;
    int key = ppkey/Npseudo_unique;
    return TinyVector<int,2>(ipp,key);
  }
  TinyVector<int,3> FromDKey(int ddkey) const{
    int pm = ddkey % 2;
    int rest = ddkey / 2;
    int s = rest % Nbath;
    int iy = rest/Nbath;
    return TinyVector<int,3>(pm,s,iy);
  }
  bool IsLoopInvolved(int iloop, int key) const{
    return (key % to_p_3(iloop+1)) >= to_p_3(iloop);
  }
  int ToLoop(int iloop, int xkey) const{
    // return -1, +1, or 0, depending on what is in xkey for this loop. Namely, xkey is build from loops and
    // stands for frequency, for example om+y_0-y_1+y_3 would stand for [1,1,-1,0,1,...] and this would return 1 for iloop=0,1,4 and -1 for 2, and 0 for 3...
    int ik = (xkey % to_p_3(iloop+1)) / to_p_3(iloop);
    if (ik>1) return -1;
    else return ik;
  }
  double ToFrDebug(int key, const Array<double,1>& freq, int Nloops) const{
    double x=0;
    int rest = key;
    cout<<"input key="<<key<<" freq="<<freq<<" Nloops="<<Nloops<<endl;
    for (int iloop=Nloops-1; iloop>=0; --iloop){
      int ik = rest/to_p_3(iloop);
      rest = rest % to_p_3(iloop);
      cout<<"iloop="<<iloop<<" ik="<<ik<<" rest="<<rest<<" ";
      if (ik==2) ik=-1;
      x += freq(iloop)*ik;
      cout<<"freq["<<iloop<<"]="<<freq(iloop)<<"*"<<ik<<" x="<<x<<endl;
    }
    return x;
  }
  double ToFr(int key, const Array<double,1>& freq, int Nloops) const{
    double x=0;
    int rest = key;
    for (int iloop=Nloops-1; iloop>=0; --iloop){
      int ik = rest/to_p_3(iloop);
      rest = rest % to_p_3(iloop);
      if (ik==2) ik=-1;
      x += freq(iloop)*ik;
    }
    return x;
  }
  Array<int,1> FromXKey(int key, int Nloops) const{
    Array<int,1> ss(Nloops);
    int rest = key;
    for (int iloop=Nloops-1; iloop>=0; --iloop){
      int ik = rest/to_p_3(iloop);
      rest = rest % to_p_3(iloop);
      ss(iloop) = (ik==2 ? -1 : ik);
    }
    return ss;
  }
};


Array<int,1> GivePMArray(int X, int n){
  Array<int,1> fsgn(n);
  for (int i=0; i<n; ++i){
    if (X & (1<<i)) fsgn(i)=1;
    else fsgn(i)=-1;
  }
  return fsgn;
}
Array<int,1> GiveZOArray(int X, int n){
  Array<int,1> fcng(n);
  for (int i=0; i<n; ++i){
    if (X & (1<<i)) fcng(i)=1;
    else fcng(i)=0;
  }
  return fcng;
}
deque<int> GiveSNArray(int X, int n){
  deque<int> creal;
  for (int i=0; i<=n; ++i){
    if (X & (1<<i)){
      creal.push_back(i);
    }
  }
  return creal;
}

void PrintM(ostream& out, const Array<int,2>& wts){
  for (int i=0; i<wts.extent(0); ++i){
    out<<"   ";
    for (int j=0; j<wts.extent(1); ++j)
      out<<setw(2)<<wts(i,j)<<", ";
    out<<endl;
  }
}

/*
void DetermineFermiSignsWithGS(Array<int,2>& fsgn, Array<bool,2>& ftype_conj, Array<deque<int>,2>& ftype_real,
			       Array<int,2>& wtb, const Array<int,1>& gs)
{
  // To determine the sign in fermi functions on the real axis, one needs to rearange the frequency loops so that
  // the Green'f function in the functional for this particular pseudo-particle m appears as the external frequency, i.e., G_m(omega),
  // while the other Green's functions will necessary contain more than one frequency argument G_n(omega+epsilon+...).
  // In this arrangement, the rest of the pseudoparticles will contain frequencies that go through hybridization lines. We choose each
  // hybridization line i to have its unique frequency epsilon_i. Than other pseudoparticles will contain G_n(omega+epsilon_1-epsilon_2+...)
  // It turns out that the only terms that survive projection have the same sign of epsilon_i in all pseudoparticle arguments (or they don't appear in some).
  // The fermi function has than the same sign as the arguments in the pseudo Green's function. For example, if epsilon_2 appears as G_n(omega-epsilon_2)*f(-epsilon_2),
  // the sign in the fermi function has to be negative. On the other hand, if G_n(omega+epsilon_2)*f(epsilon_2), the sign in the fermi function has to be positive.
  // 
  //  Example: consider diagram with hybridizations between [(0,5),(4,1),(2,6),(7,3)]
  //  hybridization from 0->5 is associated with y_0, 4-> with y_1, etc.
  //  We have the following frequency setup that can be used to evaluate
  //   Phi = Delta_0(y_0)*Delta_1(y_1)*Delta_2(y_2)*Delta_3(y_3)*G_0(om+y_2)*G_1(om+y_1+y_2)*G_2(om+y_1)*G_3(om+y_1+y_3)*G_4(om+y_3)*G_5(om+y_0+y_3)*G_6(om+y_0+y_2+y_3)*G_7(om+y_0+y_2)
  //  which can be presented with the following table:
  //            om, y_0,y_1,y_2,y_3             om, y_0,y_1,y_2,y_3          om, y_0,y_1,y_2,y_3
  //   Delta(0)  0   1   0   0   0               0   1   0   0   0	      0   1   0   0   0 
  //   Delta(1)  0   0   1   0   0	         0   0   1   0   0	      0   0   1   0   0 
  //   Delta(2)  0   0   0   1   0	         0   0   0   1   0	      0   0   0   1   0 
  //   Delta(3)  0   0   0   0   1	         0   0   0   0   1	      0   0   0   0   1 
  //   G_PP(0)   1   0   0   1   0	         1   0   0   0   0	      1   0  -1   0   0 
  //   G_PP(1)   1   0   1   1   0	         1   0   1   0   0	      1   0   0   0   0 
  //   G_PP(2)   1   0   1   0   0	         1   0   1  -1   0	      1   0   0  -1   0 
  //   G_PP(3)   1   0   1   0   1	         1   0   1  -1   1	      1   0   0  -1   1 
  //   G_PP(4)   1   0   0   0   1	         1   0   0  -1   1	      1   0  -1  -1   1 
  //   G_PP(5)   1   1   0   0   1	         1   1   0  -1   1	      1   1  -1  -1   1 
  //   G_PP(6)   1   1   0   1   1	         1   1   0   0   1	      1   1  -1   0   1 
  //   G_PP(7)   1   1   0   1   0 	         1   1   0   0   0 	      1   1  -1   0   0 
  //
  //  The two alternative setups on the right contain either G_0(om) and G_1(om), respectively. In the table with G_0(om) we see that y_2 always appears with negative sign in G, while
  //  y_0,y_1,y_3 appear with positive. It turns out the self-energy Sigma_0 will contain f(y_0)*f(y_1)*f(y_3)*f(-y_2).
  //  The last table on the right corresponds to Sigma_1, and contains y_0,y_3 with positive and y_1,y_2 with negative. The self-energy Sigma_1 will contain f(y_0)*f(-y_1)*f(-y_2)*f(y_3)
  //  In the code below, the input is wtb(3*norder,norder+1), which is the above described table (either one of them).
  //  From this table we construct 2^norder alternative tables, which are called wts, and contain only the right-lower corner of wts(2*norder,norder) parts, because the first column
  //  is always unchanged, and the top norder rows are also always unchanged.
  //  If we count sign in wts and check which G_m has no component in (2*norder,norder) part of the table (has only G_m(omega)), we can determine the sign of f(y_i) and know to which
  //  self-energy these signs correspond to.
  //  Although we always go through 2^norder options (because each y_i can be unchanged or added/subtracted om column), only 2*norder options make one row with zeros (except at om).
  //  It is precisely those 2*norder tables that correspond to PP-self-energy choice of frequency that determine signs of fermi functions.
  //wtb(3*norder,norder+1);
  //#define _DEBUG
  int norder = wtb.extent(1)-1;
  Array<int,1> is(norder); // is(i) is the sign of y_i variable in pseudoparticle green's functions. Note that all G_pp have the same sign of y_i variable.
  is=0;
  for (int ip=1; ip<=norder; ++ip){
    int sn=0, dn=0;
    for (int i=0; i<2*norder; ++i){
      if (wtb(norder+i,ip)!=0){
	sn += wtb(norder+i,ip);
	dn += 1;
      }
    }
    //cout<<"ip="<<ip<<" sn="<<sn<<" dn="<<dn<<endl;
    double s = (sn+0.0)/dn;
    if (fabs(fabs(s)-1)>1e-6){
      cout<<"ERROR at ip="<<ip<<" sign is not uniform "<<s<<endl;
      cout<<"wtb="<<wtb<<endl;
    }
    is(ip-1)=static_cast<int>(round(s));
  }
  //cout<<"s="<<is<<endl;
  //cout<<"wtb="<<wtb<<endl;
  Array<int,1> sb(norder); // sign for the baths. These are fixed and equal in wtb and all wts.
  for (int ib=0; ib<norder; ++ib) sb(ib) = wtb(ib,ib+1);
  int all_loops_length = sum(abs(wtb(Range(norder,3*norder-1),Range(1,norder))));
  bool want2change_wtb=true;
  bool succeeded=false;
#ifdef _DEBUG
  cout<<"sb="<<sb<<" all_loops_length="<<all_loops_length<<endl;
#endif  
  //Array<bool,2> ftype_conj(norder,2*norder);
  //Array<deque<int>,2> ftype_real(norder,2*norder);
  ftype_conj.resize(norder,2*norder);
  ftype_real.resize(norder,2*norder);
  ftype_conj = false;
  
  Array<int,1> ss(norder);
  //Array<int,2> fsgn(2*norder,norder); 
  fsgn.resize(2*norder,norder); 
  Array<int,2> wts(2*norder,norder), wts_stored(2*norder,norder);
  wts = 0;
  for (int i=0; i<(1<<norder); ++i){
#ifdef _DEBUG    
    cout<<i<<"   ";
#endif    
    for (int ip=0; ip<norder; ++ip){
      if (i & (1<<ip)){
#ifdef _DEBUG	
	cout<<-is(ip)<<" ";
#endif	
	wts(Range(0,2*norder-1),ip) = wtb(Range(norder,3*norder-1),ip+1)-is(ip);// wtb(Range(norder,3*norder-1),0);
	ss(ip)=-is(ip);
      }else{
#ifdef _DEBUG	
	cout<<0<<" ";
#endif	
	wts(Range(0,2*norder-1),ip) = wtb(Range(norder,3*norder-1),ip+1);
	ss(ip)=is(ip);
	//cout<<" wts="<<wts<<endl;
      }
    }
    // wts is constructed for one of 2^n combinations
    int pseudo=-1; // now checking if wts corresponds to one of tables that represent self-energy, i.e., G_pseudo(omega).
    for (int j=0; j<2*norder; ++j){
      if (sum(abs(wts(j,Range::all())))==0){
	pseudo=j; // yes, we have pseudo determined.
	break;
      }
    }
#ifdef _DEBUG    
    cout<<endl;
#endif    
    if (pseudo>=0){
      fsgn(pseudo,Range::all()) = ss(Range::all());
      for (int ib=0; ib<norder; ++ib)
	if (sb(ib)*ss(ib)<0) ftype_conj(ib,pseudo)=true;

      for (int ib=0; ib<norder; ++ib){
	for (int ip=0; ip<2*norder; ++ip)
	  if (ip!=pseudo and wts(ip,ib)==0) ftype_real(ib,pseudo).push_back(ip);
      }
      int my_loops_length = sum(abs(wts));
      if (i==0){
	if (gs(pseudo)==1){
#ifdef _DEBUG  	  
	  cout<<"We will keep wtb!"<<endl;
#endif	  
	  want2change_wtb=false;
	}else{
#ifdef _DEBUG  	  
	  cout<<"We want to replace wtb"<<endl;
#endif	  
	  want2change_wtb=true;
	}
      }else if(want2change_wtb){
	if (gs(pseudo)==1 and my_loops_length==all_loops_length){
#ifdef _DEBUG  
	  cout<<"This is competitor and wtb is changed!"<<endl;
#endif
	  //wtb(Range(norder,3*norder-1),Range(1,norder)) = wts(Range::all(),Range::all());
	  wts_stored = wts;
	  succeeded=true;
	}
      }
#ifdef _DEBUG      
      PrintM(cout,wts);
      cout<<"loops_length["<<pseudo<<"]="<<my_loops_length<<endl;
      cout<<"fsgn["<<pseudo<<"]="<<ss<<endl;
      for (int ib=0; ib<norder; ++ib)
	if (ftype_conj(ib,pseudo)) cout<<"ftype_conj["<<pseudo<<","<<ib<<"]="<<ftype_conj(ib,pseudo)<<endl;
      for (int ib=0; ib<norder; ++ib)
	if (ftype_real(ib,pseudo).size()>0)
	  cout<<"ftype_real["<<pseudo<<","<<ib<<"]="<< ftype_real(ib,pseudo)<<endl;
#endif      
    }
  }
  //cout<<"fsgn="<<fsgn<<endl;
  if(want2change_wtb and succeeded)
    wtb(Range(norder,3*norder-1),Range(1,norder)) = wts_stored(Range::all(),Range::all());
}
*/

void DetermineFermiSigns_fsgn_only(Array<int,1>& xfsgn, const Array<int,2>& wtb)
{
  // To determine the sign in fermi functions on the real axis, one needs to rearange the frequency loops so that
  // the Green'f function in the functional for this particular pseudo-particle m appears as the external frequency, i.e., G_m(omega),
  // while the other Green's functions will necessary contain more than one frequency argument G_n(omega+epsilon+...).
  // In this arrangement, the rest of the pseudoparticles will contain frequencies that go through hybridization lines. We choose each
  // hybridization line i to have its unique frequency epsilon_i. Than other pseudoparticles will contain G_n(omega+epsilon_1-epsilon_2+...)
  // It turns out that the only terms that survive projection have the same sign of epsilon_i in all pseudoparticle arguments (or they don't appear in some).
  // The fermi function has than the same sign as the arguments in the pseudo Green's function. For example, if epsilon_2 appears as G_n(omega-epsilon_2)*f(-epsilon_2),
  // the sign in the fermi function has to be negative. On the other hand, if G_n(omega+epsilon_2)*f(epsilon_2), the sign in the fermi function has to be positive.
  // 
  //  Example: consider diagram with hybridizations between [(0,5),(4,1),(2,6),(7,3)]
  //  hybridization from 0->5 is associated with y_0, 4-> with y_1, etc.
  //  We have the following frequency setup that can be used to evaluate
  //   Phi = Delta_0(y_0)*Delta_1(y_1)*Delta_2(y_2)*Delta_3(y_3)*G_0(om+y_2)*G_1(om+y_1+y_2)*G_2(om+y_1)*G_3(om+y_1+y_3)*G_4(om+y_3)*G_5(om+y_0+y_3)*G_6(om+y_0+y_2+y_3)*G_7(om+y_0+y_2)
  //  which can be presented with the following table:
  //            om, y_0,y_1,y_2,y_3             om, y_0,y_1,y_2,y_3          om, y_0,y_1,y_2,y_3
  //   Delta(0)  0   1   0   0   0               0   1   0   0   0	      0   1   0   0   0 
  //   Delta(1)  0   0   1   0   0	         0   0   1   0   0	      0   0   1   0   0 
  //   Delta(2)  0   0   0   1   0	         0   0   0   1   0	      0   0   0   1   0 
  //   Delta(3)  0   0   0   0   1	         0   0   0   0   1	      0   0   0   0   1 
  //   G_PP(0)   1   0   0   1   0	         1   0   0   0   0	      1   0  -1   0   0 
  //   G_PP(1)   1   0   1   1   0	         1   0   1   0   0	      1   0   0   0   0 
  //   G_PP(2)   1   0   1   0   0	         1   0   1  -1   0	      1   0   0  -1   0 
  //   G_PP(3)   1   0   1   0   1	         1   0   1  -1   1	      1   0   0  -1   1 
  //   G_PP(4)   1   0   0   0   1	         1   0   0  -1   1	      1   0  -1  -1   1 
  //   G_PP(5)   1   1   0   0   1	         1   1   0  -1   1	      1   1  -1  -1   1 
  //   G_PP(6)   1   1   0   1   1	         1   1   0   0   1	      1   1  -1   0   1 
  //   G_PP(7)   1   1   0   1   0 	         1   1   0   0   0 	      1   1  -1   0   0 
  //
  //  The two alternative setups on the right contain either G_0(om) and G_1(om), respectively. In the table with G_0(om) we see that y_2 always appears with negative sign in G, while
  //  y_0,y_1,y_3 appear with positive. It turns out the self-energy Sigma_0 will contain f(y_0)*f(y_1)*f(y_3)*f(-y_2).
  //  The last table on the right corresponds to Sigma_1, and contains y_0,y_3 with positive and y_1,y_2 with negative. The self-energy Sigma_1 will contain f(y_0)*f(-y_1)*f(-y_2)*f(y_3)
  //  In the code below, the input is wtb(3*norder,norder+1), which is the above described table (either one of them).
  //  From this table we construct 2^norder alternative tables, which are called wts, and contain only the right-lower corner of wts(2*norder,norder) parts, because the first column
  //  is always unchanged, and the top norder rows are also always unchanged.
  //  If we count sign in wts and check which G_m has no component in (2*norder,norder) part of the table (has only G_m(omega)), we can determine the sign of f(y_i) and know to which
  //  self-energy these signs correspond to.
  //  Although we always go through 2^norder options (because each y_i can be unchanged or added/subtracted om column), only 2*norder options make one row with zeros (except at om).
  //  It is precisely those 2*norder tables that correspond to PP-self-energy choice of frequency that determine signs of fermi functions.
  //#define _DEBUG
  int norder = wtb.extent(1)-1;
  Array<int,1> is(norder); // is(i) is the sign of y_i variable in pseudoparticle green's functions. Note that all G_pp have the same sign of y_i variable.
  is=0;
  for (int ip=1; ip<=norder; ++ip){
    int sn=0, dn=0;
    for (int i=0; i<2*norder; ++i){
      if (wtb(norder+i,ip)!=0){
	sn += wtb(norder+i,ip);
	dn += 1;
      }
    }
    //cout<<"ip="<<ip<<" sn="<<sn<<" dn="<<dn<<endl;
    double s = (sn+0.0)/dn;
    if (fabs(fabs(s)-1)>1e-6){
      cout<<"ERROR at ip="<<ip<<" sign is not uniform "<<s<<endl;
      cout<<"wtb="<<wtb<<endl;
    }
    is(ip-1)=static_cast<int>(round(s));
  }
  //Array<int,1> sb(norder); // sign for the baths. These are fixed and equal in wtb and all wts.
  //for (int ib=0; ib<norder; ++ib) sb(ib) = wtb(ib,ib+1);
#ifdef _DEBUG
  //cout<<"sb="<<sb<<endl;//" all_loops_length="<<all_loops_length<<endl;
#endif  

  xfsgn.resize(2*norder); xfsgn=0;
  //xftype_conj.resize(norder); xftype_conj=0;
  //xftype_real.resize(norder,2*norder); xftype_real=0;
  //
  Array<int,1> ss(norder);
  Array<int,2> wts(2*norder,norder);
  wts = 0;
  for (int i=0; i<(1<<norder); ++i){ // here we go over 2^n possible linear combinations, which try all possible loop combinations
#ifdef _DEBUG    
    cout<<i<<"   ";
#endif    
    for (int ip=0; ip<norder; ++ip){
      if (i & (1<<ip)){  // some selection of ip's will be selected to flip sign, and any number of ip's will be tried as i goes through all possible numbers up to 2^n
#ifdef _DEBUG	
	cout<<-is(ip)<<" ";
#endif	
	wts(Range(0,2*norder-1),ip) = wtb(Range(norder,3*norder-1),ip+1)-is(ip); // this flips one column, namely, wtb(:,ip+1) so that all nonzeros become zeros and all zeros become nonzeros.
	                   // it is equivalent to reversing the loop: if it was going in the clockwise direction, it will now go in counterclockwise direction.
	ss(ip)=-is(ip);
      }else{
#ifdef _DEBUG	
	cout<<0<<" ";
#endif	
	wts(Range(0,2*norder-1),ip) = wtb(Range(norder,3*norder-1),ip+1); // othersie just copy the pseudo-part of wtb table. We always skip the header (bath part).
	ss(ip)=is(ip);
	//cout<<" wts="<<wts<<endl;
      }
    }
    // wts is constructed for one of 2^n combinations
    int pseudo=-1; // now checking if wts corresponds to one of tables that represent self-energy, i.e., G_pseudo(omega).
    // It means that pseudo line wts(ip,:) should contain only omega, i.e., should be 0,0,0,0 because wtb=1,0,0,0,0
    for (int j=0; j<2*norder; ++j){
      if (sum(abs(wts(j,Range::all())))==0){
	pseudo=j; // yes, we have pseudo determined.
	break;
      }
    }
#ifdef _DEBUG    
    cout<<endl;
#endif    
    if (pseudo>=0){
      // (ss+1)/2
      //-1 -> 0
      // 1 -> 1
      //fsgn(pseudo,Range::all()) = ss(Range::all());
      //for (int ib=0; ib<norder; ++ib)
      //  if (sb(ib)*ss(ib)<0) ftype_conj(ib,pseudo)=true;
      //for (int ib=0; ib<norder; ++ib)
      // 	for (int ip=0; ip<2*norder; ++ip)
      // 	  if (ip!=pseudo and wts(ip,ib)==0) ftype_real(ib,pseudo).push_back(ip);
      int fsg=0; // this will pack array ss composed of 0 & 1 into its binary representation. Instead of 0,1,0,0,1,0 we will get number 10.
      for (int l=0; l<norder; ++l)
	if (ss(l)==1) fsg += (1<<l);
      xfsgn(pseudo) = fsg;
    }
  }
  //cout<<"fsgn="<<fsgn<<endl;
}
void DetermineFermiSigns(Array<int,1>& xfsgn, Array<int,1>& xftype_conj, Array<int,2>& xftype_real, const Array<int,2>& wtb)
{
  // To determine the sign in fermi functions on the real axis, one needs to rearange the frequency loops so that
  // the Green'f function in the functional for this particular pseudo-particle m appears as the external frequency, i.e., G_m(omega),
  // while the other Green's functions will necessary contain more than one frequency argument G_n(omega+epsilon+...).
  // In this arrangement, the rest of the pseudoparticles will contain frequencies that go through hybridization lines. We choose each
  // hybridization line i to have its unique frequency epsilon_i. Than other pseudoparticles will contain G_n(omega+epsilon_1-epsilon_2+...)
  // It turns out that the only terms that survive projection have the same sign of epsilon_i in all pseudoparticle arguments (or they don't appear in some).
  // The fermi function has than the same sign as the arguments in the pseudo Green's function. For example, if epsilon_2 appears as G_n(omega-epsilon_2)*f(-epsilon_2),
  // the sign in the fermi function has to be negative. On the other hand, if G_n(omega+epsilon_2)*f(epsilon_2), the sign in the fermi function has to be positive.
  // 
  //  Example: consider diagram with hybridizations between [(0,5),(4,1),(2,6),(7,3)]
  //  hybridization from 0->5 is associated with y_0, 4-> with y_1, etc.
  //  We have the following frequency setup that can be used to evaluate
  //   Phi = Delta_0(y_0)*Delta_1(y_1)*Delta_2(y_2)*Delta_3(y_3)*G_0(om+y_2)*G_1(om+y_1+y_2)*G_2(om+y_1)*G_3(om+y_1+y_3)*G_4(om+y_3)*G_5(om+y_0+y_3)*G_6(om+y_0+y_2+y_3)*G_7(om+y_0+y_2)
  //  which can be presented with the following table:
  //            om, y_0,y_1,y_2,y_3             om, y_0,y_1,y_2,y_3          om, y_0,y_1,y_2,y_3
  //   Delta(0)  0   1   0   0   0               0   1   0   0   0	      0   1   0   0   0 
  //   Delta(1)  0   0   1   0   0	         0   0   1   0   0	      0   0   1   0   0 
  //   Delta(2)  0   0   0   1   0	         0   0   0   1   0	      0   0   0   1   0 
  //   Delta(3)  0   0   0   0   1	         0   0   0   0   1	      0   0   0   0   1 
  //   G_PP(0)   1   0   0   1   0	         1   0   0   0   0	      1   0  -1   0   0 
  //   G_PP(1)   1   0   1   1   0	         1   0   1   0   0	      1   0   0   0   0 
  //   G_PP(2)   1   0   1   0   0	         1   0   1  -1   0	      1   0   0  -1   0 
  //   G_PP(3)   1   0   1   0   1	         1   0   1  -1   1	      1   0   0  -1   1 
  //   G_PP(4)   1   0   0   0   1	         1   0   0  -1   1	      1   0  -1  -1   1 
  //   G_PP(5)   1   1   0   0   1	         1   1   0  -1   1	      1   1  -1  -1   1 
  //   G_PP(6)   1   1   0   1   1	         1   1   0   0   1	      1   1  -1   0   1 
  //   G_PP(7)   1   1   0   1   0 	         1   1   0   0   0 	      1   1  -1   0   0 
  //
  //  The two alternative setups on the right contain either G_0(om) and G_1(om), respectively. In the table with G_0(om) we see that y_2 always appears with negative sign in G, while
  //  y_0,y_1,y_3 appear with positive. It turns out the self-energy Sigma_0 will contain f(y_0)*f(y_1)*f(y_3)*f(-y_2).
  //  The last table on the right corresponds to Sigma_1, and contains y_0,y_3 with positive and y_1,y_2 with negative. The self-energy Sigma_1 will contain f(y_0)*f(-y_1)*f(-y_2)*f(y_3)
  //  In the code below, the input is wtb(3*norder,norder+1), which is the above described table (either one of them).
  //  From this table we construct 2^norder alternative tables, which are called wts, and contain only the right-lower corner of wts(2*norder,norder) parts, because the first column
  //  is always unchanged, and the top norder rows are also always unchanged.
  //  If we count sign in wts and check which G_m has no component in (2*norder,norder) part of the table (has only G_m(omega)), we can determine the sign of f(y_i) and know to which
  //  self-energy these signs correspond to.
  //  Although we always go through 2^norder options (because each y_i can be unchanged or added/subtracted om column), only 2*norder options make one row with zeros (except at om).
  //  It is precisely those 2*norder tables that correspond to PP-self-energy choice of frequency that determine signs of fermi functions.
  //wtb(3*norder,norder+1);
  //#define _DEBUG
  int norder = wtb.extent(1)-1;
  Array<int,1> is(norder); // is(i) is the sign of y_i variable in pseudoparticle green's functions. Note that all G_pp have the same sign of y_i variable.
  is=0;
  for (int ip=1; ip<=norder; ++ip){
    int sn=0, dn=0;
    for (int i=0; i<2*norder; ++i){
      if (wtb(norder+i,ip)!=0){
	sn += wtb(norder+i,ip);
	dn += 1;
      }
    }
    //cout<<"ip="<<ip<<" sn="<<sn<<" dn="<<dn<<endl;
    double s = (sn+0.0)/dn;
    if (fabs(fabs(s)-1)>1e-6){
      cout<<"ERROR at ip="<<ip<<" sign is not uniform "<<s<<endl;
      cout<<"wtb="<<wtb<<endl;
    }
    is(ip-1)=static_cast<int>(round(s));
  }
  //cout<<"s="<<is<<endl;
  //cout<<"wtb="<<wtb<<endl;
  Array<int,1> sb(norder); // sign for the baths. These are fixed and equal in wtb and all wts.
  for (int ib=0; ib<norder; ++ib) sb(ib) = wtb(ib,ib+1);
  //int all_loops_length = sum(abs(wtb(Range(norder,3*norder-1),Range(1,norder))));
#ifdef _DEBUG
  cout<<"sb="<<sb<<endl;//" all_loops_length="<<all_loops_length<<endl;
#endif  
  //ftype_conj.resize(norder,2*norder);
  //ftype_real.resize(norder,2*norder);
  //ftype_conj = false;
  //for (int i=0; i<norder; ++i)
  //  for (int j=0; j<2*norder; ++j)
  //    ftype_real(i,j).clear();
  //fsgn.resize(2*norder,norder); 

  xfsgn.resize(2*norder); xfsgn=0;
  xftype_conj.resize(norder); xftype_conj=0;
  xftype_real.resize(norder,2*norder); xftype_real=0;
  //
  Array<int,1> ss(norder);
  Array<int,2> wts(2*norder,norder);
  wts = 0;
  for (int i=0; i<(1<<norder); ++i){ // here we go over 2^n possible linear combinations, which try all possible loop combinations
#ifdef _DEBUG    
    cout<<i<<"   ";
#endif    
    for (int ip=0; ip<norder; ++ip){
      if (i & (1<<ip)){  // some selection of ip's will be selected to flip sign, and any number of ip's will be tried as i goes through all possible numbers up to 2^n
#ifdef _DEBUG	
	cout<<-is(ip)<<" ";
#endif	
	wts(Range(0,2*norder-1),ip) = wtb(Range(norder,3*norder-1),ip+1)-is(ip); // this flips one column, namely, wtb(:,ip+1) so that all nonzeros become zeros and all zeros become nonzeros.
	                   // it is equivalent to reversing the loop: if it was going in the clockwise direction, it will now go in counterclockwise direction.
	ss(ip)=-is(ip);
      }else{
#ifdef _DEBUG	
	cout<<0<<" ";
#endif	
	wts(Range(0,2*norder-1),ip) = wtb(Range(norder,3*norder-1),ip+1); // othersie just copy the pseudo-part of wtb table. We always skip the header (bath part).
	ss(ip)=is(ip);
	//cout<<" wts="<<wts<<endl;
      }
    }
    // wts is constructed for one of 2^n combinations
    int pseudo=-1; // now checking if wts corresponds to one of tables that represent self-energy, i.e., G_pseudo(omega).
    // It means that pseudo line wts(ip,:) should contain only omega, i.e., should be 0,0,0,0 because wtb=1,0,0,0,0
    for (int j=0; j<2*norder; ++j){
      if (sum(abs(wts(j,Range::all())))==0){
	pseudo=j; // yes, we have pseudo determined.
	break;
      }
    }
#ifdef _DEBUG    
    cout<<endl;
#endif    
    if (pseudo>=0){
      // (ss+1)/2
      //-1 -> 0
      // 1 -> 1
      //fsgn(pseudo,Range::all()) = ss(Range::all());
      //for (int ib=0; ib<norder; ++ib)
      //  if (sb(ib)*ss(ib)<0) ftype_conj(ib,pseudo)=true;
      //for (int ib=0; ib<norder; ++ib)
      // 	for (int ip=0; ip<2*norder; ++ip)
      // 	  if (ip!=pseudo and wts(ip,ib)==0) ftype_real(ib,pseudo).push_back(ip);
      
      int fsg=0; // this will pack array ss composed of 0 & 1 into its binary representation. Instead of 0,1,0,0,1,0 we will get number 10.
      for (int l=0; l<norder; ++l)
	if (ss(l)==1) fsg += (1<<l);
      xfsgn(pseudo) = fsg;
      // this will pack conjugation into binary representation. if sb(ib)*ss(ib)<0, which means that column for pseudo and column for bath have opposite sign, we need to
      // conjugate. This means that pseudo would have arguments Delta(y_i)*G(om - y_i...) or Delta(-y_i)*G(om + y_i...)
      for (int ib=0; ib<norder; ++ib)
	if (sb(ib)*ss(ib)<0)
	  xftype_conj(ib) += (1<<pseudo);
      // finally we determine which terms should have real part of G rathe than complex G.
      // If G (apart from the one we are eliminating by functional derivatibe ==pseudo)
      //  does not contain y_i in its argument, i.e., Delta(y_i)*G(om... +0*y_i...), it should appear real.
      for (int ib=0; ib<norder; ++ib){
	int reals=0;
	for (int ip=0; ip<2*norder; ++ip)
	  if (ip!=pseudo and wts(ip,ib)==0) reals += (1<<ip); // this ip is packed into reals.
	xftype_real(ib,pseudo) = reals;
      }
      
      //int my_loops_length = sum(abs(wts));
#ifdef _DEBUG      
      PrintM(cout,wts);
      //cout<<"loops_length["<<pseudo<<"]="<<my_loops_length<<endl;
      cout<<"fsgn["<<pseudo<<"]="<<ss<<endl;
      //for (int ib=0; ib<norder; ++ib)
      //  if (ftype_conj(ib,pseudo)) cout<<"ftype_conj["<<pseudo<<","<<ib<<"]="<<ftype_conj(ib,pseudo)<<endl;
      //for (int ib=0; ib<norder; ++ib)
      //  cout<<"ftype_real["<<pseudo<<","<<ib<<"]="<< ftype_real(ib,pseudo)<<endl;
#endif      
    }
  }
  //cout<<"fsgn="<<fsgn<<endl;
}
/*
void DetermineFermiSigns(Array<int,2>& fsgn, Array<bool,2>& ftype_conj, Array<deque<int>,2>& ftype_real, const Array<int,2>& wtb,
			 Array<int,1>& xfsgn, Array<int,1>& xftype_conj, Array<int,2>& xftype_real)
{
  // To determine the sign in fermi functions on the real axis, one needs to rearange the frequency loops so that
  // the Green'f function in the functional for this particular pseudo-particle m appears as the external frequency, i.e., G_m(omega),
  // while the other Green's functions will necessary contain more than one frequency argument G_n(omega+epsilon+...).
  // In this arrangement, the rest of the pseudoparticles will contain frequencies that go through hybridization lines. We choose each
  // hybridization line i to have its unique frequency epsilon_i. Than other pseudoparticles will contain G_n(omega+epsilon_1-epsilon_2+...)
  // It turns out that the only terms that survive projection have the same sign of epsilon_i in all pseudoparticle arguments (or they don't appear in some).
  // The fermi function has than the same sign as the arguments in the pseudo Green's function. For example, if epsilon_2 appears as G_n(omega-epsilon_2)*f(-epsilon_2),
  // the sign in the fermi function has to be negative. On the other hand, if G_n(omega+epsilon_2)*f(epsilon_2), the sign in the fermi function has to be positive.
  // 
  //  Example: consider diagram with hybridizations between [(0,5),(4,1),(2,6),(7,3)]
  //  hybridization from 0->5 is associated with y_0, 4-> with y_1, etc.
  //  We have the following frequency setup that can be used to evaluate
  //   Phi = Delta_0(y_0)*Delta_1(y_1)*Delta_2(y_2)*Delta_3(y_3)*G_0(om+y_2)*G_1(om+y_1+y_2)*G_2(om+y_1)*G_3(om+y_1+y_3)*G_4(om+y_3)*G_5(om+y_0+y_3)*G_6(om+y_0+y_2+y_3)*G_7(om+y_0+y_2)
  //  which can be presented with the following table:
  //            om, y_0,y_1,y_2,y_3             om, y_0,y_1,y_2,y_3          om, y_0,y_1,y_2,y_3
  //   Delta(0)  0   1   0   0   0               0   1   0   0   0	      0   1   0   0   0 
  //   Delta(1)  0   0   1   0   0	         0   0   1   0   0	      0   0   1   0   0 
  //   Delta(2)  0   0   0   1   0	         0   0   0   1   0	      0   0   0   1   0 
  //   Delta(3)  0   0   0   0   1	         0   0   0   0   1	      0   0   0   0   1 
  //   G_PP(0)   1   0   0   1   0	         1   0   0   0   0	      1   0  -1   0   0 
  //   G_PP(1)   1   0   1   1   0	         1   0   1   0   0	      1   0   0   0   0 
  //   G_PP(2)   1   0   1   0   0	         1   0   1  -1   0	      1   0   0  -1   0 
  //   G_PP(3)   1   0   1   0   1	         1   0   1  -1   1	      1   0   0  -1   1 
  //   G_PP(4)   1   0   0   0   1	         1   0   0  -1   1	      1   0  -1  -1   1 
  //   G_PP(5)   1   1   0   0   1	         1   1   0  -1   1	      1   1  -1  -1   1 
  //   G_PP(6)   1   1   0   1   1	         1   1   0   0   1	      1   1  -1   0   1 
  //   G_PP(7)   1   1   0   1   0 	         1   1   0   0   0 	      1   1  -1   0   0 
  //
  //  The two alternative setups on the right contain either G_0(om) and G_1(om), respectively. In the table with G_0(om) we see that y_2 always appears with negative sign in G, while
  //  y_0,y_1,y_3 appear with positive. It turns out the self-energy Sigma_0 will contain f(y_0)*f(y_1)*f(y_3)*f(-y_2).
  //  The last table on the right corresponds to Sigma_1, and contains y_0,y_3 with positive and y_1,y_2 with negative. The self-energy Sigma_1 will contain f(y_0)*f(-y_1)*f(-y_2)*f(y_3)
  //  In the code below, the input is wtb(3*norder,norder+1), which is the above described table (either one of them).
  //  From this table we construct 2^norder alternative tables, which are called wts, and contain only the right-lower corner of wts(2*norder,norder) parts, because the first column
  //  is always unchanged, and the top norder rows are also always unchanged.
  //  If we count sign in wts and check which G_m has no component in (2*norder,norder) part of the table (has only G_m(omega)), we can determine the sign of f(y_i) and know to which
  //  self-energy these signs correspond to.
  //  Although we always go through 2^norder options (because each y_i can be unchanged or added/subtracted om column), only 2*norder options make one row with zeros (except at om).
  //  It is precisely those 2*norder tables that correspond to PP-self-energy choice of frequency that determine signs of fermi functions.
  //wtb(3*norder,norder+1);
  //#define _DEBUG
  int norder = wtb.extent(1)-1;
  Array<int,1> is(norder); // is(i) is the sign of y_i variable in pseudoparticle green's functions. Note that all G_pp have the same sign of y_i variable.
  is=0;
  for (int ip=1; ip<=norder; ++ip){
    int sn=0, dn=0;
    for (int i=0; i<2*norder; ++i){
      if (wtb(norder+i,ip)!=0){
	sn += wtb(norder+i,ip);
	dn += 1;
      }
    }
    //cout<<"ip="<<ip<<" sn="<<sn<<" dn="<<dn<<endl;
    double s = (sn+0.0)/dn;
    if (fabs(fabs(s)-1)>1e-6){
      cout<<"ERROR at ip="<<ip<<" sign is not uniform "<<s<<endl;
      cout<<"wtb="<<wtb<<endl;
    }
    is(ip-1)=static_cast<int>(round(s));
  }
  //cout<<"s="<<is<<endl;
  //cout<<"wtb="<<wtb<<endl;
  Array<int,1> sb(norder); // sign for the baths. These are fixed and equal in wtb and all wts.
  for (int ib=0; ib<norder; ++ib) sb(ib) = wtb(ib,ib+1);
  int all_loops_length = sum(abs(wtb(Range(norder,3*norder-1),Range(1,norder))));
#ifdef _DEBUG
  cout<<"sb="<<sb<<" all_loops_length="<<all_loops_length<<endl;
#endif  
  //Array<bool,2> ftype_conj(norder,2*norder);
  //Array<deque<int>,2> ftype_real(norder,2*norder);
  ftype_conj.resize(norder,2*norder);
  ftype_real.resize(norder,2*norder);
  ftype_conj = false;

  for (int i=0; i<norder; ++i)
    for (int j=0; j<2*norder; ++j)
      ftype_real(i,j).clear();
    
  xfsgn.resize(2*norder); xfsgn=0;
  xftype_conj.resize(norder); xftype_conj=0;
  xftype_real.resize(norder,2*norder); xftype_real=0;
  
  
  Array<int,1> ss(norder);
  //Array<int,2> fsgn(2*norder,norder); 
  fsgn.resize(2*norder,norder); 
  Array<int,2> wts(2*norder,norder);
  wts = 0;
  for (int i=0; i<(1<<norder); ++i){
#ifdef _DEBUG    
    cout<<i<<"   ";
#endif    
    for (int ip=0; ip<norder; ++ip){
      if (i & (1<<ip)){
#ifdef _DEBUG	
	cout<<-is(ip)<<" ";
#endif	
	wts(Range(0,2*norder-1),ip) = wtb(Range(norder,3*norder-1),ip+1)-is(ip);// wtb(Range(norder,3*norder-1),0);
	ss(ip)=-is(ip);
      }else{
#ifdef _DEBUG	
	cout<<0<<" ";
#endif	
	wts(Range(0,2*norder-1),ip) = wtb(Range(norder,3*norder-1),ip+1);
	ss(ip)=is(ip);
	//cout<<" wts="<<wts<<endl;
      }
    }
    // wts is constructed for one of 2^n combinations
    int pseudo=-1; // now checking if wts corresponds to one of tables that represent self-energy, i.e., G_pseudo(omega).
    for (int j=0; j<2*norder; ++j){
      if (sum(abs(wts(j,Range::all())))==0){
	pseudo=j; // yes, we have pseudo determined.
	break;
      }
    }
#ifdef _DEBUG    
    cout<<endl;
#endif    
    if (pseudo>=0){
      // (ss+1)/2
      //-1 -> 0
      // 1 -> 1
      fsgn(pseudo,Range::all()) = ss(Range::all());
      
      for (int ib=0; ib<norder; ++ib)
        if (sb(ib)*ss(ib)<0) ftype_conj(ib,pseudo)=true;
      
      for (int ib=0; ib<norder; ++ib)
      	for (int ip=0; ip<2*norder; ++ip)
      	  if (ip!=pseudo and wts(ip,ib)==0) ftype_real(ib,pseudo).push_back(ip);
      
      int fsg=0;
      for (int l=0; l<norder; ++l)
	if (ss(l)==1) fsg += (1<<l);
      xfsgn(pseudo) = fsg;
      
      for (int ib=0; ib<norder; ++ib)
	if (sb(ib)*ss(ib)<0)
	  xftype_conj(ib) += (1<<pseudo);
      
      for (int ib=0; ib<norder; ++ib){
	int reals=0;
	for (int ip=0; ip<2*norder; ++ip)
	  if (ip!=pseudo and wts(ip,ib)==0) reals += (1<<ip);
	xftype_real(ib,pseudo) = reals;
      }
      
      //int my_loops_length = sum(abs(wts));
#ifdef _DEBUG      
      PrintM(cout,wts);
      //cout<<"loops_length["<<pseudo<<"]="<<my_loops_length<<endl;
      cout<<"fsgn["<<pseudo<<"]="<<ss<<endl;
      for (int ib=0; ib<norder; ++ib)
	if (ftype_conj(ib,pseudo)) cout<<"ftype_conj["<<pseudo<<","<<ib<<"]="<<ftype_conj(ib,pseudo)<<endl;
      for (int ib=0; ib<norder; ++ib)
	cout<<"ftype_real["<<pseudo<<","<<ib<<"]="<< ftype_real(ib,pseudo)<<endl;
#endif      
    }
  }
  //cout<<"fsgn="<<fsgn<<endl;
}
*/
#endif // _frequency_loop_real_axis
