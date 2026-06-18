#ifndef _SKELETONS
#define _SKELETONS

class Skeletons{
  // Contains information to quickly find if diagram is skeleton.
  //  -- momenta contain auxiliary momentum of every hybridization line. If two pseudo propagators have the same auxiliary momentum it is non-skeleton.
  //  -- loops contain a list of momentum loops for the current diagram.
public:
  static Array<int,1> Freq;
  vector<deque<int>> loops;
  Array<int,1> momenta;
  //bool debug;
  ////
  Skeletons(){};
  Skeletons(const Array<int,2>& starting_diag) : loops(starting_diag.extent(0)), momenta(2*starting_diag.extent(0))//, debug(false)
  {
    DetermineLoops(starting_diag);
    //PrintLoops();
    InitializeMomenta();
  }
  void Resize_Initialize(const Array<int,2>& starting_diag){//, bool _debug_=false){
    //debug = _debug_;
    int Norder = starting_diag.extent(0);
    loops.resize(Norder);
    momenta.resize(2*Norder);
    DetermineLoops(starting_diag);
    //if (debug) PrintLoops();
    InitializeMomenta();
  }
  Skeletons& operator=(const Skeletons& skl)
  {
    //cout<<"Copying skeletons from "<< (&skl) <<" to *this="<< this<<" sizes="<<momenta.size()<<","<<skl.momenta.size()<<","<<loops.size()<<","<<skl.loops.size()<<endl;
    //cout<<"skl.loops[2].size()="<<skl.loops[2].size() <<" *="<< (&skl.loops[2][0])<<endl;
    if (momenta.size()!=skl.momenta.size()) momenta.resize(skl.momenta.size());
    //cout<<"momenta="<<momenta<<" *="<< (&momenta) << endl;
    //cout<<"skl.momenta="<<skl.momenta<<" *="<< (&skl.momenta) << endl;
    momenta = skl.momenta.copy();
    if (loops.size()!=skl.loops.size())  loops.resize(skl.loops.size());
    //cout<<"loops="<<loops<<" *="<< (&loops) << endl;
    //cout<<"skl.loops="<<skl.loops<<" *="<< (&skl.loops) << endl;
    //cout<<"skl.loops[2].size()="<<skl.loops[2].size() <<endl;
    for (int i=0; i<loops.size(); ++i){
      //cout<<i<<" "<<skl.loops[i].size()<<" loops[i]="<< &(loops[i][0]) <<" skl.loops[i]="<<(&skl.loops[i][0])<<endl;
      //loops[i].resize(skl.loops[i].size());
      loops[i] = skl.loops[i];
    }
    //cout<<"at the end loops="<<loops<<" *="<< (&loops) << endl;
    return *this;
  }
  void copy(const Skeletons& skl)
  {
    if (momenta.size()!=skl.momenta.size()) momenta.resize(skl.momenta.size());
    momenta = skl.momenta.copy();
    if (loops.size()!=skl.loops.size())  loops.resize(skl.loops.size());
    for (int i=0; i<loops.size(); ++i) loops[i] = skl.loops[i];
  }
  bool T_Is_Reducible(const Array<int,2>& diag)
  {
    int n = diag.extent(0);
    Array<int,1> momenta(2*n);
    momenta=0;
    //cout<<"Handling diag=" << diag << endl;
    for (int ip=0; ip<n; ++ip){
      int w = Freq(ip);
      int p0 = diag(ip,0);
      int p1 = diag(ip,1);
      if (p0>p1) swap(p0,p1);
      // Now p0 < p1
      if ( abs(p1-p0)<=n ){
	int m = p1-p0+1;
	for (int i=0; i<m-1; ++i){
	  int loops_ip_i = p1-i;
	  int loops_ip_i_p_1 = p1-i-1;
	  //cout << "Adding to k["<<loops_ip_i_p_1 << "]+=w["<<ip<<"] w="<<w<<endl;
	  momenta(loops_ip_i_p_1) += w;
	}
      }else{
	int m = p0-p1+2*n+1;
	for (int i=0; i<m-1; ++i){
	  int loops_ip_i = (p0+2*n-i)%(2*n);
	  int loops_ip_i_p_1 = (p0+2*n-i-1)%(2*n);
	  //cout << "Adding to k["<<loops_ip_i_p_1 << "]+=w["<<ip<<"] w="<<w<<endl;
	  momenta(loops_ip_i_p_1) += w;
	}
      }
    }
    //cout<<"momenta="<<momenta<<endl;
    for (int i=0; i<momenta.extent(0); ++i)
      for (int j=i+1; j<momenta.extent(0); ++j)
	if (momenta(i)==momenta(j)) return true;
    return false;
  }
  void DetermineLoops(const Array<int,2>& diag)
  {
    int n = diag.extent(0);
    for (int ip=0; ip<n; ++ip){
      int p0 = diag(ip,0);
      int p1 = diag(ip,1);
      if ( abs(p1-p0)<=n ){
	if (p1>p0){
	  loops[ip].resize(p1-p0+1);
	  for (int i=0; i<p1-p0+1; ++i) loops[ip][i] = p1-i;
	}else{
	  loops[ip].resize(p0-p1+1);
	  for (int i=0; i<p0-p1+1; ++i) loops[ip][i] = p0-i;
	}
      }else{
	if (p0>p1) swap(p0,p1);
	loops[ip].resize(p0-p1+2*n+1);
	for (int i=0; i<loops[ip].size(); ++i) loops[ip][i] = (p0-i+2*n)%(2*n);
      }
    }
  }
  bool AreLoopsCompatibleWithDiag(const Array<int,2>& diag){
    int n = diag.extent(0);
    for (int ip=0; ip<n; ++ip){
      int p0 = diag(ip,0);
      int p1 = diag(ip,1);
      int m=loops[ip].size();
      bool correct = (loops[ip][0]==p0 and loops[ip][m-1]==p1) or (loops[ip][0]==p1 and loops[ip][m-1]==p0);
      if (not correct){
	cout<<"ERROR loop["<<ip<<"]="<<loops[ip]<<" while d["<<ip<<"]=("<<diag(ip,0)<<","<<diag(ip,1)<<") is incompatible"<<endl;
	cout<<"diag="<<diag<<" loops="<<loops<<endl;
	return false;
      }
    }
    return true;
  }
  void PrintLoops(){
    cout<<"loops="<<endl;
    for (int ip=0; ip<loops.size(); ++ip){
      for (int j=0; j<loops[ip].size(); ++j)
	cout<<loops[ip][j]<<",";
      cout<<endl;
    }
  }
  void InitializeMomenta()
  {
    int n = loops.size();
    momenta.resize(2*n);
    momenta=0;
    for (int ip=0; ip<loops.size(); ++ip){
      int w = Freq(ip);
      for (int j=0;  j<loops[ip].size()-1; ++j){
	int dl = loops[ip][j]-loops[ip][j+1];
	if (dl==1 or dl==1-2*n){
	  momenta(loops[ip][j+1]) += w;
#ifdef _DEBUG
	  cout << "Adding G[" << loops[ip][j+1] << "]+=w[" << ip << "]=" << w << endl;
#endif	
	} else if (dl==-1 or dl==2*n-1){
	  momenta(loops[ip][j]) -= w;
#ifdef _DEBUG	
	  cout << "Adding G[" << loops[ip][j] << "]+=-w[" << ip << "]" << w << endl;
#endif	
	} else{
	  cout << "ERROR: Not yet implemented for loops!";
	}
      }
    }
  }
  void AddAtTheStart(deque<int>& loop, int cn, int co, int freq){
#ifdef _DEBUG  
    cout << "AddAtTheStart: cn=" << cn << ", co=" << co << endl;
#endif
    int twon = momenta.extent(0);
    if (cn>co){
      for (int j=co+1; j<cn+1; ++j) loop.push_front(j); //loop.insert(loop.begin(),j);
      for (int j=co; j<cn; ++j){
#ifdef _DEBUG
	cout << "k[" << j << "]+=" << freq << endl;
#endif
	momenta(j) += freq;
      }
    }else{
      for (int j=co+1; j<cn+1+twon; ++j) loop.push_front(j%twon); //loop.insert(loop.begin(),j%twon);
      for (int j=co; j<cn+twon; ++j){
	int jj = j % twon;
#ifdef _DEBUG	
	cout << "k[" << jj << "]+=" << freq << endl;
#endif      
	momenta(jj) += freq;
      }
    }
  }
  void AddAtTheEnd(deque<int>& loop, int cn, int co, int freq){
#ifdef _DEBUG  
    cout << "AddAtTheEnd: cn=" << cn << ", co=" <<  co << endl;
#endif  
    int twon = momenta.extent(0);
    if (cn<co){
      for (int j=co-1; j>cn-1; --j){
	loop.push_back(j); // loop.insert(loop.end(),j);
#ifdef _DEBUG
	cout << "k[" << j << "]+=" << freq << endl;
#endif      
	momenta(j) += freq;
      }
    } else{
      for (int j=co-1; j>cn-1-twon; --j){
	int jj = (j+twon)%(twon);
	loop.push_back(jj); // loop.insert(loop.end(), jj );
#ifdef _DEBUG
	cout << "k[" << jj << "]+=" << freq << endl;
#endif      
	momenta(jj) += freq;
      }
    }
  }
  void DeleteAtTheStart(deque<int>& loop, int cn, int co, int freq)
  {
#ifdef _DEBUG  
    cout << "DeletaAtTheStart: cn="<< cn << ", co="<< co << endl;
#endif  
    int twon = momenta.extent(0);
    if (cn<co){
      for (int j=co; j>cn; --j) loop.pop_front();
      for (int j=co-1; j>cn-1; --j){
	momenta(j) -= freq;
#ifdef _DEBUG      
	cout << "k[" << j << "]-=" << freq << endl;
#endif      
      }
    } else {
      for (int j=co; j>cn-twon; --j) loop.pop_front();
      for (int j=co-1; j>cn-1-twon; --j){
	int jj = (j+twon)%twon;
	momenta(jj) -= freq;
#ifdef _DEBUG      
	cout << "k[" << jj << "]-=" << freq << endl;
#endif
      }
    }
  }
  void DeleteAtTheEnd(deque<int>& loop, int cn, int co, int freq){
#ifdef _DEBUG  
    cout << "DeleteAtTheEnd co="<< co << ", cn=" <<  cn << endl;
#endif
    int twon = momenta.extent(0);
    if (cn>co){
      for (int j=co; j<cn; ++j){
	momenta(j) -= freq;
#ifdef _DEBUG      
	cout << "k[" << j << "]-=" << freq << endl;
#endif	
	loop.pop_back();
      }
    } else {
      for (int j=co-twon; j<cn; ++j){
	int jj = (j+twon) % twon;
#ifdef _DEBUG      
	cout << "k[" << jj << "]-=" << freq << endl;
#endif      
	momenta(jj) -= freq;
	loop.pop_back();
      }
    }
  }
  bool IsIn(int a, const deque<int>& loop) const
  {
    for (int i=0; i<loop.size(); ++i)
      if (loop[i]==a) return true;
    return false;
  }
  Array<int,2> Exchange(const TinyVector<int,4>& ih_jh, const Array<int,2>& diag, int& extra_sign) 
  {
    //#define _DEBUG    
    //cout<<"diag="<<diag<<" loops="<<loops<<endl;
    // here ih_jh is (1..n,0..1,1..n,0..1), hybridization places we want to exchange.
    Array<int,2> ndiag(diag.extent(0),diag.extent(1)); ndiag =  diag.copy();
    ndiag(ih_jh[0],ih_jh[1]) = diag(ih_jh[2],ih_jh[3]);
    ndiag(ih_jh[2],ih_jh[3]) = diag(ih_jh[0],ih_jh[1]);
    //if (debug){
    //  cout<<"diag="<<diag<<" ih_jh="<<ih_jh<<" ndiag="<<ndiag<<endl;
    //  PrintLoops();
    //}
    for (int i=0; i<2; ++i){
      int iloop = ih_jh[2*i];
      int ins   = ih_jh[2*i+1];
      int cn = ndiag(iloop,ins);
      int co =  diag(iloop,ins);
#ifdef _DEBUG    
      cout << "cn=" << cn << ", co=" << co << ", iloop=" << iloop << " ins=" << ins << " loops[iloop][0]="<<loops[iloop][0]<<" loops[iloop]="<<loops[iloop]<<endl;
#endif    
      if (loops[iloop][0]==co){ // we are changing the start
	if (not IsIn(cn, loops[iloop]))
	  AddAtTheStart(loops[iloop], cn, co, Freq(iloop));
	else
	  DeleteAtTheStart(loops[iloop], cn, co, Freq(iloop));
      }else{
	if (not IsIn(cn, loops[iloop]))
	  AddAtTheEnd(loops[iloop], cn, co, Freq(iloop));
	else{
	  //cout<<"Deleting at the end loops[iloop]="<<loops[iloop]<<" cn="<<cn<<" co="<<co<<endl;
	  DeleteAtTheEnd(loops[iloop], cn, co, Freq(iloop));
	}
      }
    }
    extra_sign=1;
    if (ndiag(ih_jh[0],0) > ndiag(ih_jh[0],1)) { swap( ndiag(ih_jh[0],0), ndiag(ih_jh[0],1) ); extra_sign *= -1;}
    if (ndiag(ih_jh[2],0) > ndiag(ih_jh[2],1)) { swap( ndiag(ih_jh[2],0), ndiag(ih_jh[2],1) ); extra_sign *= -1;}
    return ndiag;
  }
  bool IsReducible(){
    for (int i=0; i<momenta.extent(0); ++i)
      for (int j=i+1; j<momenta.extent(0); ++j)
	if (momenta(i)==momenta(j)){
	  //if (debug) cout<<"reducible"<<endl;
	  return true;
	}
    //if (debug) cout<<"irreducible"<<endl;
    return false;
  }
};

TinyVector<int,4> TwoRandomLoops(int n, const TinyVector<double,4>& rand){
  /*
  int li = static_cast<int>(drand()*n);
  int lj = static_cast<int>(drand()*(n-1));
  if (lj>=li) lj += 1;
  int ci = static_cast<int>(drand()*2);
  int cj = static_cast<int>(drand()*2);
  return TinyVector<int,4>(li,ci,lj,cj);
  */
  int li = static_cast<int>(rand[0]*n);
  int lj = static_cast<int>(rand[1]*(n-1));
  if (lj>=li) lj += 1;
  int ci = static_cast<int>(rand[2]*2);
  int cj = static_cast<int>(rand[3]*2);
  return TinyVector<int,4>(li,ci,lj,cj);
}
template<class RNG>
void FindNextSkeleton(int& sgn_dlt, Array<int,2>& diag, Skeletons& skl, RNG& drand)
{
  int Norder = diag.extent(0);
  if (Norder==1){
    sgn_dlt=1;
    return; // There is only one type of diagrams for NCA
  }
  Array<int,2> ndiag(Norder,2);
  bool Reducible = true;
  while (Reducible){
    TinyVector<int,4> ij;
    //ij = TwoRandomLoops(Norder);
    //ij = TwoRandomLoops(Norder, drand);
    ij = TwoRandomLoops(Norder, TinyVector<double,4>(drand(),drand(),drand(),drand()));
#ifdef _DEBUG    
    cout<<"  Exchanging=" << ij << endl;
#endif    
    int extra_sign=1;
    ndiag = skl.Exchange(ij, diag, extra_sign);
    Reducible = skl.IsReducible();
    sgn_dlt *= -extra_sign; // normally sign changes due to exchange of two deltas. But their direction (or order of two numbers in diag) can also change, hence extra sign.
    diag = ndiag.copy();
  }
}
#endif //_SKELETONS
