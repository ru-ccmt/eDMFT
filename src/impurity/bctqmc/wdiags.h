
inline string diag_p_2String(const Array<int,2>& ndiag, const Array<int,1>& npseudo)
{
  return fmt::format("{0};{1}", fmt::join(ndiag.data(), ndiag.data()+ndiag.size(), ","), fmt::join(npseudo.data(), npseudo.data()+npseudo.size(), ","));
}

vector<unsigned char> String2diagarray(string s0){
  deque<unsigned char> diagp;
  size_t ip=-1, iq;
  while ( ip+1 < s0.size()+1 ){
    string s1 = s0.substr(ip+1, 5);
    //cout<<"s1="<<s1<<endl;
    int ic = stoi(s1, &iq);
    ip+=iq+1;
    diagp.push_back( ic );
    //cout<<"ip="<<ip<<" l="<<s0.size()<<endl;
  }
  int Norder = diagp.size()/4;
  vector<unsigned char> diag_pseud(4*Norder);
  for (int i=0; i<4*Norder; ++i) diag_pseud[i] = diagp[i];
  return diag_pseud;
}


auto String2diagp(string s0)
{
  struct retVals {
    Array<int,2> diag;
    Array<int,1> pseud;
  };
  //cout<<s0<<endl;
  deque<int> diagp;
  size_t ip=-1, iq;
  while ( ip+1 < s0.size()+1 ){
    string s1 = s0.substr(ip+1, 5);
    //cout<<"s1="<<s1<<endl;
    int ic = stoi(s1, &iq);
    ip+=iq+1;
    diagp.push_back( ic );
    //cout<<"ip="<<ip<<" l="<<s0.size()<<endl;
  }
  int Norder = diagp.size()/4;
  //cout <<"Norder="<<Norder << endl;
  retVals r;
  r.diag.resize(Norder,2);
  r.pseud.resize(2*Norder);
  for (int i=0; i<Norder; ++i){
    r.diag(i,0) = diagp[2*i];
    r.diag(i,1) = diagp[2*i+1];
  }
  for (int i=0; i<2*Norder; ++i)
    r.pseud(i) = diagp[i+2*Norder];
  return r;
}

Array<int,2> sorted_diag(const Array<int,2>& diag)
{
  int n = diag.extent(0);
  // initialize original index locations
  vector<int> idx(n);
  iota(idx.begin(), idx.end(), 0);
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
	      [&diag](int i1, int i2) {return diag(i1,0) < diag(i2,0);});
  
  Array<int,2> ndiag(n,2); ndiag = diag.copy();
  for (int i=0; i<n; ++i) ndiag(i,Range::all()) = diag(idx[i],Range::all());
  return ndiag;
}
void sorted_diag2(Array<int,2>& diag)
{
  int n = diag.extent(0);
  // initialize original index locations
  vector<int> idx(n);
  iota(idx.begin(), idx.end(), 0);
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
	      [&diag](int i1, int i2) {return diag(i1,0) < diag(i2,0);});
  
  Array<int,2> ndiag(n,2); ndiag = diag;
  for (int i=0; i<n; ++i) {ndiag(i,0)=diag(idx[i],0); ndiag(i,1)=diag(idx[i],1);}
  for (int i=0; i<n; ++i) {diag(i,0)=ndiag(i,0); diag(i,1)=ndiag(i,1);}
}

Array<int,2> unflip_diag(const Array<int,2>& diag)
{
  int n = diag.extent(0);
  Array<int,2> ndiag(n,2);
  ndiag = diag.copy();
  for (int i=0; i<n; ++i)
    if (ndiag(i,0)>ndiag(i,1)) swap(ndiag(i,0),ndiag(i,1));
  return ndiag;
}
auto Different3(const Array<int,2>& diag, const Array<int,1>& cpseudo, double P_acc, map<string,double>& stored){
  struct retVals {        
    //double diag_weight;
    bool new_diag;
    string s;
  };
  int n = diag.extent(0);
  Array<int,2> ndiag(n,2);
  for (int j=0; j<n; ++j){
    if (diag(j,0)<diag(j,1)){
      ndiag(j,0)=diag(j,0);
      ndiag(j,1)=diag(j,1);
    }else{
      ndiag(j,0)=diag(j,1);
      ndiag(j,1)=diag(j,0);
    }
  }
  sorted_diag2(ndiag);
  string s0 = diag_p_2String(ndiag,cpseudo);
  if (stored.size()==0){
    stored[s0] = P_acc;
    return retVals {true,s0}; // new_diag
  }
  auto itt = stored.find(s0);
  if (itt!=stored.end()){
    itt->second += P_acc;
    return retVals {false,s0}; // old_diag
  }
  Array<int,1> npseudo(2*n);
  Array<int,2> sdiag(n,2);
  for (int ir=1; ir<2*n; ++ir){
    for(int j=0; j<n; ++j){
      sdiag(j,0) = (ndiag(j,0)+ir)%(2*n);
      sdiag(j,1) = (ndiag(j,1)+ir)%(2*n);
      if (sdiag(j,0)>sdiag(j,1)) swap(sdiag(j,0),sdiag(j,1));
    }
    sorted_diag2(sdiag);
    for(int j=0; j<2*n; ++j) npseudo(j) = cpseudo((j-ir+2*n)%(2*n));
    string s = diag_p_2String(sdiag, npseudo);
    auto itt = stored.find(s);
    if (itt!=stored.end()){
      itt->second += P_acc;
      return retVals {false,s}; // old_diag
    }
  }
  stored[s0] = P_acc;
  return retVals {true,s0}; // new_diag
}
int SymmetryFactor2(const pair<Array<int,2>,Array<int,1>>& diag_pseudo){
  Array<int,2> diag = diag_pseudo.first;
  Array<int,1> cpseudo = diag_pseudo.second;
  
  int n = diag.extent(0);
  Array<int,2> ndiag(n,2);
  Array<int,1> npseudo(2*n);
  int sym=1;
  for (int r=1; r<2*n; ++r){
    for (int j=0; j<n; ++j){
      ndiag(j,Range::all()) = (diag(j,Range::all())+r)%(2*n);
      if (ndiag(j,0)>ndiag(j,1)) swap(ndiag(j,0),ndiag(j,1));
    }
    for (int j=0; j<2*n; ++j){
      npseudo(j) = cpseudo((j-r+2*n)%(2*n));
    }
    //Array<int,2> sdiag = sorted_diag(ndiag);
    //if (sdiag == diag and npseudo == cpseudo) sym+=1;
    sorted_diag2(ndiag);
    if (ndiag == diag and npseudo == cpseudo) sym+=1;
  }
  return sym;
}
TinyVector<int,4> GiveVert(const Array<int,2>& diag, int p, int two_n){
  return TinyVector<int,4>( (diag(p,0)-1+two_n)%two_n, diag(p,0), (diag(p,1)-1+two_n)%two_n, diag(p,1) );
}
void FindPossibleDeltas(deque<TinyVector<int,4>>& possible, const TinyVector<int,4>& mms, const vector<TinyVector<int,4>>& non_zero)
{
  for (int ii=0; ii<non_zero.size(); ++ii){
    TinyVector<int,4> ms = non_zero[ii];
    // over all 4 vertices of hybridization line
    // either this vertex has not been set or it is set to conincide with this Delta entry
    bool avail = (mms[0]<0 or (mms[0]>=0 and ms[0]==mms[0]) ) and (mms[1]<0 or (mms[1]>=0 and ms[1]==mms[1]) ) and (mms[2]<0 or (mms[2]>=0 and ms[2]==mms[2]) ) and (mms[3]<0 or (mms[3]>=0 and ms[3]==mms[3]) );
    if (avail) possible.push_back(ms);
  }
}
void recursive_n(int p, Array<int,1>& pseudo1, int n, const Array<int,2>& diag, const vector<TinyVector<int,4>>& non_zero, const Array<int,4>& Dlt, deque<Array<int,1>>& pseudos){
  //bool debug=true;
  TinyVector<int,4> vert1 = GiveVert(diag, p, 2*n);
  TinyVector<int,4> mms(pseudo1(vert1[0]),pseudo1(vert1[1]),pseudo1(vert1[2]),pseudo1(vert1[3]));
  bool all_positive = mms[0]>0 and mms[1]>0 and mms[2]>0 and mms[3]>0;
  deque<TinyVector<int,4>> possible2;
  if (all_positive){
    possible2.push_back( mms );
  }else{
    FindPossibleDeltas(possible2, mms, non_zero);
  }
  for (int ii2=0; ii2<possible2.size(); ++ii2){
    TinyVector<int,4> ms = possible2[ii2];
    int dd = Dlt(ms);
    if (dd == 0) break;
    Array<int,1> pseudo2(2*n);
    pseudo2 = pseudo1.copy();
    // for (int l=0; l<4; ++l){
    //   int pl = pseudo1(vert1[l]);
    //   if (pl>0 and pl!=ms[l]){
    // 	cout<<"ERROR pl="<<pl<<" while ms[l]="<<ms[l]<<" pseudo1="<<pseudo1<<" ms="<<ms<<" vert1="<<vert1<<endl;
    // 	exit(1);
    //   }
    // }
    pseudo2(vert1[0]) = ms[0];
    pseudo2(vert1[1]) = ms[1];
    pseudo2(vert1[2]) = ms[2];
    pseudo2(vert1[3]) = ms[3];
    if (p==n-1){
      pseudos.push_back( pseudo2 );
      //if (debug) cout << "Final-> " << ii2 << " " << ms << " " << pseudo2 << endl;
      return;
    }else{
      //if (debug) cout << "pseudo"<<(p+1)<<"="<<pseudo2<<" pseudo"<<p<<"="<< pseudo1 << endl;
      recursive_n(p+1, pseudo2, /*Delta_type,*/ n, diag, non_zero, Dlt, /*Delta_types,*/ pseudos);
    }
  }
}

void Determine_all_diagrams_of_this_topology(deque<Array<int,1>>& pseudos, Array<int,2>& Delta_types, const Array<int,2>& diag, const vector<TinyVector<int,4>>& non_zero, const Array<int,4>& Dlt)
{
  pseudos.clear();
  int n = diag.extent(0);
  Array<int,1> pseudo1(2*n); pseudo1 = -1;
  int p=0;
  TinyVector<int,4> vert0 = GiveVert(diag, p, 2*n);
  for (int ii1=0; ii1<non_zero.size(); ++ii1){
    TinyVector<int,4> ms = non_zero[ii1];
    for (int i=0; i<4; ++i) pseudo1(vert0[i]) = ms[i];
    //cout<<"vert0="<<vert0<<" ms=" << ms<<" pseudo1="<<pseudo1<<endl;
    if (n==1){
      if (vert0(0)==vert0(3) and ms[0]!=ms[3]) continue;
      if (vert0(1)==vert0(2) and ms[1]!=ms[2]) continue;
      int dd = Dlt(ms);
      if (dd != 0){
	pseudos.push_back( pseudo1.copy() );
	//cout << "Final-> " << ii1 << " " << ms << " " << pseudo1 << endl;
      }
    }else{
      recursive_n(1, pseudo1, n, diag, non_zero, Dlt, pseudos);
    }
  }
  Delta_types.resize(pseudos.size(),n);
  for (int ip=0; ip<pseudos.size(); ++ip){
    //Delta_types[ip].resize(n);
    Array<int,1>& pseudo1 = pseudos[ip];
    for (int p=0; p<n; ++p){
      TinyVector<int,4> vert = GiveVert(diag, p, 2*n);
      TinyVector<int,4> ms(pseudo1(vert[0]),pseudo1(vert[1]),pseudo1(vert[2]),pseudo1(vert[3]));
      Delta_types(ip,p) = Dlt(ms);
    }
    //cout<<"Final2->"<<ip<<" p="<<pseudo1<<" Dlt="<<Delta_types[ip]/*<<" from Dt="<<Delta_types[ip]*/<<endl;
  }
}


void recursive_n2(int p, Array<int,1>& pseudo1, int n, const Array<int,2>& diag, const vector<TinyVector<int,4>>& non_zero, const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip, deque<Array<int,1>>& pseudos){
  //bool debug=true;
  int two_n = 2*n;
  TinyVector<int,4> vert1 = GiveVert(diag, p, two_n);
  TinyVector<int,4> mms(pseudo1(vert1[0]),pseudo1(vert1[1]),pseudo1(vert1[2]),pseudo1(vert1[3]));
  // check that successive atomic state differ by hybridization line either comming in or going out.
  // If this is not the case, we can immediately terminate such resursion.
  for (int j=0; j<two_n-1; ++j){
    TinyVector<int,2> nns(pseudo1(j),pseudo1(j+1));
    if (nns[0]>=0 and nns[1]>=0)
      if (Fi_Fip(nns)==0)return;
  }
  {
    TinyVector<int,2> nns(pseudo1(two_n-1),pseudo1(0));
    if (nns[0]>=0 and nns[1]>=0)
      if (Fi_Fip(nns)==0) return;
  }
  bool all_positive = mms[0]>0 and mms[1]>0 and mms[2]>0 and mms[3]>0;
  deque<TinyVector<int,4>> possible2;
  if (all_positive)
    possible2.push_back( mms );
  else
    FindPossibleDeltas(possible2, mms, non_zero);

  //cout<<"p="<<p<<" vert1="<<vert1<<" mms="<<mms<<" pseudo="<<pseudo1<<" possible2={";
  //for (int i=0; i<possible2.size(); ++i) cout<< possible2[i]<<";";
  //cout<<"}"<<endl;
  
  for (int ii2=0; ii2<possible2.size(); ++ii2){
    TinyVector<int,4> ms = possible2[ii2];
    int dd = Dlt(ms);
    if (dd == 0) break;
    Array<int,1> pseudo2(2*n);
    pseudo2 = pseudo1.copy();
    pseudo2(vert1[0]) = ms[0];
    pseudo2(vert1[1]) = ms[1];
    pseudo2(vert1[2]) = ms[2];
    pseudo2(vert1[3]) = ms[3];
    if (p==n-1){
      pseudos.push_back( pseudo2 );
      //if (debug) cout << "Final-> " << ii2 << " " << ms << " " << pseudo2 << endl;
      return;
    }else{
      //if (debug) cout << "pseudo"<<(p+1)<<"="<<pseudo2<<" pseudo"<<p<<"="<< pseudo1 << endl;
      recursive_n2(p+1, pseudo2, /*Delta_type,*/ n, diag, non_zero, Dlt, Fi_Fip, pseudos);
    }
  }
}
void Determine_all_diagrams_of_this_topology2(deque<Array<int,1>>& pseudos, Array<int,2>& Delta_types, const Array<int,2>& diag, const vector<TinyVector<int,4>>& non_zero, const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip)
{
  pseudos.clear();
  int n = diag.extent(0);
  Array<int,1> pseudo1(2*n); pseudo1 = -1;
  int p=0;
  TinyVector<int,4> vert0 = GiveVert(diag, p, 2*n);
  for (int ii1=0; ii1<non_zero.size(); ++ii1){
    TinyVector<int,4> ms = non_zero[ii1];
    for (int i=0; i<4; ++i) pseudo1(vert0[i]) = ms[i];
    //cout<<"vert0="<<vert0<<" ms=" << ms<<" pseudo1="<<pseudo1<<endl;
    if (n==1){
      if (vert0(0)==vert0(3) and ms[0]!=ms[3]) continue;
      if (vert0(1)==vert0(2) and ms[1]!=ms[2]) continue;
      int dd = Dlt(ms);
      if (dd != 0){
	pseudos.push_back( pseudo1.copy() );
	//cout << "Final-> " << ii1 << " " << ms << " " << pseudo1 << endl;
      }
    }else{
      recursive_n2(1, pseudo1, n, diag, non_zero, Dlt, Fi_Fip, pseudos);
    }
  }
  Delta_types.resize(pseudos.size(),n);
  for (int ip=0; ip<pseudos.size(); ++ip){
    //Delta_types[ip].resize(n);
    Array<int,1>& pseudo1 = pseudos[ip];
    for (int p=0; p<n; ++p){
      TinyVector<int,4> vert = GiveVert(diag, p, 2*n);
      TinyVector<int,4> ms(pseudo1(vert[0]),pseudo1(vert[1]),pseudo1(vert[2]),pseudo1(vert[3]));
      Delta_types(ip,p) = Dlt(ms);
    }
    //cout<<"Final2->"<<ip<<" p="<<pseudo1<<" Dlt="<<Delta_types[ip]/*<<" from Dt="<<Delta_types[ip]*/<<endl;
  }
}
void recursive_n3(int p, Array<int,1>& pseudo1, int n, const Array<int,1>& order, const Array<int,2>& diag, const vector<TinyVector<int,4>>& non_zero, const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip, deque<Array<int,1>>& pseudos){
  //#define _DEBUG  
  bool debug=false;
  int two_n = 2*n;
  TinyVector<int,4> vert1 = GiveVert(diag, order(p), two_n);
  TinyVector<int,4> mms(pseudo1(vert1[0]),pseudo1(vert1[1]),pseudo1(vert1[2]),pseudo1(vert1[3]));
  // check that successive atomic state differ by hybridization line either comming in or going out.
  // If this is not the case, we can immediately terminate such resursion.
  for (int j=0; j<two_n-1; ++j){
    TinyVector<int,2> nns(pseudo1(j),pseudo1(j+1));
    if (nns[0]>=0 and nns[1]>=0)
      if (Fi_Fip(nns)==0)return;
  }
  {
    TinyVector<int,2> nns(pseudo1(two_n-1),pseudo1(0));
    if (nns[0]>=0 and nns[1]>=0)
      if (Fi_Fip(nns)==0) return;
  }
  bool all_positive = mms[0]>0 and mms[1]>0 and mms[2]>0 and mms[3]>0;
  deque<TinyVector<int,4>> possible2;
  if (all_positive)
    possible2.push_back( mms );
  else
    FindPossibleDeltas(possible2, mms, non_zero);
#ifdef _DEBUG
  cout<<"p="<<p<<" idi="<<order(p)<<" vert1="<<vert1<<" mms="<<mms<<" pseudo="<<pseudo1<<" possible2={";
  for (int i=0; i<possible2.size(); ++i) cout<< possible2[i]<<";";
  cout<<"}"<<endl;
#endif
  for (int ii2=0; ii2<possible2.size(); ++ii2){
    TinyVector<int,4> ms = possible2[ii2];
    int dd = Dlt(ms);
    if (dd == 0) break;
    Array<int,1> pseudo2(2*n);
    pseudo2 = pseudo1.copy();
    pseudo2(vert1[0]) = ms[0];
    pseudo2(vert1[1]) = ms[1];
    pseudo2(vert1[2]) = ms[2];
    pseudo2(vert1[3]) = ms[3];
    if (p==n-1){
      pseudos.push_back( pseudo2 );
#ifdef _DEBUG      
      cout << "Final-> " << ii2 << " " << ms << " " << pseudo2 << endl;
#endif      
      return;
    }else{
#ifdef _DEBUG      
      cout << "    pseudo"<<(p+1)<<"="<<pseudo2<<" pseudo"<<p<<"="<< pseudo1 << endl;
#endif      
      recursive_n3(p+1, pseudo2, n, order, diag, non_zero, Dlt, Fi_Fip, pseudos);
    }
  }
}
void Determine_diagrams_of_this_topology_with(const TinyVector<int,2>& p_known, deque<Array<int,1>>& pseudos, Array<int,2>& Delta_types,
					      const Array<int,2>& diag, const vector<TinyVector<int,4>>& non_zero, const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip)
{
  //#define _DEBUG  
  pseudos.clear();
  int n = diag.extent(0);
  int two_n=2*n;
  Array<int,1> pseudo1(2*n); pseudo1 = -1;

  pseudo1(p_known[0])=p_known[1]; // This is the one we know for sure

  Array<int,1> order(n);
  for (int i=0; i<n; ++i) order(i) = i;// the order in which we are going to visit hybridizations
  if (n>1){
    for (int idi=0; idi<n; ++idi)
      if (diag(idi,0)==p_known[0] or diag(idi,1)==p_known[0]){
	swap(order(0),order(idi));
	break;
      }
    for (int idi=0; idi<n; ++idi)
      if (diag(idi,0)==(p_known[0]+1)%two_n or diag(idi,1)==(p_known[0]+1)%two_n){
	swap(order(1),order(idi));
	break;
      }
  }
#ifdef _DEBUG  
  cout<<"order="<<order<<" pseudo1="<<pseudo1<<endl;
#endif    
  int p=0;
  TinyVector<int,4> vert0 = GiveVert(diag, order(p), two_n);
  TinyVector<int,4> mms(pseudo1(vert0[0]),pseudo1(vert0[1]),pseudo1(vert0[2]),pseudo1(vert0[3])); // here we know that mms[1] was set in the input, the rest are -1
  deque<TinyVector<int,4>> possible;
  FindPossibleDeltas(possible, mms, non_zero);
#ifdef _DEBUG
  cout<<"vert0="<<vert0<<" mms="<<mms<<" possible="<<endl;
  for (int i=0; i<possible.size(); ++i)
    cout<<"   possible["<<i<<"]="<<possible[i]<<endl;
#endif  
  for (int ii1=0; ii1<possible.size(); ++ii1){
    TinyVector<int,4> ms = possible[ii1];
    int dd = Dlt(ms);
    if (dd == 0) break;
    Array<int,1> pseudo2(two_n);
    pseudo2 = pseudo1.copy();
    pseudo2(vert0[0]) = ms[0];
    pseudo2(vert0[1]) = ms[1];
    pseudo2(vert0[2]) = ms[2];
    pseudo2(vert0[3]) = ms[3];
    if (n==1){
      if (vert0(0)==vert0(3) and ms[0]!=ms[3]) continue;
      if (vert0(1)==vert0(2) and ms[1]!=ms[2]) continue;
      if (dd != 0){
	pseudos.push_back( pseudo2.copy() );
	//cout << "Final-> " << ii1 << " " << ms << " " << pseudo2 << endl;
      }
    }else{
      //cout<<"ii1="<<ii1<<" ms="<<ms<<" dd="<<dd<<" pseudo2="<<pseudo2<<endl;
      recursive_n3(1, pseudo2, n, order, diag, non_zero, Dlt, Fi_Fip, pseudos);
    }
  }
  Delta_types.resize(pseudos.size(),n);
  for (int ip=0; ip<pseudos.size(); ++ip){
    //Delta_types[ip].resize(n);
    Array<int,1>& pseudo1 = pseudos[ip];
    for (int p=0; p<n; ++p){
      TinyVector<int,4> vert = GiveVert(diag, p, 2*n);
      TinyVector<int,4> ms(pseudo1(vert[0]),pseudo1(vert[1]),pseudo1(vert[2]),pseudo1(vert[3]));
      Delta_types(ip,p) = Dlt(ms);
    }
#ifdef _DEBUG    
    cout<<"Final2->"<<ip<<" p="<<pseudo1<<" Dlt="<<Delta_types(ip,Range::all())/*<<" from Dt="<<Delta_types[ip]*/<<endl;
#endif    
  }
}

auto GiveCommonRandomPseudoPropagator(const TinyVector<int,2>& Morder, const Array<int,1>& cpseudo, const Array<int,1>& o2n, double xrand)
{
  struct retVals{
    TinyVector<int,2> p_known;
    TinyVector<int,2> q_known;
  };
  int i_known=-1;
  if (Morder[0]!=Morder[1]){
    deque<int> possible_j;
    for (int j=0; j<2*Morder[0]; ++j){
      int jt = o2n(j);
      int kt = o2n((j+1)%(2*Morder[0])); // two pseudo-particle propagators are changed when we remove a time at certain vertex, hence also kn should be positive
      bool nothing_inbetween = ((kt-jt==1) or (kt-jt==-2*Morder[1]+1)) and (jt>=0 and kt>=0);
      if (nothing_inbetween) possible_j.push_back(j);
    }
    i_known = possible_j[static_cast<int>(xrand*possible_j.size())];
  }else{
    i_known = xrand*2*Morder[0];
  }
  retVals r;
  r.p_known = i_known,cpseudo(i_known);
  if (Morder[0]!=Morder[1])
    r.q_known = o2n(i_known), cpseudo(i_known); //@@
  else
    r.q_known = i_known,cpseudo(i_known);
  return r;
}
auto GiveCommonRandomPseudoPropagator_of_size(int N, const TinyVector<int,2>& Morder, const Array<int,1>& cpseudo, const Array<int,1>& o2n, double xrand)
{//  We are searching for pseudo-propagators which are common in both diagrams (when we have two orders simultaneously present, i.e., Morder[0]!=Morder[1]).
 //  We are requesting at most N common propagators, which are returned in the following way:
 //    d0_known contains indices of propagators in the first diagram
 //    d1_known contains indices of propagaotrs in the second diagram
 //    p_known contains pseudoparticle labels that are common in both diagrams and share the same timings
 //  There should always be at least one common propagator when diagrams are different (Morder[0]!=Morder[1])
  struct retVals{
    vector<int> p_known;  // pseudoparticle indices
    vector<int> d0_known; // indices in the old diagram
    vector<int> d1_known; // indices in the new diagram
  };
  retVals r;
  if (Morder[0]!=Morder[1]){
    deque<int> possible_j;
    for (int j=0; j<2*Morder[0]; ++j){
      int jt = o2n(j);
      int kt = o2n((j+1)%(2*Morder[0])); // two pseudo-particle propagators are changed when we remove a time at certain vertex, hence also kn should be positive
      bool nothing_inbetween = ((kt-jt==1) or (kt-jt==-2*Morder[1]+1)) and (jt>=0 and kt>=0);
      if (nothing_inbetween) possible_j.push_back(j);
    }
    if (possible_j.size()==0) return r;
    int i_known = static_cast<int>(xrand*possible_j.size());
    int m=possible_j.size();
    int n=min(N,m);
    r.d0_known.resize(n);
    r.d1_known.resize(n);
    r.p_known.resize(n);
    for (int i=0; i<n; ++i){
      int k = possible_j[(i_known+i)%m];
      r.d0_known[i] = k;
      r.d1_known[i] = o2n(k);
      r.p_known[i] = cpseudo(k);
    }
  }else{
    int i_known = xrand*2*Morder[0];
    int m=2*Morder[0];
    int n=min(N,m);
    r.d0_known.resize(n);
    r.d1_known.resize(n);
    r.p_known.resize(n);
    for (int i=0; i<n; ++i){
      int k = (i_known+i)%m;
      r.d0_known[i] = k;
      r.d1_known[i] = k;
      r.p_known[i] = cpseudo(k);
    }
  }
  return r;
}



void Determine_diagrams_of_this_topology_with(const deque<int> p_known, int i_known_start,
					      deque<Array<int,1>>& pseudos, Array<int,2>& Delta_types,
					      const Array<int,2>& diag, const vector<TinyVector<int,4>>& non_zero, const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip)
{
  //#define _DEBUG  
  pseudos.clear();
  int n = diag.extent(0);
  int two_n=2*n;
  Array<int,1> pseudo1(2*n); pseudo1 = -1;

  for (int i=0; i<p_known.size(); ++i)
    pseudo1((i_known_start+i)%two_n)=p_known[i]; // This is the one we know for sure

  Array<int,1> order(n);
  iota(order.data(), order.data()+n, 0);
  if (n>1){
    vector<int> vert_count(n);
    for (int i=0; i<n; ++i) vert_count[i]=0;
    for (int idi=0; idi<n; ++idi){
      TinyVector<int,4> vert = GiveVert(diag, idi, two_n);
      for (int j=0; j<p_known.size(); ++j){
	int i_known_position = (i_known_start+j)%two_n;
	for (int k=0; k<4; ++k)
	  if (vert[k]==i_known_position)
	    vert_count[idi]+=1;
      }
    }
    stable_sort(order.data(), order.data()+n,
		[&vert_count](size_t i1, size_t i2) {return vert_count[i1] > vert_count[i2];});
#ifdef _DEBUG  
    cout<<"vert_count="<<vert_count<<" order="<<order<<" pseudo1="<<pseudo1<<endl;
#endif    
  }
  int p=0;
  TinyVector<int,4> vert0 = GiveVert(diag, order(p), two_n);
  TinyVector<int,4> mms(pseudo1(vert0[0]),pseudo1(vert0[1]),pseudo1(vert0[2]),pseudo1(vert0[3])); // here we know that mms[1] was set in the input, the rest are -1
  deque<TinyVector<int,4>> possible;
  FindPossibleDeltas(possible, mms, non_zero);
#ifdef _DEBUG
  cout<<"vert0="<<vert0<<" mms="<<mms<<" possible="<<endl;
  for (int i=0; i<possible.size(); ++i)
    cout<<"   possible["<<i<<"]="<<possible[i]<<endl;
#endif  
  for (int ii1=0; ii1<possible.size(); ++ii1){
    TinyVector<int,4> ms = possible[ii1];
    int dd = Dlt(ms);
    if (dd == 0) break;
    Array<int,1> pseudo2(two_n);
    pseudo2 = pseudo1.copy();
    pseudo2(vert0[0]) = ms[0];
    pseudo2(vert0[1]) = ms[1];
    pseudo2(vert0[2]) = ms[2];
    pseudo2(vert0[3]) = ms[3];
    if (n==1){
      if (vert0(0)==vert0(3) and ms[0]!=ms[3]) continue;
      if (vert0(1)==vert0(2) and ms[1]!=ms[2]) continue;
      if (dd != 0){
	pseudos.push_back( pseudo2.copy() );
	//cout << "Final-> " << ii1 << " " << ms << " " << pseudo2 << endl;
      }
    }else{
      //cout<<"ii1="<<ii1<<" ms="<<ms<<" dd="<<dd<<" pseudo2="<<pseudo2<<endl;
      recursive_n3(1, pseudo2, n, order, diag, non_zero, Dlt, Fi_Fip, pseudos);
    }
  }
  Delta_types.resize(pseudos.size(),n);
  for (int ip=0; ip<pseudos.size(); ++ip){
    //Delta_types[ip].resize(n);
    Array<int,1>& pseudo1 = pseudos[ip];
    for (int p=0; p<n; ++p){
      TinyVector<int,4> vert = GiveVert(diag, p, 2*n);
      TinyVector<int,4> ms(pseudo1(vert[0]),pseudo1(vert[1]),pseudo1(vert[2]),pseudo1(vert[3]));
      Delta_types(ip,p) = Dlt(ms);
    }
#ifdef _DEBUG    
    cout<<"Final2->"<<ip<<" p="<<pseudo1<<" Dlt="<<Delta_types(ip,Range::all())/*<<" from Dt="<<Delta_types[ip]*/<<endl;
#endif    
  }
}

inline bool IsIn(int a, const vector<int>& d_known)
{
  return std::find(d_known.begin(), d_known.end(), a) != d_known.end();
}
void Determine_diagrams_of_this_topology_with_option1(const vector<int>& d_known, const vector<int>& p_known,
						      deque<Array<int,1>>& pseudos, Array<int,2>& Delta_types,
						      const Array<int,2>& diag, const vector<TinyVector<int,4>>& non_zero, const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip)
{
  //#define _DEBUG  
  pseudos.clear();
  int n = diag.extent(0);
  int two_n=2*n;
  Array<int,1> pseudo1(2*n); pseudo1 = -1;

  for (int i=0; i<d_known.size(); ++i) pseudo1(d_known[i])=p_known[i]; // This is the one we know for sure
  // In which order should we process hybridization propagators?
  // It is best to do it with those that are already fixed, because recursion will have much less to do.
  Array<int,1> order(n);
  iota(order.data(), order.data()+n, 0);
  if (n>1){
    vector<int> vert_count(n);
    // this will contain how many pseudo-propagators are fixed for each hybridization.
    // Each hybridization has four pseudo-propagators attached, and their could be between 0 and 4 of them fixed by income data.
    for (int i=0; i<n; ++i) vert_count[i]=0; 
    for (int idi=0; idi<n; ++idi){  // over all hybridization lines (index idi)
      TinyVector<int,4> vert = GiveVert(diag, idi, two_n); // which vertices this hybridization line uses?
      for (int k=0; k<4; ++k)  // over all vertices that are involved in hybridization line with index idi
	if (IsIn(vert[k],d_known)) vert_count[idi]+=1; // is this vertex present in d_known indices
    }
    stable_sort(order.data(), order.data()+n,
		[&vert_count](size_t i1, size_t i2) {return vert_count[i1] > vert_count[i2];});
#ifdef _DEBUG  
    cout<<"vert_count="<<vert_count<<" order="<<order<<" pseudo1="<<pseudo1<<endl;
#endif    
  }
  int p=0;
  TinyVector<int,4> vert0 = GiveVert(diag, order(p), two_n);
  TinyVector<int,4> mms(pseudo1(vert0[0]),pseudo1(vert0[1]),pseudo1(vert0[2]),pseudo1(vert0[3])); // here we know that mms[1] was set in the input, the rest are -1
  deque<TinyVector<int,4>> possible;
  FindPossibleDeltas(possible, mms, non_zero);
#ifdef _DEBUG
  cout<<"p="<<p<<" idi="<<order(p)<<" vert0="<<vert0<<" mms="<<mms<<" pseudo="<<pseudo1<<" possible={";
  for (int i=0; i<possible.size(); ++i) cout<< possible[i]<<";";
  cout<<"}"<<endl;
#endif  
  for (int ii1=0; ii1<possible.size(); ++ii1){
    TinyVector<int,4> ms = possible[ii1];
    int dd = Dlt(ms);
    if (dd == 0) break;
    Array<int,1> pseudo2(two_n);
    pseudo2 = pseudo1.copy();
    for (int k=0; k<4; ++k) pseudo2(vert0[k]) = ms[k];
    //#ifdef _DEBUG
    //cout<<p<<") vert="<<vert0<<" ms="<<ms<<" vert_count2="<<vert_count2<<endl;
    //#endif    
    if (n==1){
      if (vert0(0)==vert0(3) and ms[0]!=ms[3]) continue;
      if (vert0(1)==vert0(2) and ms[1]!=ms[2]) continue;
      if (dd != 0){
	pseudos.push_back( pseudo2.copy() );
	//cout << "Final-> " << ii1 << " " << ms << " " << pseudo2 << endl;
      }
    }else{
      //cout<<"ii1="<<ii1<<" ms="<<ms<<" dd="<<dd<<" pseudo2="<<pseudo2<<endl;
      recursive_n3(1, pseudo2, n, order, diag, non_zero, Dlt, Fi_Fip, pseudos);
    }
  }
  Delta_types.resize(pseudos.size(),n);
  for (int ip=0; ip<pseudos.size(); ++ip){
    //Delta_types[ip].resize(n);
    Array<int,1>& pseudo1 = pseudos[ip];
    for (int p=0; p<n; ++p){
      TinyVector<int,4> vert = GiveVert(diag, p, 2*n);
      TinyVector<int,4> ms(pseudo1(vert[0]),pseudo1(vert[1]),pseudo1(vert[2]),pseudo1(vert[3]));
      Delta_types(ip,p) = Dlt(ms);
    }
#ifdef _DEBUG    
    cout<<"Final2->"<<ip<<" p="<<pseudo1<<" Dlt="<<Delta_types(ip,Range::all())/*<<" from Dt="<<Delta_types[ip]*/<<endl;
#endif    
  }
}

void recursive_n4(int p, Array<int,1>& pseudo1, int n, Array<int,1>& vert_count1, const Array<int,2>& connected_hybrid, const Array<int,2>& diag, const vector<TinyVector<int,4>>& non_zero, const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip, deque<Array<int,1>>& pseudos){
  //#define _DEBUG  
  bool debug=false;
  int two_n = 2*n;

  int idi_next = max_element(vert_count1.data(), vert_count1.data()+n)-vert_count1.data();
  vert_count1(idi_next)=-1; // We eliminate this hybridization by setting count to negative number.
  
  TinyVector<int,4> vert1 = GiveVert(diag, idi_next, two_n);
  TinyVector<int,4> mms(pseudo1(vert1[0]),pseudo1(vert1[1]),pseudo1(vert1[2]),pseudo1(vert1[3]));
  // check that successive atomic state differ by hybridization line either comming in or going out.
  // If this is not the case, we can immediately terminate such resursion.
  for (int j=0; j<two_n-1; ++j){
    TinyVector<int,2> nns(pseudo1(j),pseudo1(j+1));
    if (nns[0]>=0 and nns[1]>=0)
      if (Fi_Fip(nns)==0)return;
  }
  {
    TinyVector<int,2> nns(pseudo1(two_n-1),pseudo1(0));
    if (nns[0]>=0 and nns[1]>=0)
      if (Fi_Fip(nns)==0) return;
  }
  bool all_positive = mms[0]>0 and mms[1]>0 and mms[2]>0 and mms[3]>0;
  deque<TinyVector<int,4>> possible2;
  if (all_positive)
    possible2.push_back( mms );
  else
    FindPossibleDeltas(possible2, mms, non_zero);
#ifdef _DEBUG
  cout<<"p="<<p<<" idi="<<idi_next<<" vert_count1="<<vert_count1<<" vert1="<<vert1<<" mms="<<mms<<" pseudo="<<pseudo1<<" possible2={";
  for (int i=0; i<possible2.size(); ++i) cout<< possible2[i]<<";";
  cout<<"}"<<endl;
#endif  
  for (int ii2=0; ii2<possible2.size(); ++ii2){
    TinyVector<int,4> ms = possible2[ii2];
    int dd = Dlt(ms);
    if (dd == 0) break;
    Array<int,1> pseudo2(2*n);
    pseudo2 = pseudo1.copy();
    for (int k=0; k<4; ++k) pseudo2(vert1[k]) = ms[k];
    //#ifdef _DEBUG
    //cout<<p<<") vert="<<vert1<<" ms="<<ms<<" vert_count2="<<vert_count2<<endl;
    //#endif    
    if (p==n-1){
      pseudos.push_back( pseudo2 );
#ifdef _DEBUG      
      cout << "Final-> " << ii2 << " " << ms << " " << pseudo2 << endl;
#endif      
      return;
    }else{
      Array<int,1> vert_count2(n);
      vert_count2 = vert_count1.copy();
      for (int k=0; k<4; ++k){
	int it = vert1[k];
	int j1=connected_hybrid(it,0), j2=connected_hybrid(it,1);
	if (vert_count2(j1)>=0) vert_count2(j1)+=1;
	if (vert_count2(j2)>=0) vert_count2(j2)+=1;
      }
#ifdef _DEBUG      
      cout << "    pseudo"<<(p+1)<<"="<<pseudo2<<" pseudo"<<p<<"="<< pseudo1 << endl;
#endif      
      recursive_n4(p+1, pseudo2, n, vert_count2, connected_hybrid, diag, non_zero, Dlt, Fi_Fip, pseudos);
    }
  }
}

void Determine_diagrams_of_this_topology_with_option2(const vector<int>& d_known, const vector<int>& p_known,
						      deque<Array<int,1>>& pseudos, Array<int,2>& Delta_types,
						      const Array<int,2>& diag, const vector<TinyVector<int,4>>& non_zero, const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip)
{
  //#define _DEBUG  
  pseudos.clear();
  int n = diag.extent(0);
  int two_n=2*n;
  Array<int,1> pseudo1(2*n); pseudo1 = -1;

  for (int i=0; i<d_known.size(); ++i) pseudo1(d_known[i])=p_known[i]; // This is the one we know for sure
  // Each pseudoparticle propagator at time it is connected to exactly two hybridization lines, one at the beginning and one at the end.
  // connected_hybrid(it,0) and connected_hybrid(it,1) gives index to these two hybridization lines.
  Array<int,2> connected_hybrid(2*n,2);
  connected_hybrid=-1;
  for (int idi=0; idi<n; ++idi){
    TinyVector<int,4> vert = GiveVert(diag, idi, two_n);
    for (int i=0; i<4; ++i){
      int k=vert[i];
      if (connected_hybrid(k,0)<0)
	connected_hybrid(k,0)=idi;
      else
	connected_hybrid(k,1)=idi;
    }
  }
  // In which order should we process hybridization propagators?
  // It is best to do it with those that are already fixed, because recursion will have much less to do.
  Array<int,1> vert_count1(n);
  vert_count1=0;
  for (int i=0; i<d_known.size(); ++i){
    int it = d_known[i];
    vert_count1(connected_hybrid(it,0)) += 1;
    vert_count1(connected_hybrid(it,1)) += 1;
  }
#ifdef _DEBUG  
  cout<<"connected_hybrid="<<connected_hybrid<<endl;
  cout<<"vert_count1="<<vert_count1<<endl;
#endif
  int p=0;
  //TinyVector<int,4> vert0 = GiveVert(diag, order(p), two_n);

  // which hybridization line are we now addressing? The one with largest number of set pseudoparticles...
  int idi_next = max_element(vert_count1.data(), vert_count1.data()+n)-vert_count1.data();
  vert_count1(idi_next)=-1; // We eliminate this hybridization by setting count to negative number.
  TinyVector<int,4> vert0 = GiveVert(diag, idi_next, two_n);
  
  TinyVector<int,4> mms(pseudo1(vert0[0]),pseudo1(vert0[1]),pseudo1(vert0[2]),pseudo1(vert0[3])); // here we know that mms[1] was set in the input, the rest are -1
  deque<TinyVector<int,4>> possible;
  FindPossibleDeltas(possible, mms, non_zero);
#ifdef _DEBUG
  cout<<"p="<<p<<" idi="<<idi_next<<" vert_count1="<<vert_count1<<" vert0="<<vert0<<" mms="<<mms<<" pseudo="<<pseudo1<<" possible={";
  for (int i=0; i<possible.size(); ++i) cout<< possible[i]<<";";
  cout<<"}"<<endl;
#endif  
  for (int ii1=0; ii1<possible.size(); ++ii1){
    TinyVector<int,4> ms = possible[ii1];
    int dd = Dlt(ms);
    if (dd == 0) break;
    Array<int,1> pseudo2(two_n);
    pseudo2 = pseudo1.copy();
    for (int k=0; k<4; ++k) pseudo2(vert0[k]) = ms[k];
    //#ifdef _DEBUG
    //cout<<p<<") vert="<<vert0<<" ms="<<ms<<" vert_count2="<<vert_count2<<endl;
    //#endif    
    if (n==1){
      if (vert0(0)==vert0(3) and ms[0]!=ms[3]) continue;
      if (vert0(1)==vert0(2) and ms[1]!=ms[2]) continue;
      if (dd != 0){
	pseudos.push_back( pseudo2.copy() );
	//cout << "Final-> " << ii1 << " " << ms << " " << pseudo2 << endl;
      }
    }else{
      Array<int,1> vert_count2(n);
      vert_count2 = vert_count1.copy();
      for (int k=0; k<4; ++k){
	int it = vert0[k];
	int j1=connected_hybrid(it,0), j2=connected_hybrid(it,1);
	if (vert_count2(j1)>=0) vert_count2(j1)+=1;
	if (vert_count2(j2)>=0) vert_count2(j2)+=1;
      }
      //cout<<"ii1="<<ii1<<" ms="<<ms<<" dd="<<dd<<" pseudo2="<<pseudo2<<endl;
      recursive_n4(1, pseudo2, n, vert_count2, connected_hybrid, diag, non_zero, Dlt, Fi_Fip, pseudos);
    }
  }
  Delta_types.resize(pseudos.size(),n);
  for (int ip=0; ip<pseudos.size(); ++ip){
    //Delta_types[ip].resize(n);
    Array<int,1>& pseudo1 = pseudos[ip];
    for (int p=0; p<n; ++p){
      TinyVector<int,4> vert = GiveVert(diag, p, 2*n);
      TinyVector<int,4> ms(pseudo1(vert[0]),pseudo1(vert[1]),pseudo1(vert[2]),pseudo1(vert[3]));
      Delta_types(ip,p) = Dlt(ms);
    }
#ifdef _DEBUG    
    cout<<"Final2->"<<ip<<" p="<<pseudo1<<" Dlt="<<Delta_types(ip,Range::all())/*<<" from Dt="<<Delta_types[ip]*/<<endl;
#endif    
  }
}
