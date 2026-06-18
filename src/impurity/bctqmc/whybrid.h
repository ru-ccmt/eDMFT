auto WhichDelta(int itw, const Array<int,2>& diag_ci){
  struct retVals {        
    int isie;
    int id_which;   // index of the bath function, which can be between [0,Norder] and gives bath as diag_ci[id_which]
  };
  int Norder = diag_ci.extent(0);
  for (int id_w=0; id_w<Norder; ++id_w){
    // itw is between [0,2*Norder-1], i.e., any vertex that we are changing
    if (itw == diag_ci(id_w,0) or itw == diag_ci(id_w,1)){
      if (itw == diag_ci(id_w,1))
	return retVals {1, id_w};
      else
	return retVals {0, id_w};
    }
  }
  return retVals {0,-1};// should not happen
}
auto WhatDelta(double t_new, double beta, const Array<double,1>& times, int isie, const Array<int,1>& dlttyp_ci, const Array<int,2>& diag_ci, int id_whichi){
  struct retVals {
    double dt;
    int ie, is, spin, sgn;
  };
  int ad = abs(dlttyp_ci(id_whichi));
  int is = diag_ci(id_whichi,0), ie = diag_ci(id_whichi,1);//@@@
  double dt=0;
  // and time difference
  if (isie==1)           //@@@
    dt = t_new-times(is);
  else
    dt = times(ie)-t_new;
  int spin = static_cast<int>((ad-1)/2); // 1,2==up, 3,4==dn
  if (dt < 0){
    return retVals {dt+beta, ie, is, spin, -1};
  }else{
    return retVals {dt, ie, is, spin,  1};
  }
}
bool FlipDeltaLineIfNeeded(Array<int,1>& dlttyp_ci, Array<int,2>& diag_ci)
{
  bool turn=false;
  int Norder = diag_ci.extent(0);
  for (int idi=0; idi<Norder; ++idi){
    int ad = dlttyp_ci(idi);
    int sgn = 1;
    if (dlttyp_ci(idi) < 0){
      ad *= -1;
      sgn = -1;
    }
    if (ad%2==1){ // either 1 or 3 is turned in opposite direction, hence we will turn it around
      swap(diag_ci(idi,0),diag_ci(idi,1));
      dlttyp_ci(idi) = sgn*(ad+1); // means |1|->|2| and |3|->|4|
      turn = true;
    }
  }
  return turn;
}


TinyVector<pair<double,int>,2> GenerateRandomHybridization(int Norder, const TinyVector<double,2>& rand)
{
  // add hybridization
  double x1, x2;
  int i1=0, i2=0;
  x1 = rand[0]*(2*Norder);
  x2 = rand[1]*(2*Norder);
  if (x2<x1) swap(x1,x2);
  i1 = static_cast<int>(x1+1);
  i2 = static_cast<int>(x2+1)+1;
  return TinyVector<pair<double,int>,2>(make_pair(x1,i1), make_pair(x2,i2));
}

void InsertHybridizationReal(Array<int,2>& rdiag, Array<int,1>& o2n, Array<int,1>& n2o, TinyVector<int,2>& leftover,
			     const Array<int,2>& cdiag, const TinyVector<double,2>& rand)
{
  //struct retVals{
  //Array<int,2> diag;
  //Array<int,1> o2n;
  //Array<int,1> n2o;
  //retVals(int n) : diag(n+1,2), o2n(2*n), n2o(2*(n+1)){n2o=-1;};
  //};
  //retVals r(Norder);
  int Norder = cdiag.extent(0);
  rdiag.resize(Norder+1,2);
  o2n.resize(2*Norder);
  n2o.resize(2*(Norder+1));
  n2o=-1;
  Array<int,2> ndiag(Norder,2); ndiag=0;
  int i0 = rand[0]*(2*Norder);
  Array<int,1> o2n_1(2*Norder+1); o2n_1=-1;
  
  for (int i=0; i<Norder; ++i){
    for (int j=0; j<2; ++j){
      int q0=cdiag(i,j), p0; p0=q0;
      if (q0>=i0) p0+=1;
      o2n_1(p0)=q0;
      ndiag(i,j)=p0;
    }
  }
  int i1 = rand[1]*(2*Norder+1);
  for (int i=0; i<Norder; ++i){
    for (int j=0; j<2; ++j){
      int q0=ndiag(i,j), p0; p0=q0;
      if (q0>=i1) p0+=1;
      int s0=o2n_1(q0);
      if (s0>=0) o2n(s0)=p0;
      rdiag(i,j)=p0;
    }
  }
  for (int i=0; i<2*Norder; ++i) n2o(o2n(i))=i;
  
  if (i0>=i1){
    i0+=1;
    swap(i0,i1);
  }
  rdiag(Norder,0) = i0;
  rdiag(Norder,1) = i1;
  leftover = i0,i1;
  //return r;
}
void RemoveHybridizationReal(Array<int,2>& rdiag, Array<int,1>& o2n, Array<int,1>& n2o, TinyVector<int,2>& leftover,
			     const Array<int,2>& cdiag, int iremove)
{
  //struct retVals{
  //  Array<int,2> diag;
  //  Array<int,1> o2n;
  //  Array<int,1> n2o;
  //  TinyVector<int,2> leftover;
  //  retVals(int n) : diag(n-1,2), o2n(2*n), n2o(2*(n-1)) {o2n=-1;}
  //};
  //retVals r(Norder);
  int Norder = cdiag.extent(0);
  rdiag.resize(Norder-1,2);
  o2n.resize(2*Norder); o2n=-1;
  n2o.resize(2*(Norder-1));
  if (iremove>0)
    rdiag(Range(0,iremove-1),Range::all()) = cdiag(Range(0,iremove-1),Range::all());
  rdiag(Range(iremove,Norder-2),Range::all()) = cdiag(Range(iremove+1,Norder-1),Range::all());
  
  if (cdiag(iremove,0)<cdiag(iremove,1))
    leftover = cdiag(iremove,0),cdiag(iremove,1);
  else
    leftover = cdiag(iremove,1),cdiag(iremove,0);
  
  for (int j=0; j<2; ++j){
    for (int i=0; i<Norder-1; ++i){
      int q0 = rdiag(i,j), p0; p0=q0;
      if (q0 >= leftover[1]) p0 = q0-2;
      else if (q0 >= leftover[0]) p0 = q0-1;
      rdiag(i,j) = p0;
      o2n(q0) = p0;
      n2o(p0) = q0;
    }
  }
  //return r;
}

Array<int,2> InsertHybridization(const Array<int,2>& diag,const Array<double,1>& times1,Array<double,1>& times2,Array<int,1>& o2n,TinyVector<int,2>& leftover,const TinyVector<double,2>& rand,double beta)
{
  int Norder = diag.extent(0);
  o2n.resize(2*Norder);
  Array<int,2> ndiag(Norder+1,2);
  Array<double,1> _times1_(times1(Range(0,2*Norder)));
  TinyVector<double,2> tn(times1(0)*(1-rand[0]) + (times1(0)+beta)*rand[0], times1(0)*(1-rand[1]) + (times1(0)+beta)*rand[1]); // tn[0]&tn[1] are somewhere bwteen times1(0) and times1(0)+beta
  //TinyVector<double,2> tn(rand[0]*beta, rand[1]*beta); // two new times
  if (tn[1]<tn[0]) swap(tn[0],tn[1]);
  int i0 = find(tn[0], _times1_);
  if (tn[0]>times1(0)) i0+=1;
  int i1 = find(tn[1], _times1_)+1;
  if (tn[1]>times1(0)) i1+=1;
  //cout<<"tn="<<tn<<" _times1_="<<_times1_<<" i0="<<i0<<" i1="<<i1<<endl;
  for (int i=0; i<Norder; ++i){
    for (int j=0; j<2; ++j){
      int q0 = diag(i,j), p0; p0=q0;
      double tij = times1(q0);
      if (tij >= tn[0] and tij >= tn[1]) p0 += 2;
      else if (tij >= tn[0] or tij >= tn[1]) p0 += 1;
      o2n(q0)=p0;
      ndiag(i,j) = p0;
    }
  }
  ndiag(Norder,0) = i0;
  ndiag(Norder,1) = i1;
  leftover = i0,i1;
  
  for (int i=0; i<2*Norder; ++i) times2(o2n(i)) = times1(i);
  times2(i0) = tn[0];
  times2(i1) = tn[1];
  times2(2*(Norder+1)) = times2(0)+beta;
  //cout<<" tn="<<tn<<" ndiag="<<ndiag<<" o2n="<<o2n<<endl;
  return ndiag;
}
Array<int,2> RemoveHybridization(const Array<int,2>& diag,const Array<double,1>& times1,Array<double,1>& times2,Array<int,1>& o2n,TinyVector<int,2>& leftover,int ii,double beta)
{
  int Norder = diag.extent(0);
  o2n.resize(2*Norder); o2n=-1;
  Array<int,2> ndiag(Norder-1,2);
  if (ii>0){
    ndiag(Range(0,ii-1),Range::all()) = diag(Range(0,ii-1),Range::all());
  }
  ndiag(Range(ii,Norder-2),Range::all()) = diag(Range(ii+1,Norder-1),Range::all());
  
  if (diag(ii,0)<diag(ii,1))
    leftover = diag(ii,0),diag(ii,1);
  else
    leftover = diag(ii,1),diag(ii,0);
  
  //TinyVector<int,2> qs( max(diag(ii,0),diag(ii,1)), min(diag(ii,0),diag(ii,1)) );
  //TinyVector<int,2> qs = leftover[1],leftover[0];
  for (int j=0; j<2; ++j){
    for (int i=0; i<Norder-1; ++i){
      int q0 = ndiag(i,j), p0; p0=q0;
      if (q0 >= leftover[1]) p0 = q0-2;
      else if (q0 >= leftover[0]) p0 = q0-1;
      ndiag(i,j) = p0;
      o2n(q0) = p0;
    }
  }
  for (int i=0; i<2*Norder; ++i){
    int j = o2n(i);
    if (j>=0) times2(j) = times1(i);
  }
  times2(2*(Norder-1)) = times2(0)+beta;
  //cout<<" iremove="<<ii<<" ndiag="<<ndiag<<" o2n="<<o2n<<endl;
  return ndiag;
}

Array<int,2> InsertHybridization(const Array<int,2>& diag, const TinyVector<pair<double,int>,2>& xi, Array<int,1>& old2new)
{
  int Norder = diag.extent(0);
  old2new.resize(2*Norder);
  Array<int,2> ndiag(Norder+1,2);
  for (int i=0; i<Norder; ++i){
    for (int j=0; j<2; ++j){
      int q0 = diag(i,j), p0; p0=q0;
      if (p0 >= xi[0].first and p0 >= xi[1].first)
        p0 += 2;
      else if (p0 >= xi[0].first or p0 >= xi[1].first)
        p0 += 1;
      old2new(q0)=p0;
      ndiag(i,j) = p0;
    }
  }
  ndiag(Norder,0) = min(xi[0].second,xi[1].second);
  ndiag(Norder,1) = max(xi[0].second,xi[1].second);
  return ndiag;
}


Array<int,2> RemoveHybridization(const Array<int,2>& diag, int ii, Array<int,1>& old2new)
{
  int Norder = diag.extent(0);
  old2new.resize(2*Norder); old2new=-1;
  Array<int,2> ndiag(Norder-1,2);
  if (ii>0){
    ndiag(Range(0,ii-1),Range::all()) = diag(Range(0,ii-1),Range::all());
  }
  ndiag(Range(ii,Norder-2),Range::all()) = diag(Range(ii+1,Norder-1),Range::all());
  TinyVector<int,2> qs( max(diag(ii,0),diag(ii,1)), min(diag(ii,0),diag(ii,1)) );
  for (int j=0; j<2; ++j){
    for (int i=0; i<Norder-1; ++i){
      int q0 = ndiag(i,j), p0; p0=q0;
      if (q0 >= qs[0]) p0 = q0-2;
      else if (q0 >= qs[1]) p0 = q0-1;
      ndiag(i,j) = p0;
      old2new(q0) = p0;
    }
  }
  return ndiag;
}

vector<TinyVector<int,4>> Find_nonzero(const Array<int,2>& Fi, const Array<double,2>& Fm, Array<int,4>& Dlt){
  //   m3 -----> m4                                m3 -----> m4
  //         |  D[s]                                     ^  D[s]
  //         v                                           |
  //   m2 <----- m1                                m2 <----- m1
  //    Delta(m1,m2,m3,m4) = +- 1 for up           Delta(m1,m2,m3,m4) = +-2 for up
  //    Delta(m1,m2,m3,m4) = +- 3 for dn           Delta(m1,m2,m3,m4) = +-4 for dn
  //    c_s^+|m1>=+-|m2>                             c_s  |m1> =+-|m2>
  //    c_s  |m3>=+-|m4>                             c_s^+|m3> =+-|m4>
  //
  int Nst = Fi.extent(1)-1;
  deque<TinyVector<int,4>> non_zero;
  for (int s=0; s<2; ++s){
    for (int m1=0; m1<Nst; ++m1){
      int m2 = Fi(2*s,m1+1)-1;
      if (m2<0) continue;
      for (int m3=0; m3<Nst; ++m3){
	int m4 = Fi(2*s+1,m3+1)-1;
	if (m4<0) continue;
	non_zero.push_back( TinyVector<int,4>(m1,m2,m3,m4) );
	Dlt(m1,m2,m3,m4) = (1+2*s) * Fm(2*s,m1+1)*Fm(2*s+1,m3+1);
	//cout << "Adding + s="<<s << " m1=" <<m1<<" m2="<<m2<<" m3="<<m3<<" m4="<<m4<<" Dlt="<<Dlt[ii]<<endl;
      }
    }
    for (int m1=0; m1<Nst; ++m1){
      int m2 = Fi(2*s+1,m1+1)-1;
      if (m2<0) continue;
      for (int m3=0; m3<Nst; ++m3){
	int m4 = Fi(2*s,m3+1)-1;
	if (m4<0) continue;
	non_zero.push_back( TinyVector<int,4>(m1,m2,m3,m4) );
	Dlt(m1,m2,m3,m4) = (2+2*s)*Fm(2*s+1,m1+1)*Fm(2*s,m3+1);
	//cout << "Adding - s="<<s<<" m1=" <<m1<<" m2="<<m2<<" m3="<<m3<<" m4="<<m4<<" Dlt="<<Dlt[ii]<<endl;
      }
    }
  }
  vector<TinyVector<int,4>> tmp(non_zero.size());
  for (int i=0; i<tmp.size(); ++i) tmp[i] = non_zero[i];
  return tmp;
}
vector<vector<TinyVector<int,2>>> Find_nonzero2(const Array<int,2>& Fi, const Array<double,2>& Fm){
  //   m3 -----> m4                                m3 -----> m4
  //         |  D[s]                                     ^  D[s]
  //         v                                           |
  //   m2 <----- m1                                m2 <----- m1
  //    Delta(m1,m2,m3,m4) = +- 1 for up           Delta(m1,m2,m3,m4) = +-2 for up
  //    Delta(m1,m2,m3,m4) = +- 3 for dn           Delta(m1,m2,m3,m4) = +-4 for dn
  //    c_s^+|m1>=+-|m2>                             c_s  |m1> =+-|m2>
  //    c_s  |m3>=+-|m4>                             c_s^+|m3> =+-|m4>
  int Nst = Fi.extent(1)-1;
  vector<deque<TinyVector<int,2>>> non_zero2(Fi.extent(0));
  for (int s=0; s<Fi.extent(0); ++s){
    for (int m1=0; m1<Nst; ++m1){
      int m2 = Fi(s,m1+1)-1;
      if (m2<0) continue;
      non_zero2[s].push_back( TinyVector<int,2>(m1,m2) );
    }
  }
  vector<vector<TinyVector<int,2>>> tmp(non_zero2.size());
  for (int i=0; i<tmp.size(); ++i){
    tmp[i].resize(non_zero2[i].size());
    for (int j=0; j<tmp[i].size(); ++j){
      tmp[i][j] = non_zero2[i][j];
    }
  }
  return tmp;
}
