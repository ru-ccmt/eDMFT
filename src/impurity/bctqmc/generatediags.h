
//using namespace std;
//using namespace std::chrono;
//using namespace blitz;
//using json = nlohmann::json;
void copy(deque<Array<int,1>>& pseudos_c, const deque<Array<int,1>>& pseudos_t)
{
  pseudos_c.resize(pseudos_t.size());
  for (int i=0; i<pseudos_t.size(); ++i){
    if (pseudos_c[i].size()!=pseudos_t[i].size()) pseudos_c[i].resize(pseudos_t[i].size());
    pseudos_c[i] = pseudos_t[i];
  }
}
void copy(vector<Array<int,1>>& Delta_types_c, const vector<Array<int,1>>& Delta_types_t)
{
  Delta_types_c.resize(Delta_types_t.size());
  for (int i=0; i<Delta_types_t.size(); ++i){
    if (Delta_types_c[i].size()!=Delta_types_t[i].size()) Delta_types_c[i].resize(Delta_types_t[i].size());
    Delta_types_c[i] = Delta_types_t[i];
  }
}
template<class RNG>
void GenerateDiagrams2(Array<int,2>& i_diag,const vector<TinyVector<int,4>> non_zero,const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip, RNG& drand,int Ndiag=10000){
  //#define _DEBUG  
  cout.precision(12);
  int N_min_order=1, N_max_order=i_diag.extent(0);

  Array<double,1> Vol_diagrams0(max(N_max_order+1,8)); Vol_diagrams0=1;
  Vol_diagrams0 = 1,4,2,16./3.,27,198.4,1886.67,21846.607;
  double P_first_jump = 0.5;

  deque<Array<int,1>> pseudos_c0, pseudos_trial;
  Array<int,2> Delta_types_c0, Delta_types_trial;
  
  int Norder=i_diag.extent(0);
  Array<int,2> diag(Norder,2), cdiag(Norder,2), ndiag(Norder,2);
  diag = i_diag;
  Skeletons skl(diag);
  bool Reducible = skl.IsReducible();
  {
    int sgn_dlt=1;
    FindNextSkeleton(sgn_dlt, diag, skl,drand);
  }
  Skeletons tskl(diag);
  
  TinyVector<Array<int,2>,4> diag_c;
  TinyVector<Array<int,1>,4> pseud_c;
  TinyVector<Array<int,1>,4> dlttyp_c;
  TinyVector<int,4> sgn_t;
  for (int i=0; i<4; ++i){
    diag_c[i].resize(Norder,2);
    pseud_c[i].resize(2*Norder);
    dlttyp_c[i].resize(Norder);
  }
  Array<int,1> cpseudo(2*Norder), dpseudo(2*Norder), dpseudo2(2*Norder), cDlt_type(Norder);
  
  Determine_all_diagrams_of_this_topology2(pseudos_c0, Delta_types_c0, diag, non_zero, Dlt, Fi_Fip);
  int iwhich_c0 = static_cast<int>(drand()*pseudos_c0.size());
  cpseudo = pseudos_c0[iwhich_c0].copy();
  cDlt_type = Delta_types_c0(iwhich_c0,Range::all());
  
  Array<int,1> how_many_equal(pseudos_c0.size()); how_many_equal = 0;
  for (int ip=0; ip<pseudos_c0.size(); ++ip){
    dpseudo = pseudos_c0[ip]-cpseudo;
    how_many_equal(ip) = count_zeros( dpseudo );
    cout<<ip<<" p="<< pseudos_c0[ip] << " has " << how_many_equal(ip) << " equal  cDlt=" << cDlt_type << endl;
  }
  cout << "Starting with diag=" <<diag<< " and p=" << cpseudo <<" cDlt=" << cDlt_type << endl;
  cdiag = diag;
  int csgn_dlt = SignSort(cdiag);

  int idiag=0;
  double Pnorm=0, Pbin=0, Pnorm_sum=0, Pbin_sum=0;
  Array<double,1> Nmeassure(N_max_order+1), asign(N_max_order+1);
  Nmeassure=0; asign=0;
  //double Pchg=1, bt=1.0;
  map<string,double> stored2;
  deque<pair<Array<int,2>,Array<int,1>>> stored;
  deque<double> stored_count;
  deque<double> stored_sign;
  Array<int,1> o2n;
  int tNorder=Norder;
  tskl = skl;
  Array<double,1> dSm(4); dSm=0;
  double V0norm=1.0;
  TinyVector<double,4> P_acc, P_pr, P_eq;
  double P_tree=0;
  int iacc=0;
  int iwhich_t1, iwhich_t2, iwhich_t3;
  fsec t_next=zero_t, t_diff=zero_t, t_choose=zero_t, t_topol=zero_t, t_skel=zero_t;
  double how_many_pgp=0;
  for (int ii=0; ii<100000000; ++ii){
    auto t_1 = high_resolution_clock::now();
    diag_c[0] = cdiag;
    pseud_c[0] = cpseudo;
    dlttyp_c[0] = cDlt_type;
    sgn_t[0] = csgn_dlt;
    
    diag_c[1] = cdiag;
    how_many_equal.resize(pseudos_c0.size());
    int ic_to_compare = 0;
    for (int ip=0; ip<pseudos_c0.size(); ++ip){
      if (ip==iwhich_c0)
	how_many_equal(ip) = 0; // we should avoid the current diagram, because we want to change
      else{
	dpseudo = pseudos_c0[ip]-pseud_c[ic_to_compare];
	how_many_equal(ip) = count_zeros( dpseudo );
      }
#ifdef _DEBUG      
      cout<<ip<<" p="<< pseudos_c0[ip] << " has " << how_many_equal(ip) << " equal with po=" << pseud_c[ic_to_compare]<<endl;
#endif      
    }
    deque<int> candidates1 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates1 << endl;
#endif    
    iwhich_t1 = candidates1[ static_cast<int>(candidates1.size()*drand()) ];
    pseud_c[1] = pseudos_c0[iwhich_t1];
    dlttyp_c[1] = Delta_types_c0(iwhich_t1,Range::all());
    sgn_t[1] = csgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t1=" << iwhich_t1 << " p_c1=" << pseud_c[1] << endl;
#endif	  
    Norder = cdiag.extent(0);
    int sgn_dlt = csgn_dlt;
    int jump=0;
    double xi=drand();
    if (xi < 1./3.) jump = 1;
    else if ( xi < 2./3.) jump = -1;
    if (jump==1 and Norder >= N_max_order) jump=0;
    if (jump==-1 and Norder <= N_min_order) jump=0;
#ifdef _DEBUG    
    cout<<ii<<" diag="<<diag<<" cdiag="<<cdiag<<" jump="<<jump<<endl;
#endif    
    if (jump==1){
      tNorder = Norder+1;
      TinyVector<pair<double,int>,2>  hybn = GenerateRandomHybridization(Norder,TinyVector<double,2>(drand(),drand()));
      diag.resize(tNorder,2);
      diag = InsertHybridization(cdiag,hybn,o2n);
      tskl.Resize_Initialize(diag);
      if (tskl.IsReducible()) FindNextSkeleton(sgn_dlt, diag, tskl,drand);
      sgn_dlt = SignSort(diag);
      dpseudo2.resize(2*tNorder);
      diag_c[2].resize(tNorder,2);
      pseud_c[2].resize(2*tNorder);
      dlttyp_c[2].resize(tNorder);
      diag_c[3].resize(tNorder,2);
      pseud_c[3].resize(2*tNorder);
      dlttyp_c[3].resize(tNorder);
    }else if (jump==-1){
      int iremove = Norder*drand();
      tNorder = Norder-1;
      diag.resize(tNorder,2);
      diag = RemoveHybridization(cdiag, iremove,o2n);
      tskl.Resize_Initialize(diag);
      if (tskl.IsReducible()) FindNextSkeleton(sgn_dlt, diag, tskl,drand);
      sgn_dlt = SignSort(diag);
      dpseudo2.resize(2*tNorder);
      diag_c[2].resize(tNorder,2);
      pseud_c[2].resize(2*tNorder);
      dlttyp_c[2].resize(tNorder);
      diag_c[3].resize(tNorder,2);
      pseud_c[3].resize(2*tNorder);
      dlttyp_c[3].resize(tNorder);
    }else{
      tNorder = Norder;
      diag = cdiag;
      FindNextSkeleton(sgn_dlt, diag, tskl,drand);
    }
    //if (jump) cout<<"diag="<<diag<<" o2n="<<o2n<<endl;
    // All diagrams with all possible pseudoparticle arrangements are produced by recursion
    auto t_2 = high_resolution_clock::now();
    t_skel += t_2-t_1; 
    Determine_all_diagrams_of_this_topology2(pseudos_trial, Delta_types_trial, diag, non_zero, Dlt, Fi_Fip);
    auto t_3 = high_resolution_clock::now();
    t_topol += t_3-t_2;
#ifdef _DEBUG    
    cout<<ii<<" Number of possible diagrams="<<pseudos_trial.size()<<endl;
#endif
    int icase=0;
    icase = (drand() < P_first_jump) ? 1 : 2;
#ifdef _DEBUG    
    cout << "icase=" << icase << " cdiag="<<cdiag<<" diag="<<diag << " iwhich_c0=" << iwhich_c0 << " p_c0=" << cpseudo << " csgn_dlt=" << csgn_dlt << " sgn_dlt=" << sgn_dlt << endl;
#endif    
    diag_c[2] = diag;
    how_many_equal.resize(pseudos_trial.size());
    ic_to_compare = (icase==1) ? 1 : 0;
    for (int ip=0; ip<pseudos_trial.size(); ++ip){
      if (jump){
	for (int j=0; j<2*Norder; ++j){
	  int jn = o2n(j);
	  if (jn>=0) dpseudo(j) = pseudos_trial[ip](jn)-pseud_c[ic_to_compare](j);
	  else dpseudo(j) = -1;
	}
      }else{
	dpseudo = pseudos_trial[ip]- pseud_c[ic_to_compare];
      }
      how_many_equal(ip) = count_zeros( dpseudo );
#ifdef _DEBUG      
      cout<<setw(2)<<ip<<" p="<< pseudos_trial[ip] << " has " << how_many_equal(ip) << " equal with po=" << pseud_c[ic_to_compare]<<endl;
#endif      
    }
    deque<int> candidates2 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates2 << endl;
#endif
    iwhich_t2 = candidates2[ static_cast<int>(candidates2.size()*drand()) ];
    pseud_c[2] = pseudos_trial[iwhich_t2];
    dlttyp_c[2] = Delta_types_trial(iwhich_t2,Range::all());
    sgn_t[2] = sgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t2=" << iwhich_t2 << " p_c2=" << pseud_c[2] << endl;
#endif
    diag_c[3] = diag;
    ic_to_compare = 2;
    for (int ip=0; ip<pseudos_trial.size(); ++ip){
      if (ip==iwhich_t2)
	how_many_equal(ip) = 0; // we should avoid the current diagram, because we want to change
      else{
	dpseudo2 = pseudos_trial[ip]-pseud_c[ic_to_compare];
	how_many_equal(ip) = count_zeros( dpseudo2 );
      }
#ifdef _DEBUG      
      cout<<setw(2)<<ip<<" p="<< pseudos_trial[ip] << " has " << how_many_equal(ip) << " equal with po=" << pseud_c[ic_to_compare]<<endl;
#endif      
    }
    deque<int> candidates3 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates3 << endl;
#endif    
    iwhich_t3 = candidates3[ static_cast<int>(candidates3.size()*drand()) ];
    pseud_c[3] = pseudos_trial[iwhich_t3];
    dlttyp_c[3] = Delta_types_trial(iwhich_t3,Range::all());
    sgn_t[3] = sgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t3=" << iwhich_t3 << " p_c3=" << pseud_c[3] << endl;
#endif	  
    auto t_4 = high_resolution_clock::now();
    t_choose += t_4-t_3;
    
    P_eq = 1.,1.,1.,1.;
    
    double V_0=pseudos_c0.size()*Vol_diagrams0(Norder), V_1=pseudos_trial.size()*Vol_diagrams0(tNorder);
    if (jump){
      cout<<" jump Norder="<<Norder<<"->"<<tNorder<<" V_0="<<V_0<<" V_1="<<V_1<<" N0="<<pseudos_c0.size()<<" N1="<<pseudos_trial.size()<<endl;
    }
    //P_pr[0] = 1.0/(V_0*V_1*V_1);
    //P_pr[1] = 1.0/(V_0*V_1*V_1);
    //P_pr[2] = 1.0/(V_0*V_0*V_1);
    //P_pr[3] = 1.0/(V_0*V_0*V_1);
    P_pr[0] = 0.5*V_0/(V_0+V_1);
    P_pr[1] = P_pr[0];
    P_pr[2] = 0.5*V_1/(V_0+V_1);
    P_pr[3] = P_pr[2];
    P_tree=0;
    for (int ic=0; ic<4; ++ic){
      double P0 = P_eq[ic]*P_pr[ic];
      P_tree += P0;
      P_acc[ic] = P0;
    }
    P_acc *= 1/P_tree;
#ifdef _DEBUG    
    cout << "P_tree="<<P_tree<<" P_acc=" << P_acc <<endl;
#endif
    //double _asign_ = 0;
    for (int ic=0; ic<4; ++ic){
      //TinyVector<int,2> sgn_mmbw = SignMatrixElements(diag_c[ic], pseud_c[ic], Dlt);
      //sgn_t[ic] *= sgn_mmbw[0];
      //int sgn_tot = sgn_t[ic]*SgnMElements(diag_c[ic], pseud_c[ic], Dlt);
      //_asign_ += sgn_tot*P_acc[ic];
      sgn_t[ic] *=  SgnMElements(diag_c[ic], pseud_c[ic], Dlt);
      auto [new_diag,s] = Different3(diag_c[ic], pseud_c[ic], P_acc[ic], stored2);
      cout <<ii<<") Final-Diag "<< diag_c[ic] << " with p=" << pseud_c[ic] << " is equivalent to " << s << " wgh=" << P_acc[ic] << " sgn=" << sgn_t[ic] << endl;
    }
    auto t_5 = high_resolution_clock::now();
    t_diff += t_5-t_4;
    // meassure
    {
      for (int ic=0; ic<4; ++ic){
	int norder = diag_c[ic].extent(0);
	Nmeassure(norder) += P_acc[ic];
	asign(norder) += sgn_t[ic]*P_acc[ic];
	for (int it=0; it<pseud_c[ic].size(); ++it){ //### Norder
	  int atomic_state = pseud_c[ic](it);
	  dSm(atomic_state) += P_acc[ic];
	}
      }
    }
    double xir = drand();
#ifdef _DEBUG    
    cout << "xir=" << xir << endl;
#endif    
    if (xir > 1-P_acc[3] ){ //@@
      // accept c3
      iacc = 3;
      if (jump!=0){
	cdiag.resize(tNorder,2);
	cpseudo.resize(2*tNorder);
	cDlt_type.resize(tNorder);
      }
      cdiag = diag_c[3].copy();
      cpseudo = pseud_c[3].copy();
      cDlt_type = dlttyp_c[3].copy();
      iwhich_c0 = iwhich_t3;
      copy(pseudos_c0,pseudos_trial); pseudos_trial.clear();

      if (Delta_types_c0.extent(0)!=Delta_types_trial.extent(0) or
	  Delta_types_c0.extent(1)!=Delta_types_trial.extent(1) ) Delta_types_c0.resize(Delta_types_trial.shape());
      Delta_types_c0 = Delta_types_trial;
      
      Delta_types_trial.resize(0,0);
      csgn_dlt = sgn_dlt;
      skl = tskl;
      cout << ii<<") Choosing c3 icase=" << icase << " p=" << cpseudo << " diag=" <<cdiag<< " iwhich=" << iwhich_c0 << " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      if (jump!=0){
	Norder = tNorder;
	dpseudo.resize(2*Norder);
	diag_c[0].resize(Norder,2);
	pseud_c[0].resize(2*Norder);
	dlttyp_c[0].resize(Norder);
	diag_c[1].resize(Norder,2);
	pseud_c[1].resize(2*Norder);
	dlttyp_c[1].resize(Norder);	  
      }
    }else if (xir <= 1-P_acc[3] and xir > P_acc[0]+P_acc[1]){
      // accept c2
      iacc = 2;
      if (jump!=0){
	cdiag.resize(tNorder,2);
	cpseudo.resize(2*tNorder);
	cDlt_type.resize(tNorder);
      }
      cdiag = diag_c[2].copy();
      cpseudo = pseud_c[2].copy();
      cDlt_type = dlttyp_c[2].copy();
      iwhich_c0 = iwhich_t2;
      copy(pseudos_c0,pseudos_trial); pseudos_trial.clear();

      if (Delta_types_c0.extent(0)!=Delta_types_trial.extent(0) or
	  Delta_types_c0.extent(1)!=Delta_types_trial.extent(1) ) Delta_types_c0.resize(Delta_types_trial.shape());
      Delta_types_c0 = Delta_types_trial;
      Delta_types_trial.resize(0,0);
      //copy(Delta_types_c0,Delta_types_trial); Delta_types_trial.clear();
      
      csgn_dlt = sgn_dlt;
      skl = tskl;
      cout <<ii<<") Choosing c2 icase=" << icase << " p=" << cpseudo << " diag=" <<cdiag<< " iwhich=" << iwhich_c0 << " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      if (jump!=0){
	Norder = tNorder;
	dpseudo.resize(2*Norder);
	diag_c[0].resize(Norder,2);
	pseud_c[0].resize(2*Norder);
	dlttyp_c[0].resize(Norder);
	diag_c[1].resize(Norder,2);
	pseud_c[1].resize(2*Norder);
	dlttyp_c[1].resize(Norder);	  
      }
    }else if (xir <= P_acc[0]+P_acc[1] and xir > P_acc[0]){
      iacc = 1;
      cpseudo = pseud_c[1].copy();
      cDlt_type = dlttyp_c[1].copy();
      iwhich_c0 = iwhich_t1;
      tskl = skl;
      cout<<ii<<") Choosing c1 icase="<<icase<<"  p="<<cpseudo<<" diag="<<cdiag;
      cout<<" iwhich="<< iwhich_c0<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      if (jump!=0){
	diag.resize(Norder,2);
	dpseudo2.resize(2*Norder);
	diag_c[2].resize(Norder,2);
	pseud_c[2].resize(2*Norder);
	dlttyp_c[2].resize(Norder);
	diag_c[3].resize(Norder,2);
	pseud_c[3].resize(2*Norder);
	dlttyp_c[3].resize(Norder);
      }
    } else if (xir <= P_acc[0]){ // c0 -> backtracking
      iacc = 0;
      cout <<ii<<") ************************* BACKTRACKING ***************" << endl;
      //tskl = skl;
      tskl.copy(skl);
      if (jump!=0){
	diag.resize(Norder,2);
	dpseudo2.resize(2*Norder);
	diag_c[2].resize(Norder,2);
	pseud_c[2].resize(2*Norder);
	dlttyp_c[2].resize(Norder);
	diag_c[3].resize(Norder,2);
	pseud_c[3].resize(2*Norder);
	dlttyp_c[3].resize(Norder);
      }
    }
    idiag += 1;
    auto t_6 = high_resolution_clock::now();
    t_next += t_6-t_5;
    if (idiag >= Ndiag) break;
  }
  //cout << "how_many_pgp=" << how_many_pgp/Ndiag/(2*Norder) << endl;

  for (auto itt=stored2.begin(); itt!=stored2.end(); ++itt){
    auto [diag,pseud] = String2diagp(itt->first);
    stored.push_back(make_pair(diag,pseud));
    stored_count.push_back( itt->second );
    cout<<itt->first<<" "<<itt->second<<" diag="<<diag<<" pseud="<<pseud<<endl;
  }
  
  vector<int> sfc(stored.size());
  for (int i=0; i<stored.size(); ++i)
    sfc[i] = SymmetryFactor2(stored[i]);

  deque<Array<int,1>> pseudos;
  Array<int,2> Delta_types;
  diag.resize(1,2); diag = 1,0;
  Determine_all_diagrams_of_this_topology(pseudos,Delta_types, diag,non_zero,Dlt);
  double N00 = pseudos.size();
  Array<double,1> Vol_diagrams(N_max_order+1), Vol_diagrams_jump(N_max_order+1);
  Vol_diagrams=0; Vol_diagrams_jump=0;
  for (int i=0; i<stored.size(); ++i){
    Array<int,2> diag(stored[i].first);
    int Norder = diag.extent(0);
    Vol_diagrams(Norder) += 1.0/sfc[i];
    Determine_all_diagrams_of_this_topology(pseudos,Delta_types, diag,non_zero,Dlt);
    Vol_diagrams_jump(Norder) += N00/(sfc[i]*pseudos.size());
  }
  Array<double,1> rr(N_max_order+1); rr = 0;
  Array<int,1> ndg(N_max_order+1); ndg=0;
  double rrt=0, rrv=0;
  for (int i=0; i<stored.size(); ++i){
    int Norder = stored[i].first.extent(0);
    ndg(Norder)  += 1;
    double xt = sfc[i]*stored_count[i];
    rr(Norder) += xt;
    rrt += xt;
    rrv += xt*Vol_diagrams0(Norder);
  }
  for (int i=N_min_order; i<=N_max_order; ++i) rr(i) /= ndg(i);
  rrt /= stored.size();
  rrv /= stored.size();
  cout<<"rr="<<rr<<" rrt="<<rrt<<" rrv="<<rrv<<" sfc="<<sfc<<endl;

  Array<double,1> w_order(N_max_order+1); w_order=0;
  Array<double,1> deviation1(N_max_order+1); deviation1=0;
  double deviation3=0, deviation2=0;
  cout <<"stored diags="<<endl;
  for (int iorder=N_min_order; iorder<=N_max_order; ++iorder){
    for (int i=0; i<stored.size(); ++i){
      int Norder = stored[i].first.extent(0);
      if (iorder!=Norder) continue;
      double xt = stored_count[i]*sfc[i];
      double x1 = xt/rr(Norder);
      double x2 = xt/rrt;
      double yt = xt*Vol_diagrams(Norder);
      double x3 = yt/rrv;
      w_order(Norder) += x2;
      deviation1(Norder) += abs(x1-1);
      deviation2 += abs(x2-1);
      deviation3 += abs(x3-1);

      cdiag = stored[i].first;
      cpseudo = stored[i].second;
      /*
	int sgn_dlt = SignSort(cdiag);
	TinyVector<int,2> sgn_mmbw = SignMatrixElements(cdiag, cpseudo, Dlt);
	int sgn_tot = sgn_dlt*sgn_mmbw[0];
	int sgn_cur = sgn_dlt*sgn_mmbw[0]*sgn_mmbw[1];
      */
      cout<<setw(3)<<i<<" "<<setw(2)<<Norder<<" "<<setw(8)<<setprecision(3)<<xt<<" ";
      cout<<setw(2)<<sfc[i]<<"  ";
      cout<<left<<setw(8)<<setprecision(6)<<x1<<" "<<setw(8)<<x2<<" "<<setw(8)<<x3<<" "<<right;
      cout<<"d="<<stored[i].first<<" p=" << stored[i].second<<endl;
    }
  }
  for (int i=N_min_order; i<=N_max_order; ++i) w_order(i) /= ndg(i);
  for (int i=N_min_order; i<=N_max_order; ++i) cout<<"w["<<i<<"]=" <<w_order(i) << endl;
  
  deviation3 /= stored.size();
  deviation2 /= stored.size();
  cout<<"Vol_diagrams="<<Vol_diagrams<<endl;
  cout<<"Vol_diagrams jump="<<Vol_diagrams_jump<<endl;
  
  for (int i=N_min_order; i<=N_max_order; ++i) deviation1(i) /= ndg(i);
  cout<<"deviation1=[";
  for (int i=N_min_order; i<=N_max_order; ++i) cout<<deviation1(i)<<",";
  cout<<"] ";
  cout << " deviation2="<<deviation2<<" deviation3=" << deviation3 << endl;

  cout<<"<s>=[";
  for (int iorder=N_min_order; iorder<=N_max_order; ++iorder) cout<< asign(iorder)/Nmeassure(iorder) <<",";
  cout<<"]"<<endl;

  /*
  cout << "Before: Pbin="<<Pbin*(V0norm/Pnorm)<<" dSm="<< Array<double,1>((dSm*(V0norm/Pnorm))) <<endl;
  Pbin *= (V0norm/Pnorm)*Vol_diagrams;
  dSm  *= (V0norm/Pnorm)*Vol_diagrams;
  cout << "After Pbin="<< Pbin <<" dSm="<<dSm<<endl;
  */
  cout<<"t_next="<<t_next.count()<<" t_diff="<<t_diff.count()<<" t_choose="<<t_choose.count()<<" t_topol="<<t_topol.count()<<" t_skel="<<t_skel.count()<<endl;
}

template<class RNG>
void GenerateDiagrams(Array<int,2>& i_diag,const vector<TinyVector<int,4>> non_zero,const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip, RNG& drand,int Ndiag=10000){
  cout.precision(12);
  int Norder=i_diag.extent(0);
  Array<int,2> diag(Norder,2), cdiag(Norder,2), ndiag(Norder,2);
  diag = i_diag;
  double P_first_jump = 0.5;

  Skeletons skl(diag);
  bool Reducible = skl.IsReducible();
  {
    int sgn_dlt=1;
    FindNextSkeleton(sgn_dlt, diag, skl,drand);
  }

  TinyVector<Array<int,2>,4> diag_c;
  TinyVector<Array<int,1>,4> pseud_c;
  TinyVector<Array<int,1>,4> dlttyp_c;
  TinyVector<int,4> sgn_t;
  for (int i=0; i<4; ++i){
    diag_c[i].resize(Norder,2);
    pseud_c[i].resize(2*Norder);
    dlttyp_c[i].resize(Norder);
  }
  
  deque<Array<int,1>> pseudos_c0, pseudos_trial;
  Array<int,2> Delta_types_c0, Delta_types_trial;
  Array<int,1> cpseudo(2*Norder), dpseudo(2*Norder), cDlt_type(Norder);
  
  Skeletons skl_copy(diag);
  Determine_all_diagrams_of_this_topology2(pseudos_c0, Delta_types_c0, diag, non_zero, Dlt, Fi_Fip);
  int iwhich_c0 = static_cast<int>(drand()*pseudos_c0.size());
  cpseudo = pseudos_c0[iwhich_c0].copy();
  cDlt_type = Delta_types_c0(iwhich_c0,Range::all());
    

  Array<int,1> how_many_equal(pseudos_c0.size()); how_many_equal = 0;
  for (int ip=0; ip<pseudos_c0.size(); ++ip){
    dpseudo = pseudos_c0[ip]-cpseudo;
    how_many_equal(ip) = count_zeros( dpseudo );
    cout<<ip<<" p="<< pseudos_c0[ip] << " has " << how_many_equal(ip) << " equal  cDlt=" << cDlt_type << endl;
  }
  cout << "Starting with diag=" << diag << " and p=" << cpseudo <<" cDlt=" << cDlt_type << endl;

  cdiag = diag;
  int csgn_dlt = SignSort(cdiag);

  int idiag=0;
  double Pnorm=0, Pbin=0, Pnorm_sum=0, Pbin_sum=0;
  int Nmeassure=0;

  //bool Qstore=false;
  map<string,double> stored2;
  deque<pair<Array<int,2>,Array<int,1>>> stored;
  deque<double> stored_count;
  deque<double> stored_sign;
  
  Array<double,1> dSm(4); dSm=0;
  double V0norm=1.0;
  TinyVector<double,4> P_acc, P_pr, P_eq;
  TinyVector<double,2> P_tree; P_tree = 0,V0norm;
  int iacc=0;
  int iwhich_t1, iwhich_t2, iwhich_t3;
  fsec t_next=zero_t, t_diff=zero_t, t_choose=zero_t, t_topol=zero_t, t_skel=zero_t;
  double how_many_pgp=0;
  for (int ii=0; ii<10000000; ++ii){
    auto t_1 = high_resolution_clock::now();
    diag_c[0] = cdiag;
    pseud_c[0] = cpseudo;
    dlttyp_c[0] = cDlt_type;
    sgn_t[0] = csgn_dlt;
    //
    //bool Reducible = true;
    skl_copy = skl;
    diag = cdiag;
    int sgn_dlt = csgn_dlt;
    FindNextSkeleton(sgn_dlt, diag, skl,drand);
    // All diagrams with all possible pseudoparticle arrangements are produced by recursion
    auto t_2 = high_resolution_clock::now();
    t_skel += t_2-t_1; 
    Determine_all_diagrams_of_this_topology2(pseudos_trial, Delta_types_trial, diag, non_zero, Dlt, Fi_Fip);
    auto t_3 = high_resolution_clock::now();
    t_topol += t_3-t_2;
    cout<<ii<<" Number of possible diagrams="<<pseudos_trial.size()<<endl;
      
    diag_c[1] = cdiag;
    diag_c[2] = diag;
    diag_c[3] = diag;
    int icase=0;
    icase = (drand() < P_first_jump) ? 1 : 2;
#ifdef _DEBUG    
    cout << "icase=" << icase << " cdiag="<<cdiag<<" diag="<<diag<< " iwhich_c0=" << iwhich_c0 << " p_c0=" << cpseudo << " csgn_dlt=" << csgn_dlt << " sgn_dlt=" << sgn_dlt << endl;
#endif    
    how_many_equal.resize(pseudos_c0.size());
#ifdef _DEBUG    
    cout << "pseudos_c0.size=" << pseudos_c0.size() << endl;
#endif    
    int ic_to_compare = 0;
    for (int ip=0; ip<pseudos_c0.size(); ++ip){
      if (ip==iwhich_c0)
	how_many_equal(ip) = 0; // we should avoid the current diagram, because we want to change
      else{
	dpseudo = pseudos_c0[ip]-pseud_c[ic_to_compare];
	how_many_equal(ip) = count_zeros( dpseudo );
      }
#ifdef _DEBUG      
      cout<<ip<<" p="<< pseudos_c0[ip] << " has" << how_many_equal(ip) << " equal" << endl;
#endif      
    }
    deque<int> candidates1 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates1 << endl;
#endif    
    iwhich_t1 = candidates1[ static_cast<int>(candidates1.size()*drand()) ];
    pseud_c[1] = pseudos_c0[iwhich_t1];
    dlttyp_c[1] = Delta_types_c0(iwhich_t1,Range::all());
    sgn_t[1] = csgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t1=" << iwhich_t1 << " p_c1=" << pseud_c[1] << endl;
#endif	  
    how_many_equal.resize(pseudos_trial.size());
    ic_to_compare = (icase==1) ? 1 : 0;
    for (int ip=0; ip<pseudos_trial.size(); ++ip){
      dpseudo = pseudos_trial[ip]- pseud_c[ic_to_compare];
      how_many_equal(ip) = count_zeros( dpseudo );
#ifdef _DEBUG      
      cout<<ip<<" p="<< pseudos_trial[ip] << " has" << how_many_equal(ip) << " equal" << endl;
#endif      
    }
    deque<int> candidates2 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates2 << endl;
#endif     
    iwhich_t2 = candidates2[ static_cast<int>(candidates2.size()*drand()) ];
    pseud_c[2] = pseudos_trial[iwhich_t2];
    dlttyp_c[2] = Delta_types_trial(iwhich_t2,Range::all());
    sgn_t[2] = sgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t2=" << iwhich_t2 << " p_c2=" << pseud_c[2] << endl;
#endif
    ic_to_compare = 2;
    for (int ip=0; ip<pseudos_trial.size(); ++ip){
      if (ip==iwhich_t2)
	how_many_equal(ip) = 0; // we should avoid the current diagram, because we want to change
      else{
	dpseudo = pseudos_trial[ip]-pseud_c[ic_to_compare];
	how_many_equal(ip) = count_zeros( dpseudo );
      }
#ifdef _DEBUG      
      cout<<ip<<" p="<< pseudos_trial[ip] << " has" << how_many_equal(ip) << " equal" << endl;
#endif      
    }
    deque<int> candidates3 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates3 << endl;
#endif    
    iwhich_t3 = candidates3[ static_cast<int>(candidates3.size()*drand()) ];
    pseud_c[3] = pseudos_trial[iwhich_t3];
    dlttyp_c[3] = Delta_types_trial(iwhich_t3,Range::all());
    sgn_t[3] = sgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t3=" << iwhich_t3 << " p_c3=" << pseud_c[3] << endl;
#endif	  

    auto t_4 = high_resolution_clock::now();
    t_choose += t_4-t_3;
    
    P_eq = 1.,1.,1.,1.;
    
    double N0 = pseudos_c0.size(), N1 = pseudos_trial.size();
    //P_pr[0] = 1.0/(N0*N1*N1);
    //P_pr[1] = 1.0/(N0*N1*N1);
    //P_pr[2] = 1.0/(N0*N0*N1);
    //P_pr[3] = 1.0/(N0*N0*N1);
    P_pr[0] = 0.5*N0/(N0+N1);
    P_pr[1] = 0.5*N0/(N0+N1);
    P_pr[2] = 0.5*N1/(N0+N1);
    P_pr[3] = 0.5*N1/(N0+N1);
    P_tree[0]=0;
    for (int ic=0; ic<4; ++ic){
      //double P0 = P_eq[ic]*P_pr[ic]; 
      P_tree[0] += P_eq[ic]*P_pr[ic];
      P_acc[ic] = (P_eq[ic] + P_tree[1])*P_pr[ic];
    }
    P_acc *= 1/(P_tree[0]+P_tree[1]);
    //
    for (int ic=0; ic<4; ++ic){
      TinyVector<int,2> sgn_mmbw = SignMatrixElements(diag_c[ic], pseud_c[ic], Dlt);
      sgn_t[ic] *= sgn_mmbw[0];
      auto [new_diag,s] = Different3(diag_c[ic], pseud_c[ic], P_acc[ic], stored2);
      //cout << "Final-Diag "<< diag_c[ic] << " with p=" << pseud_c[ic] << " is equivalent to " << s << " and count is " << iequiv << " wgh=" << P_acc[ic] << " sgn=" << sgn_t[ic] << endl;
      //}
    }
    auto t_5 = high_resolution_clock::now();
    t_diff += t_5-t_4;
    
    // meassure
    {
      Nmeassure += 1;
      double Wphys  = P_tree[0]/(P_tree[0]+P_tree[1]);
      double Walter = P_tree[1]/(P_tree[0]+P_tree[1]);
      Pnorm += Walter;
      Pnorm_sum += Walter;
      Pbin_sum  += Wphys;
      Pbin += Wphys;
      for (int ic=0; ic<4; ++ic){
	double partialPhi = P_eq[ic]*P_pr[ic]/(P_tree[0]+P_tree[1]);
	for (int it=0; it<2*Norder; ++it){
	  int atomic_state = pseud_c[ic](it);
	  dSm(atomic_state) += partialPhi;
	}
      }
    }
    double xir = drand();
    cout << "xir=" << xir << endl;
    if (xir > 1-P_acc[3] ){ //@@
      // accept c3
      iacc = 3;
      cpseudo = pseud_c[3].copy();
      cDlt_type = dlttyp_c[3].copy();
      cdiag = diag_c[3].copy();
      iwhich_c0 = iwhich_t3;
      pseudos_c0 = pseudos_trial; pseudos_trial.clear();
      if (Delta_types_c0.extent(0)!=Delta_types_trial.extent(0) or
	  Delta_types_c0.extent(1)!=Delta_types_trial.extent(1) ) Delta_types_c0.resize(Delta_types_trial.shape());
      Delta_types_c0 = Delta_types_trial;
      //Delta_types_trial.clear();
      Delta_types_trial.resize(0,0);
      csgn_dlt = sgn_dlt;
      cout << "Choosing c3 icase=" << icase << " p=" << cpseudo << " diag=" <<cdiag<< " iwhich=" << iwhich_c0 << " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
    }else if (xir <= 1-P_acc[3] and xir > P_acc[0]+P_acc[1]){
      // accept c2
      iacc = 2;
      cpseudo = pseud_c[2].copy();
      cDlt_type = dlttyp_c[2].copy();
      cdiag = diag_c[2].copy();
      iwhich_c0 = iwhich_t2;
      pseudos_c0 = pseudos_trial; pseudos_trial.clear();
      if (Delta_types_c0.extent(0)!=Delta_types_trial.extent(0) or
	  Delta_types_c0.extent(1)!=Delta_types_trial.extent(1) ) Delta_types_c0.resize(Delta_types_trial.shape());
      Delta_types_c0 = Delta_types_trial;
      Delta_types_trial.resize(0,0);
      //Delta_types_c0 = Delta_types_trial; Delta_types_trial.clear();
      csgn_dlt = sgn_dlt;
      cout << "Choosing c2 icase=" << icase << " p=" << cpseudo << " diag=" <<cdiag<< " iwhich=" << iwhich_c0 << " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
    }else if (xir <= P_acc[0]+P_acc[1] and xir > P_acc[0]){
      iacc = 1;
      cpseudo = pseud_c[1].copy();
      cDlt_type = dlttyp_c[1].copy();
      iwhich_c0 = iwhich_t1;
      skl = skl_copy;
      cout<<"Choosing c1 icase="<<icase<<"  p="<<cpseudo<<" diag="<<cdiag;
      cout<<" iwhich="<< iwhich_c0<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
    } else if (xir <= P_acc[0]){ // c0 -> backtracking
      iacc = 0;
      cout << "************************* BACKTRACKING ***************" << endl;
      skl = skl_copy;
    }
    idiag += 1;

    auto t_6 = high_resolution_clock::now();
    t_next += t_6-t_5;
    
    if (idiag >= Ndiag) break;
  }
  //cout << "how_many_pgp=" << how_many_pgp/Ndiag/(2*Norder) << endl;
  
  //if (not Qstore){
  for (auto itt=stored2.begin(); itt!=stored2.end(); ++itt){
    auto [diag,pseud] = String2diagp(itt->first);
    stored.push_back(make_pair(diag,pseud));
    stored_count.push_back( itt->second );
    cout<<itt->first<<" "<<itt->second<<" diag="<<diag<<" pseud="<<pseud<<endl;
  }
  //}
  
  vector<int> sfc(stored.size());
  for (int i=0; i<stored.size(); ++i)
    sfc[i] = SymmetryFactor2(stored[i]);
  
  double rr = 0;
  for (int i=0; i<stored.size(); ++i){
    rr += sfc[i]*stored_count[i];
  }
  rr /= stored.size();
  
  double deviation=0;
  double Vol_diagrams=0;
  cout <<"stored diags="<<endl;
  for (int i=0; i<stored.size(); ++i){
    double x = stored_count[i]*sfc[i]/rr;
    deviation += abs(x-1);

    cdiag = stored[i].first;
    cpseudo = stored[i].second;

    int sgn_dlt = SignSort(cdiag);
    TinyVector<int,2> sgn_mmbw = SignMatrixElements(cdiag, cpseudo, Dlt);
    int sgn_tot = sgn_dlt*sgn_mmbw[0];
    int sgn_cur = sgn_dlt*sgn_mmbw[0]*sgn_mmbw[1];

    //double ssign = stored_sign[i]/stored_count[i];
    //if (fabs(ssign-sgn_tot)>1e-10){
    //  cout << "ERROR ssign=" << ssign << " while sgn_tot=" << sgn_tot << " should be the same but are not" << endl;
    //}
    cout<<setw(2)<< i << " " << setw(8) << setprecision(3) << stored_count[i] << " "<< setw(2)<< setprecision(5) << sfc[i] << "  "<< left << setw(8) << setprecision(6) << x << right << " " <<stored[i].first<<" p=" << stored[i].second;
    cout<<endl;
    //cout<<" ssgn=" << left << setw(3)<< ssign << " sgn_t=" << setw(2) << sgn_tot << " sgn_c=" << setw(2) << sgn_cur << " sgn_dlt=" << setw(2) << sgn_dlt << " sgn_mm=" << setw(2) << sgn_mmbw[0] << " sgn_bw=" << setw(2) << sgn_mmbw[1] << right << endl;

    Vol_diagrams += 1.0/sfc[i];
  }
  cout << "deviation=" << deviation/stored.size() << endl;
  
  cout << "Before: Pbin="<<Pbin*(V0norm/Pnorm)<<" dSm="<< Array<double,1>((dSm*(V0norm/Pnorm))) <<endl;
  Pbin *= (V0norm/Pnorm)*Vol_diagrams;
  dSm  *= (V0norm/Pnorm)*Vol_diagrams;

  cout << "After Pbin="<< Pbin <<" dSm="<<dSm<<endl;
  cout<<"t_next="<<t_next.count()<<" t_diff="<<t_diff.count()<<" t_choose="<<t_choose.count()<<" t_topol="<<t_topol.count()<<" t_skel="<<t_skel.count()<<endl;
}

template<class RNG>
void GenerateDiagrams0(Array<int,2>& i_diag,const vector<TinyVector<int,4>> non_zero,const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip, RNG& drand,int Ndiag=10000){
  //#define _DEBUG  
  cout.precision(12);
  int Norder=i_diag.extent(0);
  Array<int,2> diag(Norder,2), cdiag(Norder,2), ndiag(Norder,2);
  diag = i_diag;
  double P_first_jump = 0.5;

  Skeletons skl(diag);
  bool Reducible = skl.IsReducible();
  {
    int sgn_dlt=1;
    FindNextSkeleton(sgn_dlt, diag, skl,drand);
  }

  const int Ndd=4;
  TinyVector<Array<int,2>,Ndd> diag_c;
  TinyVector<Array<int,1>,Ndd> pseud_c;
  TinyVector<Array<int,1>,Ndd> dlttyp_c;
  TinyVector<int,Ndd> sgn_t;
  for (int i=0; i<Ndd; ++i){
    diag_c[i].resize(Norder,2);
    pseud_c[i].resize(2*Norder);
    dlttyp_c[i].resize(Norder);
  }
  
  deque<Array<int,1>> pseudos_c0, pseudos_trial;
  Array<int,2> Delta_types_c0, Delta_types_trial;
  Array<int,1> cpseudo(2*Norder), dpseudo(2*Norder), cDlt_type(Norder);
  
  Skeletons skl_copy(diag);
  TinyVector<int,2> p_known(0,0);  //@@
  Determine_diagrams_of_this_topology_with(p_known, pseudos_c0, Delta_types_c0, diag, non_zero, Dlt, Fi_Fip); //@@
  { //@@
    //Determine_all_diagrams_of_this_topology2(pseudos_c0, Delta_types_c0, diag, non_zero, Dlt, Fi_Fip);
    int iwhich_c0 = static_cast<int>(drand()*pseudos_c0.size());
    cpseudo = pseudos_c0[iwhich_c0].copy(); 
    cDlt_type = Delta_types_c0(iwhich_c0,Range::all());
  } //@@
  
  Array<int,1> how_many_equal(pseudos_c0.size()); how_many_equal = 0;
  for (int ip=0; ip<pseudos_c0.size(); ++ip){
    dpseudo = pseudos_c0[ip]-cpseudo;
    how_many_equal(ip) = count_zeros( dpseudo );
    cout<<ip<<" p="<< pseudos_c0[ip] << " has " << how_many_equal(ip) << " equal  cDlt=" << cDlt_type << endl;
  }
  cout << "Starting with diag=" << diag << " and p=" << cpseudo <<" cDlt=" << cDlt_type << endl;

  cdiag = diag;
  int csgn_dlt = SignSort(cdiag);

  int idiag=0;
  double Pnorm=0, Pbin=0, Pnorm_sum=0, Pbin_sum=0;
  int Nmeassure=0;
  //bool Qstore=false;
  map<string,double> stored2;
  deque<pair<Array<int,2>,Array<int,1>>> stored;
  deque<double> stored_count;
  deque<double> stored_sign;
  
  Array<double,1> dSm(Ndd); dSm=0;
  double V0norm=1.0;
  TinyVector<double,Ndd> P_acc, P_pr, P_eq;
  TinyVector<double,2> P_tree(0,V0norm);
  int iacc=0;
  int iwhich_t1, iwhich_t2, iwhich_t3;
  fsec t_next=zero_t, t_diff=zero_t, t_choose=zero_t, t_topol=zero_t, t_skel=zero_t;
  double how_many_pgp=0;

  for (int ii=0; ii<10000000; ++ii){
    auto t_1 = high_resolution_clock::now();
    diag_c[0] = cdiag;
    pseud_c[0] = cpseudo;
    dlttyp_c[0] = cDlt_type;
    sgn_t[0] = csgn_dlt;
    //
    skl_copy = skl;
    diag = cdiag;
    int sgn_dlt = csgn_dlt;
    FindNextSkeleton(sgn_dlt, diag, skl,drand);
    // All diagrams with all possible pseudoparticle arrangements are produced by recursion
    auto t_2 = high_resolution_clock::now();
    t_skel += t_2-t_1; 
    //Determine_all_diagrams_of_this_topology2(pseudos_trial, Delta_types_trial, diag, non_zero, Dlt, Fi_Fip);

    int i_known = drand()*2*Norder;  //@@
    TinyVector<int,2> p_known(i_known,cpseudo(i_known)); //@@
    cout<<"p_known="<<p_known[1]<<" i_known="<<p_known[0]<<" cpseudo="<<cpseudo<<endl;   //@@
    Determine_diagrams_of_this_topology_with(p_known, pseudos_c0, Delta_types_c0, cdiag, non_zero, Dlt, Fi_Fip);  //@@
#ifdef _DEBUG      
    for (int ip=0; ip<pseudos_c0.size(); ++ip) cout<<ip<<" p="<< pseudos_c0[ip] << endl; //@@
#endif      
    Determine_diagrams_of_this_topology_with(p_known, pseudos_trial, Delta_types_trial, diag, non_zero, Dlt, Fi_Fip); //@@
    cout<<"Number of diagrams in the old and new topology="<<pseudos_c0.size()<<"/"<<pseudos_trial.size()<<endl; //@@
    
    auto t_3 = high_resolution_clock::now();
    t_topol += t_3-t_2;
    cout<<ii<<" Number of diagrams in the new topology="<<pseudos_trial.size()<<endl;

    diag_c[1] = cdiag;
    int icase = 1;//(drand() < P_first_jump) ? 1 : 2;  //@@
#ifdef _DEBUG    
    cout << "icase=" << icase << " cdiag="<<cdiag<<" diag="<<diag<<" p_c0=" << cpseudo << " csgn_dlt=" << csgn_dlt << " sgn_dlt=" << sgn_dlt << endl;
#endif    
    how_many_equal.resize(pseudos_c0.size());
#ifdef _DEBUG    
    cout << "pseudos_c0.size=" << pseudos_c0.size() << endl;
#endif    
    int ic_to_compare = 0;
    for (int ip=0; ip<pseudos_c0.size(); ++ip){
      dpseudo = pseudos_c0[ip]-pseud_c[ic_to_compare];
      how_many_equal(ip) = count_zeros( dpseudo );
#ifdef _DEBUG      
      cout<<ip<<" p="<< pseudos_c0[ip] << " has" << how_many_equal(ip) << " equal" << endl;
#endif      
    }
    deque<int> candidates1 = which_have_same_equals_but_not_same(how_many_equal,2*Norder);  //@@
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates1 << endl;
#endif    
    iwhich_t1 = candidates1[ static_cast<int>(candidates1.size()*drand()) ];
    pseud_c[1] = pseudos_c0[iwhich_t1];
    dlttyp_c[1] = Delta_types_c0(iwhich_t1,Range::all());
    sgn_t[1] = csgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t1=" << iwhich_t1 << " p_c1=" << pseud_c[1] << endl;
#endif	  
    diag_c[2] = diag;
    how_many_equal.resize(pseudos_trial.size());
    ic_to_compare = 1;//(icase==1) ? 1 : 0;
    for (int ip=0; ip<pseudos_trial.size(); ++ip){
      dpseudo = pseudos_trial[ip]- pseud_c[ic_to_compare];
      how_many_equal(ip) = count_zeros( dpseudo );
#ifdef _DEBUG      
      cout<<ip<<" p="<< pseudos_trial[ip] << " has" << how_many_equal(ip) << " equal" << endl;
#endif      
    }
    deque<int> candidates2 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates2 << endl;
#endif     
    iwhich_t2 = candidates2[ static_cast<int>(candidates2.size()*drand()) ];
    pseud_c[2] = pseudos_trial[iwhich_t2];
    dlttyp_c[2] = Delta_types_trial(iwhich_t2,Range::all());
    sgn_t[2] = sgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t2=" << iwhich_t2 << " p_c2=" << pseud_c[2] << endl;
#endif
    diag_c[3] = diag;
    ic_to_compare = 2;
    for (int ip=0; ip<pseudos_trial.size(); ++ip){
      if (ip==iwhich_t2)
	how_many_equal(ip) = 0; // we should avoid the current diagram, because we want to change
      else{
	dpseudo = pseudos_trial[ip]-pseud_c[ic_to_compare];
	how_many_equal(ip) = count_zeros( dpseudo );
      }
#ifdef _DEBUG      
      cout<<ip<<" p="<< pseudos_trial[ip] << " has" << how_many_equal(ip) << " equal" << endl;
#endif      
    }
    deque<int> candidates3 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates3 << endl;
#endif    
    iwhich_t3 = candidates3[ static_cast<int>(candidates3.size()*drand()) ];
    pseud_c[3] = pseudos_trial[iwhich_t3];
    dlttyp_c[3] = Delta_types_trial(iwhich_t3,Range::all());
    sgn_t[3] = sgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t3=" << iwhich_t3 << " p_c3=" << pseud_c[3] << endl;
#endif	  
    auto t_4 = high_resolution_clock::now();
    t_choose += t_4-t_3;
    
    P_eq = 1.,1.,1.,1.;
    
    double N0 = pseudos_c0.size(), N1 = pseudos_trial.size();
    P_pr[0] = 0.5*N0/(N0+N1);
    P_pr[1] = 0.5*N0/(N0+N1);
    P_pr[2] = 0.5*N1/(N0+N1);
    P_pr[3] = 0.5*N1/(N0+N1);
    cout<<"P_pr="<<P_pr<<endl;
    P_tree[0]=0;
    for (int ic=0; ic<Ndd; ++ic){
      //double P0 = P_eq[ic]*P_pr[ic]; 
      P_tree[0] += P_eq[ic]*P_pr[ic];
      P_acc[ic] = (P_eq[ic] + P_tree[1])*P_pr[ic];
    }
    P_acc *= 1/(P_tree[0]+P_tree[1]);
#ifdef _DEBUG    
    cout << "P_tree="<<P_tree<<" P_acc=" << P_acc <<" P_pr="<<P_pr<<endl;
#endif
    for (int ic=0; ic<Ndd; ++ic)
      cout<<ii<<") ic="<<ic<<" d="<<diag_c[ic]<<" p="<<pseud_c[ic]<<" Pacc="<<P_acc[ic]<<endl;
    for (int ic=0; ic<Ndd; ++ic){
      TinyVector<int,2> sgn_mmbw = SignMatrixElements(diag_c[ic], pseud_c[ic], Dlt);
      sgn_t[ic] *= sgn_mmbw[0];
      auto [new_diag,s] = Different3(diag_c[ic], pseud_c[ic], P_acc[ic], stored2);
      if (new_diag)
	cout << "Final-Diag "<< diag_c[ic] << " with p=" << pseud_c[ic] << " is new diagram " << s << " wgh=" << P_acc[ic] << " sgn=" << sgn_t[ic] << endl;
      else
	cout << "Final-Diag "<< diag_c[ic] << " with p=" << pseud_c[ic] << " is equivalent to " << s << " wgh=" << P_acc[ic] << " sgn=" << sgn_t[ic] << endl;
    }
    auto t_5 = high_resolution_clock::now();
    t_diff += t_5-t_4;
    // meassure
    {
      Nmeassure += 1;
      double Wphys  = P_tree[0]/(P_tree[0]+P_tree[1]);
      double Walter = P_tree[1]/(P_tree[0]+P_tree[1]);
      Pnorm += Walter;
      Pnorm_sum += Walter;
      Pbin_sum  += Wphys;
      Pbin += Wphys;
      for (int ic=0; ic<Ndd; ++ic){
	double partialPhi = P_eq[ic]*P_pr[ic]/(P_tree[0]+P_tree[1]);
	for (int it=0; it<2*Norder; ++it){
	  int atomic_state = pseud_c[ic](it);
	  dSm(atomic_state) += partialPhi;
	}
      }
    }
    double xir = drand();
    cout<<"xir="<<xir<<" P_acc="<<P_acc<<endl;
    if (xir > 1-P_acc[3] ){ //@@
      // accept c3
      iacc = 3;
      cpseudo = pseud_c[3].copy();
      cDlt_type = dlttyp_c[3].copy();
      cdiag = diag_c[3].copy();
      //iwhich_c0 = iwhich_t3; //@@
      //pseudos_c0 = pseudos_trial; pseudos_trial.clear(); //@@
      //Delta_types_c0 = Delta_types_trial; Delta_types_trial.clear(); //@@
      csgn_dlt = sgn_dlt;
      cout << "Choosing c3 icase=" << icase << " p=" << cpseudo << " diag=" <<cdiag<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
    }else if (xir <= 1-P_acc[3] and xir > P_acc[0]+P_acc[1]){
      // accept c2
      iacc = 2;
      cpseudo = pseud_c[2].copy();
      cDlt_type = dlttyp_c[2].copy();
      cdiag = diag_c[2].copy();
      //iwhich_c0 = iwhich_t2; //@@
      //pseudos_c0 = pseudos_trial; pseudos_trial.clear(); //@@
      //Delta_types_c0 = Delta_types_trial; Delta_types_trial.clear(); //@@
      csgn_dlt = sgn_dlt;
      cout << "Choosing c2 icase=" << icase << " p=" << cpseudo << " diag=" <<cdiag<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
    }else if (xir <= P_acc[0]+P_acc[1] and xir > P_acc[0]){
      iacc = 1;
      cpseudo = pseud_c[1].copy();
      cDlt_type = dlttyp_c[1].copy();
      //iwhich_c0 = iwhich_t1;  //@@
      skl = skl_copy;
      cout<<"Choosing c1 icase="<<icase<<"  p="<<cpseudo<<" diag="<<cdiag<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
    } else if (xir <= P_acc[0]){ // c0 -> backtracking
      iacc = 0;
      cout << "************************* BACKTRACKING ***************" << endl;
      skl = skl_copy;
    }
    idiag += 1;

    auto t_6 = high_resolution_clock::now();
    t_next += t_6-t_5;
    
    if (idiag >= Ndiag) break;
  }
  //cout << "how_many_pgp=" << how_many_pgp/Ndiag/(2*Norder) << endl;
  
  for (auto itt=stored2.begin(); itt!=stored2.end(); ++itt){
    auto [diag,pseud] = String2diagp(itt->first);
    stored.push_back(make_pair(diag,pseud));
    stored_count.push_back( itt->second );
    cout<<itt->first<<" "<<itt->second<<" diag="<<diag<<" pseud="<<pseud<<endl;
  }
  
  vector<int> sfc(stored.size());
  for (int i=0; i<stored.size(); ++i)
    sfc[i] = SymmetryFactor2(stored[i]);
  
  double rr = 0;
  for (int i=0; i<stored.size(); ++i){
    rr += sfc[i]*stored_count[i];
  }
  rr /= stored.size();
  
  double deviation=0;
  double Vol_diagrams=0;
  cout <<"stored diags="<<endl;
  for (int i=0; i<stored.size(); ++i){
    double x = stored_count[i]*sfc[i]/rr;
    deviation += abs(x-1);

    cdiag = stored[i].first;
    cpseudo = stored[i].second;

    int sgn_dlt = SignSort(cdiag);
    TinyVector<int,2> sgn_mmbw = SignMatrixElements(cdiag, cpseudo, Dlt);
    int sgn_tot = sgn_dlt*sgn_mmbw[0];
    int sgn_cur = sgn_dlt*sgn_mmbw[0]*sgn_mmbw[1];

    cout<<setw(2)<< i << " " << setw(8) << setprecision(3) << stored_count[i] << " "<< setw(2)<< setprecision(5) << sfc[i] << "  "<< left << setw(8) << setprecision(6) << x << right << " " <<stored[i].first<<" p=" << stored[i].second;
    cout<<endl;

    Vol_diagrams += 1.0/sfc[i];
  }
  cout << "deviation=" << deviation/stored.size() << endl;
  
  cout << "Before: Pbin="<<Pbin*(V0norm/Pnorm)<<" dSm="<< Array<double,1>((dSm*(V0norm/Pnorm))) <<endl;
  Pbin *= (V0norm/Pnorm)*Vol_diagrams;
  dSm  *= (V0norm/Pnorm)*Vol_diagrams;

  cout << "After Pbin="<< Pbin <<" dSm="<<dSm<<endl;
  cout<<"t_next="<<t_next.count()<<" t_diff="<<t_diff.count()<<" t_choose="<<t_choose.count()<<" t_topol="<<t_topol.count()<<" t_skel="<<t_skel.count()<<endl;
}

inline int Which_Case(double xi, const TinyVector<double,4>& P_acc)
{
  int icc=0;
  double sum_Pacc = P_acc[icc];
  //cout<<"P_acc="<<P_acc<<endl;
  //cout<<"icc="<<icc<<" sum_Pacc="<<sum_Pacc<<" xi="<<xi<<endl;
  while (icc<4 and sum_Pacc<xi){
    icc += 1;
    sum_Pacc += P_acc[icc];
    //cout<<"icc="<<icc<<" sum_Pacc="<<sum_Pacc<<" xi="<<xi<<endl;
  }
  return icc;
}
template<class RNG>
void GenerateDiagrams0_(Array<int,2>& i_diag,const vector<TinyVector<int,4>> non_zero,const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip, RNG& drand,int Ndiag=10000){
  //#define _DEBUG  
  cout.precision(12);
  int Norder=i_diag.extent(0);
  Array<int,2> diag(Norder,2), cdiag(Norder,2), ndiag(Norder,2);
  diag = i_diag;
  double P_first_jump = 0.5;

  Skeletons skl(diag);
  bool Reducible = skl.IsReducible();
  {
    int sgn_dlt=1;
    FindNextSkeleton(sgn_dlt, diag, skl,drand);
  }

  int Ndd=4;
  TinyVector<Array<int,2>,4> diag_c;
  TinyVector<Array<int,1>,4> pseud_c;
  TinyVector<Array<int,1>,4> dlttyp_c;
  TinyVector<int,4> sgn_t;
  for (int i=0; i<4; ++i){
    diag_c[i].resize(Norder,2);
    pseud_c[i].resize(2*Norder);
    dlttyp_c[i].resize(Norder);
  }
  
  deque<Array<int,1>> pseudos_c0, pseudos_trial;
  Array<int,2> Delta_types_c0, Delta_types_trial;
  Array<int,1> cpseudo(2*Norder), dpseudo(2*Norder), cDlt_type(Norder);
  
  Skeletons skl_copy(diag);
  TinyVector<int,2> p_known(0,0);  //@@
  Determine_diagrams_of_this_topology_with(p_known, pseudos_c0, Delta_types_c0, diag, non_zero, Dlt, Fi_Fip); //@@
  { //@@
    int iwhich_c0 = static_cast<int>(drand()*pseudos_c0.size());
    cpseudo = pseudos_c0[iwhich_c0].copy(); 
    cDlt_type = Delta_types_c0(iwhich_c0,Range::all());
  } //@@
  
  Array<int,1> how_many_equal(pseudos_c0.size()); how_many_equal = 0;
  for (int ip=0; ip<pseudos_c0.size(); ++ip){
    dpseudo = pseudos_c0[ip]-cpseudo;
    how_many_equal(ip) = count_zeros( dpseudo );
    cout<<ip<<" p="<< pseudos_c0[ip] << " has " << how_many_equal(ip) << " equal  cDlt=" << cDlt_type << endl;
  }
  cout << "Starting with diag=" << diag << " and p=" << cpseudo <<" cDlt=" << cDlt_type << endl;

  cdiag = diag;
  int csgn_dlt = SignSort(cdiag);

  int idiag=0;
  double Pnorm=0, Pbin=0, Pnorm_sum=0, Pbin_sum=0;
  int Nmeassure=0;
  //bool Qstore=false;
  map<string,double> stored2;
  deque<pair<Array<int,2>,Array<int,1>>> stored;
  deque<double> stored_count;
  deque<double> stored_sign;
  
  Array<double,1> dSm(4); dSm=0;
  double V0norm=1.0;
  TinyVector<double,4> P_acc, P_pr, P_eq;
  TinyVector<double,2> P_tree(0,V0norm);
  int iacc=0;
  int iwhich_t1, iwhich_t2, iwhich_t3;
  fsec t_next=zero_t, t_diff=zero_t, t_choose=zero_t, t_topol=zero_t, t_skel=zero_t;
  double how_many_pgp=0;

  for (int ii=0; ii<10000000; ++ii){
    auto t_1 = high_resolution_clock::now();
    diag_c[0] = cdiag;
    pseud_c[0] = cpseudo;
    dlttyp_c[0] = cDlt_type;
    sgn_t[0] = csgn_dlt;
    //
    skl_copy = skl;
    diag = cdiag;
    int sgn_dlt = csgn_dlt;
    if (not skl.AreLoopsCompatibleWithDiag(diag)){
      cout<<ii<<") ERROR diag&loops incompatible: diag="<<diag<<" and loops="<<skl.loops<<endl;
      exit(1);
    }
    FindNextSkeleton(sgn_dlt, diag, skl,drand);
    // All diagrams with all possible pseudoparticle arrangements are produced by recursion
    auto t_2 = high_resolution_clock::now();
    t_skel += t_2-t_1; 
    //Determine_all_diagrams_of_this_topology2(pseudos_trial, Delta_types_trial, diag, non_zero, Dlt, Fi_Fip);

    int i_known = drand()*2*Norder;  //@@
    //TinyVector<int,2> p_known(i_known,cpseudo(i_known)); //@@
    int n_cmn=max(1,Norder/2);
    deque<int> p_known(n_cmn);
    for (int i=0; i<n_cmn; ++i){p_known[i]= cpseudo((i_known+i)%(2*Norder));}
    cout<<"p_known="<<p_known<<" i_known="<<i_known<<" n_cmn="<<n_cmn<<" cpseudo="<<cpseudo<<endl;   //@@
    Determine_diagrams_of_this_topology_with(p_known, i_known, pseudos_c0, Delta_types_c0, cdiag, non_zero, Dlt, Fi_Fip);  //@@
    Determine_diagrams_of_this_topology_with(p_known, i_known, pseudos_trial, Delta_types_trial, diag, non_zero, Dlt, Fi_Fip); //@@
    cout<<"Number of diagrams in the old and new topology="<<pseudos_c0.size()<<"/"<<pseudos_trial.size()<<endl; //@@
    if (pseudos_c0.size()<1 or pseudos_trial.size()<1){
      n_cmn-=1;
      p_known.resize(n_cmn);
      for (int i=0; i<n_cmn; ++i){p_known[i]= cpseudo((i_known+i)%(2*Norder));}
      cout<<"p_known="<<p_known<<" i_known="<<i_known<<" n_cmn="<<n_cmn<<" cpseudo="<<cpseudo<<endl;   //@@
      Determine_diagrams_of_this_topology_with(p_known, i_known, pseudos_c0, Delta_types_c0, cdiag, non_zero, Dlt, Fi_Fip);  //@@
      Determine_diagrams_of_this_topology_with(p_known, i_known, pseudos_trial, Delta_types_trial, diag, non_zero, Dlt, Fi_Fip); //@@
      cout<<"Number of diagrams in the old and new topology="<<pseudos_c0.size()<<"/"<<pseudos_trial.size()<<endl; //@@
    }
#ifdef _DEBUG      
    for (int ip=0; ip<pseudos_c0.size(); ++ip) cout<<ip<<" p="<< pseudos_c0[ip] << endl; //@@
#endif      
    
    if (pseudos_c0.size()<1 or pseudos_trial.size()<1){
      cout<<"ERROR size of pseudo_c0="<<pseudos_c0.size()<<" and pseudo_trial="<<pseudos_trial.size()<<endl;
      Determine_all_diagrams_of_this_topology2(pseudos_c0, Delta_types_c0, cdiag, non_zero, Dlt, Fi_Fip);
      cout<<"cdiag="<<cdiag<<endl;
      for (int ip=0; ip<pseudos_c0.size(); ++ip)
	cout<<ip<<" p="<< pseudos_c0[ip] << endl;
      cout<<"diag="<<diag<<endl;
      Determine_all_diagrams_of_this_topology2(pseudos_trial, Delta_types_trial, diag, non_zero, Dlt, Fi_Fip);
      for (int ip=0; ip<pseudos_trial.size(); ++ip)
	cout<<ip<<" p="<< pseudos_trial[ip] << endl;
      exit(0);
    }
    Ndd=4;
    if (pseudos_c0.size()==1) Ndd-=1; // tree is smaller now with only 3 entries.
    if (pseudos_trial.size()==1) Ndd-=1; // tree is smaller now with only 3 entries.
    
    auto t_3 = high_resolution_clock::now();
    t_topol += t_3-t_2;
    cout<<ii<<" Number of diagrams in the new topology="<<pseudos_trial.size()<<" Ndd="<<Ndd<<endl;

    //cout<<"4: skl.loops[2].size()="<<skl.loops[2].size()<<endl;
    
    int idc=1;
    if (pseudos_c0.size()>1){
      diag_c[idc] = cdiag;
      //int icase = 1;//(drand() < P_first_jump) ? 1 : 2;  //@@
#ifdef _DEBUG    
      cout << " cdiag="<<cdiag<<" diag="<<diag<<" p_c0=" << cpseudo << " csgn_dlt=" << csgn_dlt << " sgn_dlt=" << sgn_dlt << endl;
#endif    
      how_many_equal.resize(pseudos_c0.size());
#ifdef _DEBUG    
      cout << "pseudos_c0.size=" << pseudos_c0.size() << endl;
#endif    
      int ic_to_compare = 0;
      for (int ip=0; ip<pseudos_c0.size(); ++ip){
	dpseudo = pseudos_c0[ip]-pseud_c[ic_to_compare];
	how_many_equal(ip) = count_zeros( dpseudo );
#ifdef _DEBUG      
	cout<<ip<<" p="<< pseudos_c0[ip] << " has" << how_many_equal(ip) << " equal" << endl;
#endif      
      }
      deque<int> candidates1 = which_have_same_equals_but_not_same(how_many_equal,2*Norder);  //@@
#ifdef _DEBUG    
      cout << "candidates=" <<  candidates1 << endl;
#endif
      iwhich_t1 = candidates1[ static_cast<int>(candidates1.size()*drand()) ];
      pseud_c[idc] = pseudos_c0[iwhich_t1];
      dlttyp_c[idc] = Delta_types_c0(iwhich_t1,Range::all());
      sgn_t[idc] = csgn_dlt;
#ifdef _DEBUG    
      cout << "iwhich_t1=" << iwhich_t1 << " p_c1=" << pseud_c[idc] << endl;
#endif
      idc+=1;
    }
    diag_c[idc] = diag;
    how_many_equal.resize(pseudos_trial.size());
    int ic_to_compare = idc-1;
    for (int ip=0; ip<pseudos_trial.size(); ++ip){
      dpseudo = pseudos_trial[ip]- pseud_c[ic_to_compare];
      how_many_equal(ip) = count_zeros( dpseudo );
#ifdef _DEBUG      
      cout<<ip<<" p="<< pseudos_trial[ip] << " has" << how_many_equal(ip) << " equal" << endl;
#endif      
    }
    deque<int> candidates2 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates2 << endl;
#endif     
    iwhich_t2 = candidates2[ static_cast<int>(candidates2.size()*drand()) ];
    pseud_c[idc] = pseudos_trial[iwhich_t2];
    dlttyp_c[idc] = Delta_types_trial(iwhich_t2,Range::all());
    sgn_t[idc] = sgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t2=" << iwhich_t2 << " p_c2=" << pseud_c[idc] << endl;
#endif
    idc+=1;
    if (pseudos_trial.size()>1){
      diag_c[idc] = diag;
      ic_to_compare = idc-1;
      for (int ip=0; ip<pseudos_trial.size(); ++ip){
	if (ip==iwhich_t2)
	  how_many_equal(ip) = 0; // we should avoid the current diagram, because we want to change
	else{
	  dpseudo = pseudos_trial[ip]-pseud_c[ic_to_compare];
	  how_many_equal(ip) = count_zeros( dpseudo );
	}
#ifdef _DEBUG      
	cout<<ip<<" p="<< pseudos_trial[ip] << " has" << how_many_equal(ip) << " equal" << endl;
#endif
      }
      deque<int> candidates3 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
      cout << "candidates=" <<  candidates3 << endl;
#endif    
      iwhich_t3 = candidates3[ static_cast<int>(candidates3.size()*drand()) ];
      pseud_c[idc] = pseudos_trial[iwhich_t3];
      dlttyp_c[idc] = Delta_types_trial(iwhich_t3,Range::all());
      sgn_t[idc] = sgn_dlt;
#ifdef _DEBUG    
      cout << "iwhich_t3=" << iwhich_t3 << " p_c3=" << pseud_c[3] << endl;
#endif
      idc+=1;
    }
    if (Ndd!=idc) cout<<"ERROR we expected Ndd="<<Ndd<<" to be equal to idc="<<idc<<endl;
    
    auto t_4 = high_resolution_clock::now();
    t_choose += t_4-t_3;

    P_eq = 1.,1.,1.,1.;

    double N0 = pseudos_c0.size(), N1 = pseudos_trial.size();
    P_pr=0,0,0,0;
    idc=0;
    if (pseudos_c0.size()>1){
      P_pr[idc  ]=0.5*N0;
      P_pr[idc+1]=0.5*N0;
      idc+=2;
    }else{
      P_pr[idc]=N0;
      idc+=1;
    }
    if (pseudos_trial.size()>1){
      P_pr[idc  ]=0.5*N1;
      P_pr[idc+1]=0.5*N1;
      idc+=2;
    }else{
      P_pr[idc]=N1;
      idc+=1;
    }
    /*
    P_pr[idc++]=N0;
    if (pseudos_c0.size()>1) P_pr[idc++]=N0;
    P_pr[idc++]=N1;
    if (pseudos_trial.size()>1) P_pr[idc++]=N1;
    */
    if (Ndd!=idc) cout<<"ERROR we expected Ndd="<<Ndd<<" to be equal to idc="<<idc<<endl;
    //cout<<"before P_pr="<<P_pr<<endl;
    double dsum=0;
    for (int i=0; i<idc; ++i) dsum += P_pr[i];
    P_pr *= 1.0/dsum;
    cout<<"P_pr="<<P_pr<<endl;
#ifdef _DEBUG
    cout<<"P_pr="<<P_pr<<" Ndd="<<Ndd<<endl;
#endif    
    //P_pr[0] = 0.5*N0/(N0+N1);
    //P_pr[1] = 0.5*N0/(N0+N1);
    //P_pr[2] = 0.5*N1/(N0+N1);
    //P_pr[3] = 0.5*N1/(N0+N1);
    
    P_tree=0,1;
    P_acc=0,0,0,0;
    for (int ic=0; ic<Ndd; ++ic){
      //double P0 = P_eq[ic]*P_pr[ic]; 
      P_tree[0] += P_eq[ic]*P_pr[ic];
      P_acc[ic] = (P_eq[ic] + P_tree[1])*P_pr[ic];
    }
    P_acc *= 1/(P_tree[0]+P_tree[1]);
    //P_tree[0] *= (Ndd/4.0);
#ifdef _DEBUG    
    cout << "P_tree="<<P_tree<<" P_acc=" << P_acc <<" P_pr="<<P_pr<<endl;
#endif
    for (int ic=0; ic<Ndd; ++ic)
      cout<<ii<<") ic="<<ic<<" d="<<diag_c[ic]<<" p="<<pseud_c[ic]<<" Pacc="<<P_acc[ic]<<endl;
    for (int ic=0; ic<Ndd; ++ic){
      TinyVector<int,2> sgn_mmbw = SignMatrixElements(diag_c[ic], pseud_c[ic], Dlt);
      sgn_t[ic] *= sgn_mmbw[0];
      auto [new_diag,s] = Different3(diag_c[ic], pseud_c[ic], P_acc[ic], stored2);
      if (new_diag)
	cout << "Final-Diag "<< diag_c[ic] << " with p=" << pseud_c[ic] << " is new diagram " << s << " wgh=" << P_acc[ic] << " sgn=" << sgn_t[ic] << endl;
      else
	cout << "Final-Diag "<< diag_c[ic] << " with p=" << pseud_c[ic] << " is equivalent to " << s << " wgh=" << P_acc[ic] << " sgn=" << sgn_t[ic] << endl;
    }
    auto t_5 = high_resolution_clock::now();
    t_diff += t_5-t_4;

    // meassure
    {
      Nmeassure += 1;
      double Wphys  = P_tree[0]/(P_tree[0]+P_tree[1]);
      double Walter = P_tree[1]/(P_tree[0]+P_tree[1]);
      Pnorm += Walter;
      Pnorm_sum += Walter;
      Pbin_sum  += Wphys;
      Pbin += Wphys;
      for (int ic=0; ic<Ndd; ++ic){
	double partialPhi = P_eq[ic]*P_pr[ic]/(P_tree[0]+P_tree[1]);
	for (int it=0; it<2*Norder; ++it){
	  int atomic_state = pseud_c[ic](it);
	  if (atomic_state>4 or atomic_state<0){ cout<<"ERROR atomic_state="<<atomic_state<<endl;exit(1);}
	  dSm(atomic_state) += partialPhi;
	}
      }
    }
    double xir = drand();
    int iacc = Which_Case(xir, P_acc);
    if (iacc>=Ndd){
      cout<<"ERROR It should not happen! iacc="<< iacc<<" P_acc="<<P_acc<<" xi="<<xir<<endl;
      exit(0);
    }
    cout <<"xir="<<xir<<" P_acc="<<P_acc<<" iacc="<<iacc<<endl;
    switch(iacc){
    case 3:
      //if (xir > 1-P_acc[3] ){ //@@
      // accept c3
      // occurs only if pseudos_c0.size()==2 and pseudos_trial.size()==2
      cpseudo = pseud_c[iacc].copy();
      cDlt_type = dlttyp_c[iacc].copy();
      cdiag = diag_c[iacc].copy();
      csgn_dlt = sgn_dlt;
      cout << "Choosing c3,  p=" << cpseudo << " diag=" <<cdiag<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      break;
      //}else if (xir <= 1-P_acc[3] and xir > P_acc[0]+P_acc[1]){
    case 2:
      // accept c2
      // occurs only if pseudos_c0.size() + pseudos_trial.size() >=3. It can be (2,2) or (2,1) or (1,2).
      // In all these cases, we are accepting the new diagram, not the old diagram.
      cpseudo = pseud_c[iacc].copy();
      cDlt_type = dlttyp_c[iacc].copy();
      cdiag = diag_c[iacc].copy();
      csgn_dlt = sgn_dlt;
      cout << "Choosing c2,  p=" << cpseudo << " diag=" <<cdiag<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      break;
      //}else if (xir <= P_acc[0]+P_acc[1] and xir > P_acc[0]){
    case 1:
      // It can occur in any case. But if (pseudos_c0.size(),pseudos_trial.size()) is either (1,1) or (1,2) we accept previous diag configuration, and
      // when it is (2,1) or (2,2) we accept old configuration, i.e., backtracking
      cpseudo = pseud_c[iacc].copy();
      cDlt_type = dlttyp_c[iacc].copy();
      if (pseudos_c0.size()==1){
	cdiag = diag_c[iacc].copy();
	csgn_dlt = sgn_dlt;
      }else{
	skl = skl_copy;
      }
      cout<<"Choosing c1,  p="<<cpseudo<<" diag="<<cdiag<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      break;
      //} else if (xir <= P_acc[0]){ // c0 -> backtracking
    case 0:
      cout << "************************* BACKTRACKING ***************" << endl;
      skl = skl_copy;
      break;
    default:
      cout<<"ERROR It should not happen"<<endl; exit(0);
    }
    idiag += 1;

    auto t_6 = high_resolution_clock::now();
    t_next += t_6-t_5;
    
    if (idiag >= Ndiag) break;
  }
  //cout << "how_many_pgp=" << how_many_pgp/Ndiag/(2*Norder) << endl;
  
  for (auto itt=stored2.begin(); itt!=stored2.end(); ++itt){
    auto [diag,pseud] = String2diagp(itt->first);
    stored.push_back(make_pair(diag,pseud));
    stored_count.push_back( itt->second );
    cout<<itt->first<<" "<<itt->second<<" diag="<<diag<<" pseud="<<pseud<<endl;
  }
  
  vector<int> sfc(stored.size());
  for (int i=0; i<stored.size(); ++i)
    sfc[i] = SymmetryFactor2(stored[i]);
  
  double rr = 0;
  for (int i=0; i<stored.size(); ++i){
    rr += sfc[i]*stored_count[i];
  }
  rr /= stored.size();
  
  double deviation=0;
  double Vol_diagrams=0;
  cout <<"stored diags="<<endl;
  for (int i=0; i<stored.size(); ++i){
    double x = stored_count[i]*sfc[i]/rr;
    deviation += abs(x-1);

    cdiag = stored[i].first;
    cpseudo = stored[i].second;

    int sgn_dlt = SignSort(cdiag);
    TinyVector<int,2> sgn_mmbw = SignMatrixElements(cdiag, cpseudo, Dlt);
    int sgn_tot = sgn_dlt*sgn_mmbw[0];
    int sgn_cur = sgn_dlt*sgn_mmbw[0]*sgn_mmbw[1];

    cout<<setw(2)<< i << " " << setw(8) << setprecision(3) << stored_count[i] << " "<< setw(2)<< setprecision(5) << sfc[i] << "  "<< left << setw(8) << setprecision(6) << x << right << " " <<stored[i].first<<" p=" << stored[i].second;
    cout<<endl;

    Vol_diagrams += 1.0/sfc[i];
  }
  cout << "deviation=" << deviation/stored.size() << endl;
  
  cout << "Before: Pbin="<<Pbin*(V0norm/Pnorm)<<" dSm="<< Array<double,1>((dSm*(V0norm/Pnorm))) <<endl;
  Pbin *= (V0norm/Pnorm)*Vol_diagrams;
  dSm  *= (V0norm/Pnorm)*Vol_diagrams;

  cout << "After Pbin="<< Pbin <<" dSm="<<dSm<<endl;
  cout<<"t_next="<<t_next.count()<<" t_diff="<<t_diff.count()<<" t_choose="<<t_choose.count()<<" t_topol="<<t_topol.count()<<" t_skel="<<t_skel.count()<<endl;
}
template<class RNG>
void GenerateDiagrams20(Array<int,2>& i_diag,const vector<TinyVector<int,4>> non_zero,const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip, RNG& drand,int Ndiag=10000){
  //#define _DEBUG  
  cout.precision(12);
  int N_min_order=1, N_max_order=i_diag.extent(0);
  Array<double,1> Vol_diagrams0(N_max_order-N_min_order+1); Vol_diagrams0=1;
  Array<double,1> Vol_diagrams1(N_max_order-N_min_order+1); Vol_diagrams1=1;
  {
    Array<double,1> Vol_d0(8); Vol_d0 = 1,4,2,16./3.,27,198.4,1886.67,21846.607;
    Array<double,1> Vol_d1(8); Vol_d1 = 1,4,2,16./3.,36,310.4,3424.00,45427.857;
    TinyVector<double,3> pv(0.15760576,0.5415707,-1.3749251);
    for (int n=N_min_order; n<=N_max_order; ++n){
      Vol_diagrams0(n-N_min_order) = (n<8) ? Vol_d0(n) : pv[0]*n*n+pv[1]*n+pv[2];
      Vol_diagrams1(n-N_min_order) = (n<8) ? Vol_d1(n) : pv[0]*n*n+pv[1]*n+pv[2];
    }
  }
  //Array<double,1> Vol_diagrams0(max(N_max_order+1,8)); Vol_diagrams0=1;
  //Vol_diagrams0 = 1,4,2,16./3.,27,198.4,1886.67,21846.607;
  double P_first_jump = 0.5;

  deque<Array<int,1>> pseudos_c0, pseudos_trial;
  Array<int,2> Delta_types_c0, Delta_types_trial;
  
  int Norder=i_diag.extent(0);
  Array<int,2> diag(Norder,2), cdiag(Norder,2), ndiag(Norder,2);
  diag = i_diag;
  Skeletons skl(diag);
  bool Reducible = skl.IsReducible();
  {
    int sgn_dlt=1;
    FindNextSkeleton(sgn_dlt, diag, skl,drand);
  }
  Skeletons tskl(diag);
  
  TinyVector<Array<int,2>,4> diag_c;
  TinyVector<Array<int,1>,4> pseud_c;
  TinyVector<Array<int,1>,4> dlttyp_c;
  TinyVector<int,4> sgn_t;
  for (int i=0; i<4; ++i){
    diag_c[i].resize(Norder,2);
    pseud_c[i].resize(2*Norder);
    dlttyp_c[i].resize(Norder);
  }
  Array<int,1> cpseudo(2*Norder), dpseudo(2*Norder), dpseudo2(2*Norder), cDlt_type(Norder);

  TinyVector<int,2> p_known(0,0);  //@@
  Determine_diagrams_of_this_topology_with(p_known, pseudos_c0, Delta_types_c0, diag, non_zero, Dlt, Fi_Fip); //@@
  { //@@
    //Determine_all_diagrams_of_this_topology2(pseudos_c0, Delta_types_c0, diag, non_zero, Dlt, Fi_Fip);
    int iwhich_c0 = static_cast<int>(drand()*pseudos_c0.size());
    cpseudo = pseudos_c0[iwhich_c0].copy(); 
    cDlt_type = Delta_types_c0(iwhich_c0,Range::all());
  } //@@
  
  Array<int,1> how_many_equal(pseudos_c0.size()); how_many_equal = 0;
  for (int ip=0; ip<pseudos_c0.size(); ++ip){
    dpseudo = pseudos_c0[ip]-cpseudo;
    how_many_equal(ip) = count_zeros( dpseudo );
    cout<<ip<<" p="<< pseudos_c0[ip] << " has " << how_many_equal(ip) << " equal  cDlt=" << cDlt_type << endl;
  }
  cout << "Starting with diag=" <<diag<< " and p=" << cpseudo <<" cDlt=" << cDlt_type << endl;
  cdiag = diag;
  int csgn_dlt = SignSort(cdiag);

  int idiag=0;
  double Pnorm=0, Pbin=0, Pnorm_sum=0, Pbin_sum=0;
  Array<double,1> Nmeassure(N_max_order+1), asign(N_max_order+1);
  Nmeassure=0; asign=0;
  map<string,double> stored2;
  deque<pair<Array<int,2>,Array<int,1>>> stored;
  deque<double> stored_count;
  deque<double> stored_sign;
  Array<int,1> o2n;
  int tNorder=Norder;
  tskl = skl;
  Array<double,1> dSm(4); dSm=0;
  double V0norm=1.0;
  TinyVector<double,4> P_acc, P_pr, P_eq;
  double P_tree=0;
  int iacc=0;
  int iwhich_t1, iwhich_t2, iwhich_t3;
  fsec t_next=zero_t, t_diff=zero_t, t_choose=zero_t, t_topol=zero_t, t_skel=zero_t;
  double how_many_pgp=0;
  for (int ii=0; ii<100000000; ++ii){
    auto t_1 = high_resolution_clock::now();
    diag_c[0] = cdiag;
    pseud_c[0] = cpseudo;
    dlttyp_c[0] = cDlt_type;
    sgn_t[0] = csgn_dlt;
    //

    //int i_known = drand()*2*Norder;  //@@
    Norder = cdiag.extent(0);
    int sgn_dlt = csgn_dlt;
    int jump=0;
    double xi=drand();
    if (xi < 1./3.) jump = 1;
    else if ( xi < 2./3.) jump = -1;
    if (jump==1 and Norder >= N_max_order) jump=0;
    if (jump==-1 and Norder <= N_min_order) jump=0;
#ifdef _DEBUG    
    cout<<ii<<" diag="<<diag<<" cdiag="<<cdiag<<" jump="<<jump<<endl;
#endif    
    if (jump==1){
      tNorder = Norder+1;
      TinyVector<pair<double,int>,2>  hybn = GenerateRandomHybridization(Norder,TinyVector<double,2>(drand(),drand()));
      diag.resize(tNorder,2);
      diag = InsertHybridization(cdiag,hybn,o2n);
      tskl.Resize_Initialize(diag);
      if (tskl.IsReducible()) FindNextSkeleton(sgn_dlt, diag, tskl,drand);
      sgn_dlt = SignSort(diag);
      dpseudo2.resize(2*tNorder);
      diag_c[2].resize(tNorder,2);
      pseud_c[2].resize(2*tNorder);
      dlttyp_c[2].resize(tNorder);
      diag_c[3].resize(tNorder,2);
      pseud_c[3].resize(2*tNorder);
      dlttyp_c[3].resize(tNorder);
    }else if (jump==-1){
      int iremove = Norder*drand();
      tNorder = Norder-1;
      diag.resize(tNorder,2);
      diag = RemoveHybridization(cdiag, iremove,o2n);
      tskl.Resize_Initialize(diag);
      if (tskl.IsReducible()) FindNextSkeleton(sgn_dlt, diag, tskl,drand);
      sgn_dlt = SignSort(diag);
      dpseudo2.resize(2*tNorder);
      diag_c[2].resize(tNorder,2);
      pseud_c[2].resize(2*tNorder);
      dlttyp_c[2].resize(tNorder);
      diag_c[3].resize(tNorder,2);
      pseud_c[3].resize(2*tNorder);
      dlttyp_c[3].resize(tNorder);
    }else{
      tNorder = Norder;
      diag = cdiag;
      FindNextSkeleton(sgn_dlt, diag, tskl,drand);
    }
    //if (jump) cout<<"diag="<<diag<<" o2n="<<o2n<<endl;
    // All diagrams with all possible pseudoparticle arrangements are produced by recursion
    auto t_2 = high_resolution_clock::now();
    t_skel += t_2-t_1;
    /*
    int i_known=-1;
    if (jump){
      deque<int> possible_j;
      for (int j=0; j<2*Norder; ++j){
	int jt = o2n(j);
	int kt = o2n((j+1)%(2*Norder)); // two pseudo-particle propagators are changed when we remove a time at certain vertex, hence also kn should be positive
	bool nothing_inbetween = ((kt-jt==1) or (kt-jt==-2*tNorder+1)) and (jt>=0 and kt>=0);
	if (nothing_inbetween) possible_j.push_back(j);
      }
      i_known = possible_j[static_cast<int>(drand()*possible_j.size())];  //@@
    }else{
      i_known = drand()*2*Norder;
    }
    TinyVector<int,2> p_known(i_known,cpseudo(i_known)); //@@
    TinyVector<int,2> q_known; //@@
    if (jump!=0)
      q_known = o2n(i_known), cpseudo(i_known); //@@
    else
      q_known = i_known,cpseudo(i_known);
    */
    auto [p_known, q_known] = GiveCommonRandomPseudoPropagator(TinyVector<int,2>(Norder,tNorder),cpseudo,o2n,drand());
    cout<<"p_known="<<p_known<<" q_known="<<q_known<<endl;   //@@
    Determine_diagrams_of_this_topology_with(p_known, pseudos_c0, Delta_types_c0, cdiag, non_zero, Dlt, Fi_Fip);  //@@
    cout<<"Number of diagrams in the old topology="<<pseudos_c0.size()<<endl; //@@
#ifdef _DEBUG      
    for (int ip=0; ip<pseudos_c0.size(); ++ip) cout<<ip<<" p="<< pseudos_c0[ip] << endl; //@@
#endif      
    diag_c[1] = cdiag;
    how_many_equal.resize(pseudos_c0.size());
    int ic_to_compare = 0;
    for (int ip=0; ip<pseudos_c0.size(); ++ip){
      //how_many_equal(ip) = 0; //@@ we should avoid the current diagram, because we want to change
      dpseudo = pseudos_c0[ip]-pseud_c[ic_to_compare];
      how_many_equal(ip) = count_zeros( dpseudo );
      //}//@@
#ifdef _DEBUG      
      cout<<ip<<" p="<< pseudos_c0[ip] << " has " << how_many_equal(ip) << " equal with po=" << pseud_c[ic_to_compare]<<endl;
#endif      
    }
    deque<int> candidates1 = which_have_same_equals_but_not_same(how_many_equal,2*Norder);  //@@    
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates1 << endl;
#endif    
    iwhich_t1 = candidates1[ static_cast<int>(candidates1.size()*drand()) ];
    pseud_c[1] = pseudos_c0[iwhich_t1];
    dlttyp_c[1] = Delta_types_c0(iwhich_t1,Range::all());
    sgn_t[1] = csgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t1=" << iwhich_t1 << " p_c1=" << pseud_c[1] << endl;
#endif
    //Determine_all_diagrams_of_this_topology2(pseudos_trial, Delta_types_trial, diag, non_zero, Dlt, Fi_Fip);
    Determine_diagrams_of_this_topology_with(q_known, pseudos_trial, Delta_types_trial, diag, non_zero, Dlt, Fi_Fip); //@@
    auto t_3 = high_resolution_clock::now();
    t_topol += t_3-t_2;
#ifdef _DEBUG
    cout<<ii<<" Number of diagrams in the new topology="<<pseudos_trial.size()<<endl;
#endif
    int icase=1; //icase = (drand() < P_first_jump) ? 1 : 2;
#ifdef _DEBUG    
    cout << "icase=" << icase << " cdiag="<<cdiag<<" diag="<<diag << " p_c0=" << cpseudo << " csgn_dlt=" << csgn_dlt << " sgn_dlt=" << sgn_dlt << endl;
#endif    
    diag_c[2] = diag;
    how_many_equal.resize(pseudos_trial.size());
    ic_to_compare = (icase==1) ? 1 : 0;
    for (int ip=0; ip<pseudos_trial.size(); ++ip){
      if (jump){
	for (int j=0; j<2*Norder; ++j){
	  int jn = o2n(j);
	  if (jn>=0) dpseudo(j) = pseudos_trial[ip](jn)-pseud_c[ic_to_compare](j);
	  else dpseudo(j) = -1;
	}
      }else{
	dpseudo = pseudos_trial[ip]- pseud_c[ic_to_compare];
      }
      how_many_equal(ip) = count_zeros( dpseudo );
#ifdef _DEBUG      
      cout<<setw(2)<<ip<<" p="<< pseudos_trial[ip] << " has " << how_many_equal(ip) << " equal with po=" << pseud_c[ic_to_compare]<<endl;
#endif      
    }
    deque<int> candidates2 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates2 << endl;
#endif
    iwhich_t2 = candidates2[ static_cast<int>(candidates2.size()*drand()) ];
    pseud_c[2] = pseudos_trial[iwhich_t2];
    dlttyp_c[2] = Delta_types_trial(iwhich_t2,Range::all());
    sgn_t[2] = sgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t2=" << iwhich_t2 << " p_c2=" << pseud_c[2] << endl;
#endif
    diag_c[3] = diag;
    ic_to_compare = 2;
    for (int ip=0; ip<pseudos_trial.size(); ++ip){
      if (ip==iwhich_t2)
	how_many_equal(ip) = 0; // we should avoid the current diagram, because we want to change
      else{
	dpseudo2 = pseudos_trial[ip]-pseud_c[ic_to_compare];
	how_many_equal(ip) = count_zeros( dpseudo2 );
      }
#ifdef _DEBUG      
      cout<<setw(2)<<ip<<" p="<< pseudos_trial[ip] << " has " << how_many_equal(ip) << " equal with po=" << pseud_c[ic_to_compare]<<endl;
#endif      
    }
    deque<int> candidates3 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates3 << endl;
#endif    
    iwhich_t3 = candidates3[ static_cast<int>(candidates3.size()*drand()) ];
    pseud_c[3] = pseudos_trial[iwhich_t3];
    dlttyp_c[3] = Delta_types_trial(iwhich_t3,Range::all());
    sgn_t[3] = sgn_dlt;
#ifdef _DEBUG    
    cout << "iwhich_t3=" << iwhich_t3 << " p_c3=" << pseud_c[3] << endl;
#endif	  
    auto t_4 = high_resolution_clock::now();
    t_choose += t_4-t_3;
    
    P_eq = 1.,1.,1.,1.;
    
    double V_0=pseudos_c0.size()*Vol_diagrams0(Norder-N_min_order), V_1=pseudos_trial.size()*Vol_diagrams0(tNorder-N_min_order);
    if (jump){
      cout<<" jump Norder="<<Norder<<"->"<<tNorder<<" V_0="<<V_0<<" V_1="<<V_1<<" N0="<<pseudos_c0.size()<<" N1="<<pseudos_trial.size()<<endl;
    }
    //P_pr[0] = 1.0/(V_0*V_1*V_1);
    //P_pr[1] = 1.0/(V_0*V_1*V_1);
    //P_pr[2] = 1.0/(V_0*V_0*V_1);
    //P_pr[3] = 1.0/(V_0*V_0*V_1);
    P_pr[0] = 0.5*V_0/(V_0+V_1);
    P_pr[1] = P_pr[0];
    P_pr[2] = 0.5*V_1/(V_0+V_1);
    P_pr[3] = P_pr[2];
    P_tree=0;
    for (int ic=0; ic<4; ++ic){
      double P0 = P_eq[ic]*P_pr[ic];
      P_tree += P0;
      P_acc[ic] = P0;
    }
    P_acc *= 1/P_tree;
#ifdef _DEBUG    
    cout << "P_tree="<<P_tree<<" P_acc=" << P_acc <<endl;
#endif
    for (int ic=0; ic<4; ++ic){
      sgn_t[ic] *=  SgnMElements(diag_c[ic], pseud_c[ic], Dlt);
      auto [new_diag,s] = Different3(diag_c[ic], pseud_c[ic], P_acc[ic], stored2);
      if (new_diag)
	cout <<ii<<") Final-Diag "<< diag_c[ic] << " with p=" << pseud_c[ic] << " is new diagram " << s << " wgh=" << P_acc[ic] << " sgn=" << sgn_t[ic] << endl;
      else
	cout <<ii<<") Final-Diag "<< diag_c[ic] << " with p=" << pseud_c[ic] << " is equivalent to " << s << " wgh=" << P_acc[ic] << " sgn=" << sgn_t[ic] << endl;
    }
    auto t_5 = high_resolution_clock::now();
    t_diff += t_5-t_4;
    // meassure
    {
      for (int ic=0; ic<4; ++ic){
	int norder = diag_c[ic].extent(0);
	Nmeassure(norder) += P_acc[ic];
	asign(norder) += sgn_t[ic]*P_acc[ic];
	for (int it=0; it<pseud_c[ic].size(); ++it){ //### Norder
	  int atomic_state = pseud_c[ic](it);
	  dSm(atomic_state) += P_acc[ic];
	}
      }
    }
    double xir = drand();
#ifdef _DEBUG    
    cout << "xir=" << xir << endl;
#endif    
    if (xir > 1-P_acc[3] ){ //@@
      // accept c3
      iacc = 3;
      if (jump!=0){
	cdiag.resize(tNorder,2);
	cpseudo.resize(2*tNorder);
	cDlt_type.resize(tNorder);
      }
      cdiag = diag_c[3].copy();
      cpseudo = pseud_c[3].copy();
      cDlt_type = dlttyp_c[3].copy();
      //iwhich_c0 = iwhich_t3; //@@
      //copy(pseudos_c0,pseudos_trial); pseudos_trial.clear(); //@@
      //copy(Delta_types_c0,Delta_types_trial); Delta_types_trial.clear(); //@@
      csgn_dlt = sgn_dlt;
      skl = tskl;
      cout << ii<<") Choosing c3 icase=" << icase << " p=" << cpseudo << " diag=" <<cdiag<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      if (jump!=0){
	Norder = tNorder;
	dpseudo.resize(2*Norder);
	diag_c[0].resize(Norder,2);
	pseud_c[0].resize(2*Norder);
	dlttyp_c[0].resize(Norder);
	diag_c[1].resize(Norder,2);
	pseud_c[1].resize(2*Norder);
	dlttyp_c[1].resize(Norder);	  
      }
    }else if (xir <= 1-P_acc[3] and xir > P_acc[0]+P_acc[1]){
      // accept c2
      iacc = 2;
      if (jump!=0){
	cdiag.resize(tNorder,2);
	cpseudo.resize(2*tNorder);
	cDlt_type.resize(tNorder);
      }
      cdiag = diag_c[2].copy();
      cpseudo = pseud_c[2].copy();
      cDlt_type = dlttyp_c[2].copy();
      //iwhich_c0 = iwhich_t2; //@@
      //copy(pseudos_c0,pseudos_trial); pseudos_trial.clear(); //@@
      //copy(Delta_types_c0,Delta_types_trial); Delta_types_trial.clear(); //@@
      csgn_dlt = sgn_dlt;
      skl = tskl;
      cout <<ii<<") Choosing c2 icase=" << icase << " p=" << cpseudo << " diag=" <<cdiag<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      if (jump!=0){
	Norder = tNorder;
	dpseudo.resize(2*Norder);
	diag_c[0].resize(Norder,2);
	pseud_c[0].resize(2*Norder);
	dlttyp_c[0].resize(Norder);
	diag_c[1].resize(Norder,2);
	pseud_c[1].resize(2*Norder);
	dlttyp_c[1].resize(Norder);	  
      }
    }else if (xir <= P_acc[0]+P_acc[1] and xir > P_acc[0]){
      iacc = 1;
      cpseudo = pseud_c[1].copy();
      cDlt_type = dlttyp_c[1].copy();
      //iwhich_c0 = iwhich_t1;
      tskl = skl;
      cout<<ii<<") Choosing c1 icase="<<icase<<"  p="<<cpseudo<<" diag="<<cdiag;
      cout<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      if (jump!=0){
	diag.resize(Norder,2);
	dpseudo2.resize(2*Norder);
	diag_c[2].resize(Norder,2);
	pseud_c[2].resize(2*Norder);
	dlttyp_c[2].resize(Norder);
	diag_c[3].resize(Norder,2);
	pseud_c[3].resize(2*Norder);
	dlttyp_c[3].resize(Norder);
      }
    } else if (xir <= P_acc[0]){ // c0 -> backtracking
      iacc = 0;
      cout <<ii<<") ************************* BACKTRACKING ***************" << endl;
      //tskl = skl;
      tskl.copy(skl);
      if (jump!=0){
	diag.resize(Norder,2);
	dpseudo2.resize(2*Norder);
	diag_c[2].resize(Norder,2);
	pseud_c[2].resize(2*Norder);
	dlttyp_c[2].resize(Norder);
	diag_c[3].resize(Norder,2);
	pseud_c[3].resize(2*Norder);
	dlttyp_c[3].resize(Norder);
      }
    }
    idiag += 1;
    auto t_6 = high_resolution_clock::now();
    t_next += t_6-t_5;
    if (idiag >= Ndiag) break;
  }
  //cout << "how_many_pgp=" << how_many_pgp/Ndiag/(2*Norder) << endl;

  for (auto itt=stored2.begin(); itt!=stored2.end(); ++itt){
    auto [diag,pseud] = String2diagp(itt->first);
    stored.push_back(make_pair(diag,pseud));
    stored_count.push_back( itt->second );
    cout<<itt->first<<" "<<itt->second<<" diag="<<diag<<" pseud="<<pseud<<endl;
  }
  
  vector<int> sfc(stored.size());
  for (int i=0; i<stored.size(); ++i)
    sfc[i] = SymmetryFactor2(stored[i]);

  deque<Array<int,1>> pseudos;
  Array<int,2> Delta_types;
  diag.resize(1,2); diag = 1,0;
  Determine_all_diagrams_of_this_topology(pseudos,Delta_types, diag,non_zero,Dlt);
  double N00 = pseudos.size();
  Array<double,1> Vol_diagrams(N_max_order-N_min_order+1), Vol_diagrams_jump(N_max_order-N_min_order+1);
  Vol_diagrams=0; Vol_diagrams_jump=0;
  for (int i=0; i<stored.size(); ++i){
    Array<int,2> diag(stored[i].first);
    int Norder = diag.extent(0);
    Vol_diagrams(Norder-N_min_order) += 1.0/sfc[i];
    Determine_all_diagrams_of_this_topology(pseudos,Delta_types, diag,non_zero,Dlt);
    Vol_diagrams_jump(Norder-N_min_order) += N00/(sfc[i]*pseudos.size());
  }
  
  Array<double,1> rr(N_max_order+1); rr = 0;
  Array<int,1> ndg(N_max_order+1); ndg=0;
  double rrt=0, rrv=0;
  for (int i=0; i<stored.size(); ++i){
    int Norder = stored[i].first.extent(0);
    ndg(Norder)  += 1;
    double xt = sfc[i]*stored_count[i];
    rr(Norder) += xt;
    rrt += xt;
    //rrv += xt*Vol_diagrams0(Norder-N_min_order);
  }
  for (int i=N_min_order; i<=N_max_order; ++i) rr(i) /= ndg(i);
  rrt /= stored.size();
  //rrv /= stored.size();
  cout<<"rr="<<rr<<" rrt="<<rrt<</*" rrv="<<rrv<<*/" sfc="<<sfc<<endl;

  Array<double,1> w_order(N_max_order+1); w_order=0;
  Array<double,1> deviation1(N_max_order+1); deviation1=0;
  double /*deviation3=0,*/ deviation2=0;
  cout <<"stored diags="<<endl;
  for (int iorder=N_min_order; iorder<=N_max_order; ++iorder){
    for (int i=0; i<stored.size(); ++i){
      int Norder = stored[i].first.extent(0);
      if (iorder!=Norder) continue;
      double xt = stored_count[i]*sfc[i];
      double x1 = xt/rr(Norder);
      double x2 = xt/rrt;
      double yt = xt*Vol_diagrams(Norder-N_min_order);
      //double x3 = yt/rrv;
      w_order(Norder) += x2;
      deviation1(Norder) += abs(x1-1);
      deviation2 += abs(x2-1);
      //deviation3 += abs(x3-1);

      cdiag = stored[i].first;
      cpseudo = stored[i].second;
      /*
	int sgn_dlt = SignSort(cdiag);
	TinyVector<int,2> sgn_mmbw = SignMatrixElements(cdiag, cpseudo, Dlt);
	int sgn_tot = sgn_dlt*sgn_mmbw[0];
	int sgn_cur = sgn_dlt*sgn_mmbw[0]*sgn_mmbw[1];
      */
      cout<<setw(3)<<i<<" "<<setw(2)<<Norder<<" "<<setw(8)<<setprecision(3)<<xt<<" ";
      cout<<setw(2)<<sfc[i]<<"  ";
      cout<<left<<setw(8)<<setprecision(6)<<x1<<" "<<setw(8)<<x2<<" "<<setw(8)/*<<x3*/<<" "<<right;
      cout<<"d="<<stored[i].first<<" p=" << stored[i].second<<endl;
    }
  }
  for (int i=N_min_order; i<=N_max_order; ++i) w_order(i) /= ndg(i);
  for (int i=N_min_order; i<=N_max_order; ++i) cout<<"w["<<i<<"]=" <<w_order(i) << endl;
  
  //deviation3 /= stored.size();
  deviation2 /= stored.size();
  cout<<"Vol_diagrams="<<Vol_diagrams<<endl;
  cout<<"Vol_diagrams jump="<<Vol_diagrams_jump<<endl;
  
  for (int i=N_min_order; i<=N_max_order; ++i) deviation1(i) /= ndg(i);
  cout<<"deviation1=[";
  for (int i=N_min_order; i<=N_max_order; ++i) cout<<deviation1(i)<<",";
  cout<<"] ";
  cout << " deviation2="<<deviation2<</*" deviation3=" << deviation3 <<*/ endl;

  cout<<"<s>=[";
  for (int iorder=N_min_order; iorder<=N_max_order; ++iorder) cout<< asign(iorder)/Nmeassure(iorder) <<",";
  cout<<"]"<<endl;

  /*
  cout << "Before: Pbin="<<Pbin*(V0norm/Pnorm)<<" dSm="<< Array<double,1>((dSm*(V0norm/Pnorm))) <<endl;
  Pbin *= (V0norm/Pnorm)*Vol_diagrams;
  dSm  *= (V0norm/Pnorm)*Vol_diagrams;
  cout << "After Pbin="<< Pbin <<" dSm="<<dSm<<endl;
  */
  cout<<"t_next="<<t_next.count()<<" t_diff="<<t_diff.count()<<" t_choose="<<t_choose.count()<<" t_topol="<<t_topol.count()<<" t_skel="<<t_skel.count()<<endl;
}

template<class RNG>
void GenerateDiagrams20_(Array<int,2>& i_diag,const vector<TinyVector<int,4>> non_zero,const Array<int,4>& Dlt, const Array<int,2>& Fi_Fip, RNG& drand,int Ndiag=10000){
  //#define _DEBUG  
  cout.precision(12);
  int N_min_order=1, N_max_order=i_diag.extent(0);
  Array<double,1> Vol_diagrams0(N_max_order-N_min_order+1); Vol_diagrams0=1;
  Array<double,1> Vol_diagrams1(N_max_order-N_min_order+1); Vol_diagrams1=1;
  {
    Array<double,1> Vol_d0(8); Vol_d0 = 1,4,2,16./3.,27,198.4,1886.67,21846.607;
    Array<double,1> Vol_d1(8); Vol_d1 = 1,4,2,16./3.,36,310.4,3424.00,45427.857;
    TinyVector<double,3> pv(0.15760576,0.5415707,-1.3749251);
    for (int n=N_min_order; n<=N_max_order; ++n){
      Vol_diagrams0(n-N_min_order) = (n<8) ? Vol_d0(n) : pv[0]*n*n+pv[1]*n+pv[2];
      Vol_diagrams1(n-N_min_order) = (n<8) ? Vol_d1(n) : pv[0]*n*n+pv[1]*n+pv[2];
    }
  }
  //Array<double,1> Vol_diagrams0(max(N_max_order+1,8)); Vol_diagrams0=1;
  //Vol_diagrams0 = 1,4,2,16./3.,27,198.4,1886.67,21846.607;
  double P_first_jump = 0.5;

  deque<Array<int,1>> pseudos_c0, pseudos_trial;
  Array<int,2> Delta_types_c0, Delta_types_trial;
  
  int Norder=i_diag.extent(0);
  Array<int,2> diag(Norder,2), cdiag(Norder,2), ndiag(Norder,2);
  diag = i_diag;
  Skeletons skl(diag);
  bool Reducible = skl.IsReducible();
  {
    int sgn_dlt=1;
    FindNextSkeleton(sgn_dlt, diag, skl,drand);
  }
  Skeletons tskl(diag);
  
  TinyVector<Array<int,2>,4> diag_c;
  TinyVector<Array<int,1>,4> pseud_c;
  TinyVector<Array<int,1>,4> dlttyp_c;
  TinyVector<int,4> sgn_t;
  for (int i=0; i<4; ++i){
    diag_c[i].resize(Norder,2);
    pseud_c[i].resize(2*Norder);
    dlttyp_c[i].resize(Norder);
  }
  Array<int,1> cpseudo(2*Norder), dpseudo(2*Norder), dpseudo2(2*Norder), cDlt_type(Norder);

  vector<int> p_known({0});  //$$
  vector<int> d_known({0});  //$$
  //Determine_diagrams_of_this_topology_with(p_known, pseudos_c0, Delta_types_c0, diag, non_zero, Dlt, Fi_Fip); //@@
  Determine_diagrams_of_this_topology_with_option1(d_known,p_known,pseudos_c0,Delta_types_c0,diag,non_zero,Dlt,Fi_Fip);
  
  { //@@
    //Determine_all_diagrams_of_this_topology2(pseudos_c0, Delta_types_c0, diag, non_zero, Dlt, Fi_Fip);
    int iwhich_c0 = static_cast<int>(drand()*pseudos_c0.size());
    cpseudo = pseudos_c0[iwhich_c0].copy(); 
    cDlt_type = Delta_types_c0(iwhich_c0,Range::all());
  } //@@
  
  Array<int,1> how_many_equal(pseudos_c0.size()); how_many_equal = 0;
  for (int ip=0; ip<pseudos_c0.size(); ++ip){
    dpseudo = pseudos_c0[ip]-cpseudo;
    how_many_equal(ip) = count_zeros( dpseudo );
    cout<<ip<<" p="<< pseudos_c0[ip] << " has " << how_many_equal(ip) << " equal  cDlt=" << cDlt_type << endl;
  }
  cout << "Starting with diag=" <<diag<< " and p=" << cpseudo <<" cDlt=" << cDlt_type << endl;
  cdiag = diag;
  int csgn_dlt = SignSort(cdiag);

  int idiag=0;
  int Ndd=4; //$$
  double Pnorm=0, Pbin=0, Pnorm_sum=0, Pbin_sum=0;
  Array<double,1> Nmeassure(N_max_order+1), asign(N_max_order+1);
  Nmeassure=0; asign=0;
  map<string,double> stored2;
  deque<pair<Array<int,2>,Array<int,1>>> stored;
  deque<double> stored_count;
  deque<double> stored_sign;
  Array<int,1> o2n;
  int tNorder=Norder;
  tskl = skl;
  Array<double,1> dSm(4); dSm=0;
  double V0norm=1.0;
  TinyVector<double,4> P_acc, P_pr, P_eq;
  double P_tree=0;
  int iacc=0;
  int iwhich_t1, iwhich_t2, iwhich_t3;
  fsec t_next=zero_t, t_diff=zero_t, t_choose=zero_t, t_topol=zero_t, t_skel=zero_t;
  double how_many_pgp=0;
  for (int ii=0; ii<100000000; ++ii){
    auto t_1 = high_resolution_clock::now();
    if (diag_c[0].extent(0)!=cdiag.extent(0)){
      diag_c[0].resize(Norder,2);
      pseud_c[0].resize(2*Norder);
      dlttyp_c[0].resize(Norder);
    }
    diag_c[0] = cdiag;
    pseud_c[0] = cpseudo;
    dlttyp_c[0] = cDlt_type;
    sgn_t[0] = csgn_dlt;
    //
    Norder = cdiag.extent(0);
    int sgn_dlt = csgn_dlt;
    int jump=0;
    double xi=drand();
    if (xi < 1./3.) jump = 1;
    else if ( xi < 2./3.) jump = -1;
    if (jump==1 and Norder >= N_max_order) jump=0;
    if (jump==-1 and Norder <= N_min_order) jump=0;
#ifdef _DEBUG    
    cout<<ii<<") diag="<<diag<<" cdiag="<<cdiag<<" jump="<<jump<<endl;
#endif    
    if (jump==1){
      tNorder = Norder+1;
      TinyVector<pair<double,int>,2>  hybn = GenerateRandomHybridization(Norder,TinyVector<double,2>(drand(),drand()));
      diag.resize(tNorder,2);
      diag = InsertHybridization(cdiag,hybn,o2n);
      tskl.Resize_Initialize(diag);
      if (tskl.IsReducible()) FindNextSkeleton(sgn_dlt, diag, tskl,drand);
      sgn_dlt = SignSort(diag);
    }else if (jump==-1){
      int iremove = Norder*drand();
      tNorder = Norder-1;
      diag.resize(tNorder,2);
      diag = RemoveHybridization(cdiag, iremove,o2n);
      tskl.Resize_Initialize(diag);
      if (tskl.IsReducible()) FindNextSkeleton(sgn_dlt, diag, tskl,drand);
      sgn_dlt = SignSort(diag);
    }else{
      tNorder = Norder;
      diag = cdiag;
      FindNextSkeleton(sgn_dlt, diag, tskl,drand);
    }
    TinyVector<int,2> Morder(Norder,tNorder);
    //if (jump) cout<<"diag="<<diag<<" o2n="<<o2n<<endl;
    // All diagrams with all possible pseudoparticle arrangements are produced by recursion
    auto t_2 = high_resolution_clock::now();
    t_skel += t_2-t_1;
    //auto [p_known, q_known] = GiveCommonRandomPseudoPropagator(TinyVector<int,2>(Norder,tNorder),cpseudo,o2n,drand());

    int n_cmn=max(1,Norder/2);
    double xrand=drand();
    if (jump!=0) cout<<"Morder="<<Morder<<" n_cmn="<<n_cmn<<" o2n="<<o2n<<endl;
    auto[p_known,d0_known,d1_known]= GiveCommonRandomPseudoPropagator_of_size(n_cmn,Morder,cpseudo,o2n,xrand);
    n_cmn = p_known.size();
    if (n_cmn==0){
      int k=0;
      for (k=0; k<Norder; ++k) if (o2n(k)>=0) break;
      d0_known=vector<int>({k});
      d1_known=vector<int>({o2n(k)});
      p_known=vector<int>({cpseudo(k)});
    }
    cout<<"n_cmn="<<n_cmn<<" d0_known="<<d0_known<<" d1_known="<<d1_known<<" p_known="<<p_known<<endl;   //@@
    Determine_diagrams_of_this_topology_with_option1(d0_known,p_known,pseudos_c0,Delta_types_c0,cdiag,non_zero,Dlt,Fi_Fip);
    Determine_diagrams_of_this_topology_with_option1(d1_known,p_known,pseudos_trial,Delta_types_trial,diag,non_zero,Dlt,Fi_Fip);
    cout<<"Number of diagrams in the old and new topology="<<pseudos_c0.size()<<"/"<<pseudos_trial.size()<<endl; //@@
    while (pseudos_c0.size()<1 or pseudos_trial.size()<1){
      n_cmn-=1;
      auto[p_known,d0_known,d1_known]= GiveCommonRandomPseudoPropagator_of_size(n_cmn,Morder,cpseudo,o2n,xrand);
      cout<<"d0_known="<<d0_known<<" d1_known="<<d1_known<<" p_known="<<p_known<<endl;
      Determine_diagrams_of_this_topology_with_option1(d0_known,p_known,pseudos_c0,Delta_types_c0,cdiag,non_zero,Dlt,Fi_Fip);
      Determine_diagrams_of_this_topology_with_option1(d1_known,p_known,pseudos_trial,Delta_types_trial,diag,non_zero,Dlt,Fi_Fip);
      cout<<"Number of diagrams in the old and new topology="<<pseudos_c0.size()<<"/"<<pseudos_trial.size()<<endl; //@@
    }
    auto t_3 = high_resolution_clock::now();
    t_topol += t_3-t_2;
#ifdef _DEBUG      
    for (int ip=0; ip<pseudos_c0.size(); ++ip) cout<<ip<<" p0="<< pseudos_c0[ip] << endl; //@@
    for (int ip=0; ip<pseudos_trial.size(); ++ip) cout<<ip<<" p1="<< pseudos_trial[ip] << endl; //@@
#endif
    int Ndd=4;
    if (pseudos_c0.size()==1) Ndd-=1; // tree is smaller now with only 3 entries.
    if (pseudos_trial.size()==1) Ndd-=1; // tree is smaller now with only 3 entries.
    
    int idc=1;
    if (pseudos_c0.size()>1){
      if (diag_c[idc].extent(0)!=Norder){
	diag_c[idc].resize(Norder,2);
	pseud_c[idc].resize(2*Norder);
	dlttyp_c[idc].resize(Norder);
      }
      if (dpseudo.size()!=2*Norder) dpseudo.resize(2*Norder);
      //
      diag_c[idc] = cdiag;
      how_many_equal.resize(pseudos_c0.size());
      int ic_to_compare = 0;
      for (int ip=0; ip<pseudos_c0.size(); ++ip){
	dpseudo = pseudos_c0[ip]-pseud_c[ic_to_compare];
	how_many_equal(ip) = count_zeros( dpseudo );
#ifdef _DEBUG      
	cout<<ip<<" p="<< pseudos_c0[ip] << " has " << how_many_equal(ip) << " equal with po=" << pseud_c[ic_to_compare]<<endl;
#endif      
      }
      deque<int> candidates1 = which_have_same_equals_but_not_same(how_many_equal,2*Norder);  //@@    
#ifdef _DEBUG    
      cout << "candidates=" <<  candidates1 << endl;
#endif    
      iwhich_t1 = candidates1[ static_cast<int>(candidates1.size()*drand()) ];
      pseud_c[idc] = pseudos_c0[iwhich_t1];
      dlttyp_c[idc] = Delta_types_c0(iwhich_t1,Range::all());
      sgn_t[idc] = csgn_dlt;
#ifdef _DEBUG    
      cout << "iwhich_t1=" << iwhich_t1 << " p_c1=" << pseud_c[1] << endl;
#endif
      idc +=1;
    }
    //int icase=1; //icase = (drand() < P_first_jump) ? 1 : 2;
#ifdef _DEBUG    
    cout<<"cdiag="<<cdiag<<" diag="<<diag << " p_c0=" << cpseudo << " csgn_dlt=" << csgn_dlt << " sgn_dlt=" << sgn_dlt << endl;
#endif
    if (diag_c[idc].extent(0)!=tNorder){
      diag_c[idc].resize(tNorder,2);
      pseud_c[idc].resize(2*tNorder);
      dlttyp_c[idc].resize(tNorder);
    }
    if (dpseudo.size()!=2*Norder) dpseudo.resize(2*Norder);
    diag_c[idc] = diag;
    how_many_equal.resize(pseudos_trial.size());
    int ic_to_compare = idc-1;
    for (int ip=0; ip<pseudos_trial.size(); ++ip){
      if (jump){
	for (int j=0; j<2*Norder; ++j){
	  int jn = o2n(j);
	  if (jn>=0) dpseudo(j) = pseudos_trial[ip](jn)-pseud_c[ic_to_compare](j);
	  else dpseudo(j) = -1;
	}
      }else{
	dpseudo = pseudos_trial[ip]- pseud_c[ic_to_compare];
      }
      how_many_equal(ip) = count_zeros( dpseudo );
#ifdef _DEBUG      
      cout<<setw(2)<<ip<<" p="<< pseudos_trial[ip] << " has " << how_many_equal(ip) << " equal with po=" << pseud_c[ic_to_compare]<<endl;
#endif      
    }
    deque<int> candidates2 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
    cout << "candidates=" <<  candidates2 << endl;
#endif
    iwhich_t2 = candidates2[ static_cast<int>(candidates2.size()*drand()) ];
    pseud_c[idc] = pseudos_trial[iwhich_t2];
    dlttyp_c[idc] = Delta_types_trial(iwhich_t2,Range::all());
    sgn_t[idc] = sgn_dlt;
    idc+=1;
#ifdef _DEBUG    
    cout << "iwhich_t2=" << iwhich_t2 << " p_c2=" << pseud_c[2] << endl;
#endif
    if (pseudos_trial.size()>1){
      if (diag_c[idc].extent(0)!=tNorder){
	diag_c[idc].resize(tNorder,2);
	pseud_c[idc].resize(2*tNorder);
	dlttyp_c[idc].resize(tNorder);
      }
      if (dpseudo2.size()!=2*tNorder) dpseudo2.resize(2*tNorder);
      diag_c[idc] = diag;
      ic_to_compare=idc-1;
      for (int ip=0; ip<pseudos_trial.size(); ++ip){
	if (ip==iwhich_t2)
	  how_many_equal(ip) = 0; // we should avoid the current diagram, because we want to change
	else{
	  dpseudo2 = pseudos_trial[ip]-pseud_c[ic_to_compare];
	  how_many_equal(ip) = count_zeros( dpseudo2 );
	}
#ifdef _DEBUG      
	cout<<setw(2)<<ip<<" p="<< pseudos_trial[ip] << " has " << how_many_equal(ip) << " equal with po=" << pseud_c[ic_to_compare]<<endl;
#endif      
      }
      deque<int> candidates3 = which_have_same_equals(how_many_equal);
#ifdef _DEBUG    
      cout << "candidates=" <<  candidates3 << endl;
#endif    
      iwhich_t3 = candidates3[ static_cast<int>(candidates3.size()*drand()) ];
      pseud_c[idc] = pseudos_trial[iwhich_t3];
      dlttyp_c[idc] = Delta_types_trial(iwhich_t3,Range::all());
      sgn_t[idc] = sgn_dlt;
      idc+=1;
#ifdef _DEBUG    
      cout << "iwhich_t3=" << iwhich_t3 << " p_c3=" << pseud_c[3] << endl;
#endif
    }
    if (Ndd!=idc) cout<<"ERROR we expected Ndd="<<Ndd<<" to be equal to idc="<<idc<<endl;
    
    auto t_4 = high_resolution_clock::now();
    t_choose += t_4-t_3;
    
    P_eq = 1.,1.,1.,1.;
    
    double V_0=pseudos_c0.size()*Vol_diagrams0(Norder-N_min_order), V_1=pseudos_trial.size()*Vol_diagrams0(tNorder-N_min_order);
    if (jump){
      cout<<" jump Norder="<<Norder<<"->"<<tNorder<<" V_0="<<V_0<<" V_1="<<V_1<<" N0="<<pseudos_c0.size()<<" N1="<<pseudos_trial.size()<<endl;
    }
    //P_pr[0] = 0.5*V_0/(V_0+V_1);
    //P_pr[1] = P_pr[0];
    //P_pr[2] = 0.5*V_1/(V_0+V_1);
    //P_pr[3] = P_pr[2];
    P_pr=0,0,0,0;
    idc=0;
    if (pseudos_c0.size()>1){
      P_pr[idc  ]=0.5*V_0;
      P_pr[idc+1]=0.5*V_0;
      idc+=2;
    }else{
      P_pr[idc]=V_0;
      idc+=1;
    }
    if (pseudos_trial.size()>1){
      P_pr[idc  ]=0.5*V_1;
      P_pr[idc+1]=0.5*V_1;
      idc+=2;
    }else{
      P_pr[idc]=V_1;
      idc+=1;
    }
    P_pr *= 1.0/sum(P_pr);
    cout<<"P_pr="<<P_pr<<endl;
    
    P_tree=0;
    P_acc=0,0,0,0;
    for (int ic=0; ic<Ndd; ++ic){
      double P0 = P_eq[ic]*P_pr[ic];
      P_tree += P0;
      P_acc[ic] = P0;
    }
    P_acc *= 1/P_tree;
#ifdef _DEBUG  
    for (int ic=0; ic<Ndd; ++ic)
      cout<<ii<<") ic="<<ic<<" d="<<diag_c[ic]<<" p="<<pseud_c[ic]<<" Pacc="<<P_acc[ic]<<endl;
#endif    
    for (int ic=0; ic<Ndd; ++ic){
      sgn_t[ic] *=  SgnMElements(diag_c[ic], pseud_c[ic], Dlt);
      auto [new_diag,s] = Different3(diag_c[ic], pseud_c[ic], P_acc[ic], stored2);
      if (new_diag)
	cout <<ii<<") Final-Diag "<< diag_c[ic] << " with p=" << pseud_c[ic] << " is new diagram " << s << " wgh=" << P_acc[ic] << " sgn=" << sgn_t[ic] << endl;
      else
	cout <<ii<<") Final-Diag "<< diag_c[ic] << " with p=" << pseud_c[ic] << " is equivalent to " << s << " wgh=" << P_acc[ic] << " sgn=" << sgn_t[ic] << endl;
    }
    auto t_5 = high_resolution_clock::now();
    t_diff += t_5-t_4;
    // meassure
    {
      for (int ic=0; ic<Ndd; ++ic){
	int norder = diag_c[ic].extent(0);
	Nmeassure(norder) += P_acc[ic];
	asign(norder) += sgn_t[ic]*P_acc[ic];
	for (int it=0; it<pseud_c[ic].size(); ++it){ //### Norder
	  int atomic_state = pseud_c[ic](it);
	  dSm(atomic_state) += P_acc[ic];
	}
      }
    }
    double xir = drand();
    int iacc = Which_Case(xir, P_acc);
#ifdef _DEBUG    
    cout <<"xir="<<xir<<" P_acc="<<P_acc<<" iacc="<<iacc<<endl;
#endif
    switch(iacc){
    case 3:
      if (jump!=0){
	cdiag.resize(tNorder,2);
	cpseudo.resize(2*tNorder);
	cDlt_type.resize(tNorder);
      }
      cdiag = diag_c[iacc].copy();
      cpseudo = pseud_c[iacc].copy();
      cDlt_type = dlttyp_c[iacc].copy();
      csgn_dlt = sgn_dlt;
      skl = tskl;
      if (jump!=0) Norder = tNorder;
      cout << ii<<") Choosing c3  p=" << cpseudo << " diag=" <<cdiag<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      break;
    case 2:
      if (jump!=0){
	cdiag.resize(tNorder,2);
	cpseudo.resize(2*tNorder);
	cDlt_type.resize(tNorder);
      }
      cdiag = diag_c[iacc].copy();
      cpseudo = pseud_c[iacc].copy();
      cDlt_type = dlttyp_c[iacc].copy();
      csgn_dlt = sgn_dlt;
      skl = tskl;
      if (jump!=0) Norder = tNorder;
      cout <<ii<<") Choosing c2  p=" << cpseudo << " diag=" <<cdiag<< " csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      break;
    case 1:
      if (pseudos_c0.size()==1){
	if (jump!=0){
	  cdiag.resize(tNorder,2);
	  cpseudo.resize(2*tNorder);
	  cDlt_type.resize(tNorder);
	}
	cdiag = diag_c[iacc].copy();
	cpseudo = pseud_c[iacc].copy();
	cDlt_type = dlttyp_c[iacc].copy();
	csgn_dlt = sgn_dlt;
	skl = tskl;
	if (jump!=0) Norder = tNorder;
      }else{
	cpseudo = pseud_c[iacc].copy();
	cDlt_type = dlttyp_c[iacc].copy();
	tskl = skl;
	if (jump!=0) diag.resize(Norder,2);
      }
      cout<<ii<<") Choosing c1  p="<<cpseudo<<" diag="<<cdiag<<" csgn_dlt=" << csgn_dlt << " cDlt=" << cDlt_type << endl;
      break;
    case 0:
      cout <<ii<<") ************************* BACKTRACKING ***************" << endl;
      tskl = skl;
      if (jump!=0) diag.resize(Norder,2);
      break;
    default:
      cout<<"ERROR It should not happen"<<endl; exit(0);
    }
    idiag += 1;
    auto t_6 = high_resolution_clock::now();
    t_next += t_6-t_5;
    if (idiag >= Ndiag) break;
  }
  //cout << "how_many_pgp=" << how_many_pgp/Ndiag/(2*Norder) << endl;

  for (auto itt=stored2.begin(); itt!=stored2.end(); ++itt){
    auto [diag,pseud] = String2diagp(itt->first);
    stored.push_back(make_pair(diag,pseud));
    stored_count.push_back( itt->second );
    cout<<itt->first<<" "<<itt->second<<" diag="<<diag<<" pseud="<<pseud<<endl;
  }
  
  vector<int> sfc(stored.size());
  for (int i=0; i<stored.size(); ++i)
    sfc[i] = SymmetryFactor2(stored[i]);

  deque<Array<int,1>> pseudos;
  Array<int,2> Delta_types;
  diag.resize(1,2); diag = 1,0;
  Determine_all_diagrams_of_this_topology(pseudos,Delta_types, diag,non_zero,Dlt);
  double N00 = pseudos.size();
  Array<double,1> Vol_diagrams(N_max_order-N_min_order+1), Vol_diagrams_jump(N_max_order-N_min_order+1);
  Vol_diagrams=0; Vol_diagrams_jump=0;
  for (int i=0; i<stored.size(); ++i){
    Array<int,2> diag(stored[i].first);
    int Norder = diag.extent(0);
    Vol_diagrams(Norder-N_min_order) += 1.0/sfc[i];
    Determine_all_diagrams_of_this_topology(pseudos,Delta_types, diag,non_zero,Dlt);
    Vol_diagrams_jump(Norder-N_min_order) += N00/(sfc[i]*pseudos.size());
  }
  
  Array<double,1> rr(N_max_order+1); rr = 0;
  Array<int,1> ndg(N_max_order+1); ndg=0;
  double rrt=0, rrv=0;
  for (int i=0; i<stored.size(); ++i){
    int Norder = stored[i].first.extent(0);
    ndg(Norder)  += 1;
    double xt = sfc[i]*stored_count[i];
    rr(Norder) += xt;
    rrt += xt;
    //rrv += xt*Vol_diagrams0(Norder-N_min_order);
  }
  for (int i=N_min_order; i<=N_max_order; ++i) rr(i) /= ndg(i);
  rrt /= stored.size();
  //rrv /= stored.size();
  cout<<"rr="<<rr<<" rrt="<<rrt<</*" rrv="<<rrv<<*/" sfc="<<sfc<<endl;

  Array<double,1> w_order(N_max_order+1); w_order=0;
  Array<double,1> deviation1(N_max_order+1); deviation1=0;
  double /*deviation3=0,*/ deviation2=0;
  cout <<"stored diags="<<endl;
  for (int iorder=N_min_order; iorder<=N_max_order; ++iorder){
    for (int i=0; i<stored.size(); ++i){
      int Norder = stored[i].first.extent(0);
      if (iorder!=Norder) continue;
      double xt = stored_count[i]*sfc[i];
      double x1 = xt/rr(Norder);
      double x2 = xt/rrt;
      double yt = xt*Vol_diagrams(Norder-N_min_order);
      //double x3 = yt/rrv;
      w_order(Norder) += x2;
      deviation1(Norder) += abs(x1-1);
      deviation2 += abs(x2-1);
      //deviation3 += abs(x3-1);

      cdiag = stored[i].first;
      cpseudo = stored[i].second;
      /*
	int sgn_dlt = SignSort(cdiag);
	TinyVector<int,2> sgn_mmbw = SignMatrixElements(cdiag, cpseudo, Dlt);
	int sgn_tot = sgn_dlt*sgn_mmbw[0];
	int sgn_cur = sgn_dlt*sgn_mmbw[0]*sgn_mmbw[1];
      */
      cout<<setw(3)<<i<<" "<<setw(2)<<Norder<<" "<<setw(8)<<setprecision(3)<<xt<<" ";
      cout<<setw(2)<<sfc[i]<<"  ";
      cout<<left<<setw(8)<<setprecision(6)<<x1<<" "<<setw(8)<<x2<<" "<<setw(8)/*<<x3*/<<" "<<right;
      cout<<"d="<<stored[i].first<<" p=" << stored[i].second<<endl;
    }
  }
  for (int i=N_min_order; i<=N_max_order; ++i) w_order(i) /= ndg(i);
  for (int i=N_min_order; i<=N_max_order; ++i) cout<<"w["<<i<<"]=" <<w_order(i) << endl;
  
  //deviation3 /= stored.size();
  deviation2 /= stored.size();
  cout<<"Vol_diagrams="<<Vol_diagrams<<endl;
  cout<<"Vol_diagrams jump="<<Vol_diagrams_jump<<endl;
  
  for (int i=N_min_order; i<=N_max_order; ++i) deviation1(i) /= ndg(i);
  cout<<"deviation1=[";
  for (int i=N_min_order; i<=N_max_order; ++i) cout<<deviation1(i)<<",";
  cout<<"] ";
  cout << " deviation2="<<deviation2<</*" deviation3=" << deviation3 <<*/ endl;

  cout<<"<s>=[";
  for (int iorder=N_min_order; iorder<=N_max_order; ++iorder) cout<< asign(iorder)/Nmeassure(iorder) <<",";
  cout<<"]"<<endl;

  /*
  cout << "Before: Pbin="<<Pbin*(V0norm/Pnorm)<<" dSm="<< Array<double,1>((dSm*(V0norm/Pnorm))) <<endl;
  Pbin *= (V0norm/Pnorm)*Vol_diagrams;
  dSm  *= (V0norm/Pnorm)*Vol_diagrams;
  cout << "After Pbin="<< Pbin <<" dSm="<<dSm<<endl;
  */
  cout<<"t_next="<<t_next.count()<<" t_diff="<<t_diff.count()<<" t_choose="<<t_choose.count()<<" t_topol="<<t_topol.count()<<" t_skel="<<t_skel.count()<<endl;
}
