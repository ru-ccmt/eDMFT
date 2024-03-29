// @Copyright 2007-2017 Kristjan Haule
// 


void ClusterData::Check(function1D<NState>& praStates)
{
  function1D<Number> Projection(common::max_size);
  NState mstate(common::max_size,common::max_size), pstate(common::max_size,common::max_size);
  for (int ist=0; ist<praStates.size(); ist++){
    // calculating  psi^+_{ib1} psi_{ib2}
    for (int ib2=0; ib2<N_flavors; ib2++){
      mstate.apply(FM(2*ib2+1), Fi(2*ib2+1), praStates[ist]);  // applying c_{b2}|nstate> operator
      for (int ib1=0; ib1<N_flavors; ib1++){
	pstate.apply(FM(2*ib1  ), Fi(2*ib1  ), mstate);  // applying <pstate|c^+_{b1}|mstate><mstate|c_{b2}|nstate> operator
	if (pstate.istate != praStates[ist].istate) continue;
	Number mm = pstate.Project_to_Pra(praStates[ist], Projection);

	if (ib1==ib2){
	  cout<<"istate="<<praStates[ist].istate<<" b="<<ib1<<" Projection= ";
	  for (int m=0; m<msize(ist+1); m++) cout<<Projection[m].dbl()<<" ";
	  cout<<endl;
	}
      }
    }
  }
}

void InterpolateProbabilities(vector<function1D<double> >& Prob, const function2D<double>& AProb, int cnsize,
			      const function1D<int>& cmsize_, const function1D<int>& msize_, const function1D<int>& cindx)
{
  Prob.resize(cnsize);
  // Interpolate probabilities to the large atomic base
  for (int j=0; j<cnsize; j++){
    int isize = cmsize_[j+1];
    Prob[j].resize(isize);
    Prob[j]=0;
    int indx = cindx[j+1];
    if (indx==0) continue;
    for (int m=0; m<msize_[indx]; m++) Prob[j][m] = AProb[indx-1][m];
  }
}
void ClusterData::HB3(double beta, double mu, double U, const mesh1D& iom, const function2D<double>& AProb,
		      int nom, int aom, const function2D<dcomplex>& Delta, function2D<dcomplex>& Sigma,
		      const function1D<bool>& nonzero, const function1D<bool>& nonzerou, ostream& clog,
		      double sderiv=0.1, int nom_small=300)
{
  /*
  // Interpolates probabilities on the larger basis
  vector<function1D<double> > Prob(cnsize);
  InterpolateProbabilities(Prob, AProb, cnsize, cmsize_, msize_, cindx);
  // For this calculation the basis might be larger, as we need to take into account much more local states for high-frequency. Should not cut-off much...
  function2D<int>& F_i = cF_i;
  function2D<function2D<double> >& F_M = cF_M;
  function1D<int>& msize_ = cmsize_;
  int nsize = cnsize;
  int max_size = cmax_size;
  */
  const function2D<double>& Prob = AProb;
  
  // Previous praStates can not be used, because the basis can be larger. We will generate these states from scratch
  function1D<NState> praStates(nsize);
  for (int ist=0; ist<nsize; ist++){
    NState& m = praStates[ist];
    m.istate = ist+1;
    int msize = msize_[ist+1];
    m.M.resize(msize,msize);
    m.M = 0.0;
    for (int l=0; l<msize; l++) m.M(l,l)=1.0;
    m.exponent = 0;
  }
  
  function1D<Number> Projection(max_size);
  function2D<double> n_ab(N_flavors,N_flavors);
  
  n_ab=0.0;  
  NState nstate(max_size,max_size), mstate(max_size,max_size), qstate(max_size,max_size), pstate(max_size,max_size);
  // calculating the density matrix, which is defined by <psi^+_{ib1} psi_{ib2}>
  for (int ist=0; ist<praStates.size(); ist++){
    int istate = praStates[ist].istate;
    for (int ib2=0; ib2<N_flavors; ib2++){
      mstate.apply(F_M[2*ib2+1], F_i[2*ib2+1], praStates[ist]);  // applying c_{b2}|nstate> operator
      if (mstate.istate==0) continue;
      for (int ib1=0; ib1<N_flavors; ib1++){
	pstate.apply(F_M[2*ib1], F_i[2*ib1], mstate);  // applying <pstate|c^+_{b1}|mstate><mstate|c_{b2}|nstate> operator
	if (pstate.istate != praStates[ist].istate) continue;
	Number mm = pstate.Project_to_Pra(praStates[ist], Projection);
	
	double nij = 0;
	for (int m=0; m<msize_[istate]; m++) nij += Projection[m].dbl()*Prob[istate-1][m];
	n_ab(ib1,ib2) += nij;
      }
    }
  }
  function4D<double> nn_abcd(N_flavors,N_flavors,N_flavors,N_flavors);
  nn_abcd=0.0;
  // calculating  psi^+_{ib1} psi^+_{ib2} psi_{ib3} psi_{ib4}
  for (int ist=0; ist<praStates.size(); ist++){
    int istate = praStates[ist].istate;
    for (int ib4=0; ib4<N_flavors; ib4++){
      mstate.apply(F_M[2*ib4+1], F_i[2*ib4+1], praStates[ist]);  // applying c_{b4}|istate> operator
      if (mstate.istate==0) continue;
      for (int ib3=0; ib3<N_flavors; ib3++){
	nstate.apply(F_M[2*ib3+1], F_i[2*ib3+1], mstate);        // applying c_{b3} c_{b4}|istate>
	if (nstate.istate==0) continue;
	for (int ib2=0; ib2<N_flavors; ib2++){
	  qstate.apply(F_M[2*ib2], F_i[2*ib2  ], nstate);      // applying c^+_{b2} c_{b3} c_{b4}|istate>
	  if (qstate.istate==0) continue;
	  for (int ib1=0; ib1<N_flavors; ib1++){
	    pstate.apply(F_M[2*ib1], F_i[2*ib1  ], qstate);    // applying c^+_{b1} c^+_{b2} c_{b3} c_{b4}|istate>
	    if (pstate.istate != praStates[ist].istate) continue;
	    Number mm = pstate.Project_to_Pra(praStates[ist], Projection);
	    double nn = 0;
	    for (int m=0; m<msize_[istate]; m++) nn += Projection[m].dbl()*Prob[istate-1][m];
	    nn_abcd(ib1,ib2,ib3,ib4) += nn;
	  }
	}
      }
    }
  }
  clog<<"n_ab=[";
  double n_tot=0;
  for (int ib1=0; ib1<N_flavors; ib1++){
    clog<<n_ab(ib1,ib1)<<", ";
    n_tot += n_ab(ib1,ib1);
  }
  clog<<"]"<<endl;
  clog<<"n_tot="<<n_tot<<endl;

  /*
  clog<<"nn_abcd="<<endl;
  for (int ib1=0; ib1<N_flavors; ib1++){
    for (int ib2=0; ib2<N_flavors; ib2++){
      for (int ib3=0; ib3<N_flavors; ib3++){
	for (int ib4=0; ib4<N_flavors; ib4++){
	  if (fabs(nn_abcd(ib1,ib2,ib3,ib4))>1e-7)
	    clog<<ib1<<" "<<ib2<<" "<<ib3<<" "<<ib4<<" : "<<setw(8)<<nn_abcd(ib1,ib2,ib3,ib4)<<endl;
	}
      }
    }
  }
  */
  function4D<double> Uc(N_flavors,N_flavors,N_flavors,N_flavors);
  Construct_Ucf(Uc, U);
  // calculating the first two moments of self-energy
  vector<function1D<double> > S_oo(N_ifl), S_n1(N_ifl);
  for (int ifl=0; ifl<N_ifl; ifl++){
    S_oo[ifl].resize(N_baths(ifl));
    S_n1[ifl].resize(N_baths(ifl));
  }
  for (int ifl=0; ifl<N_ifl; ifl++){
    for (int ia=0; ia<ifl_dim[ifl]; ia++){
      for (int ib=0; ib<ifl_dim[ifl]; ib++){
	int ind_ab = tfl_index[ifl](ia,ib);
	int a = fl_from_ifl[ifl][ia];
	int b = fl_from_ifl[ifl][ib];
	double s_oo=0;
	for (int i=0; i<N_flavors; i++){
	  for (int j=0; j<N_flavors; j++){
	    s_oo += (Uc(a,i,j,b)-Uc(a,i,b,j))*n_ab(i,j);
	  }
	}
	S_oo[ifl][ind_ab] = s_oo;
      }
    }
  }

  for (int ifl=0; ifl<N_ifl; ifl++){
    for (int ia=0; ia<ifl_dim[ifl]; ia++){
      for (int ib=0; ib<ifl_dim[ifl]; ib++){
	int ind_ab = tfl_index[ifl](ia,ib);
	int a = fl_from_ifl[ifl][ia];
	int b = fl_from_ifl[ifl][ib];
	double s1=0;
	for (int i=0; i<N_flavors; i++){
	  for (int j=0; j<N_flavors; j++){
	    for (int p=0; p<N_flavors; p++){
	      for (int l=0; l<N_flavors; l++){
		double UU1=0, UU2=0;
		for (int k=0; k<N_flavors; k++){
		  UU1 += (Uc(a,i,j,k)-Uc(a,i,k,j))*(Uc(b,l,p,k)-Uc(b,l,k,p));
		  UU2 += Uc(a,k,j,l)*Uc(b,k,i,p);
		}
		s1 += UU1*(nn_abcd(p,i,j,l)-n_ab(p,l)*n_ab(i,j)) + UU2*nn_abcd(p,i,j,l);
	      }
	    }
	    for (int l=0; l<N_flavors; l++){
	      double UU=0;
	      for (int k=0; k<N_flavors; k++) UU += (Uc(a,i,j,k)-Uc(a,i,k,j))*Uc(b,l,j,k);
	      s1 += UU*n_ab(i,l);
	    }
	  }
	}
	S_n1[ifl][ind_ab] = s1;
      }
    }
  }

  function1D<double> s_oo(N_unique_fl), s_n1(N_unique_fl);
  function1D<int> deg(N_unique_fl);
  deg=0;
  s_oo=0;
  s_n1=0;
  for (int ifl=0; ifl<N_ifl; ifl++){
    for (int ib=0; ib<N_baths(ifl); ib++){
      int i_unique = bfl_index[ifl][ib];
      s_oo[i_unique] += S_oo[ifl][ib];
      s_n1[i_unique] += S_n1[ifl][ib];
      deg[i_unique]++;
    }
  }
  for (int i=0; i<N_unique_fl; i++){
    s_oo[i] *= 1./deg[i];
    s_n1[i] *= 1./deg[i];
  }
  
  for (int ifl=0; ifl<N_ifl; ifl++){
    clog<<"s_oo["<<ifl<<"]=";
    for (int ib=0; ib<N_baths(ifl); ib++) clog<<S_oo[ifl][ib]<<" ";
    clog<<endl;
    clog<<"S_1["<<ifl<<"]=";
    for (int ib=0; ib<N_baths(ifl); ib++) clog<<S_n1[ifl][ib]<<" ";
    clog<<endl;
  }
  clog<<"s_oo=[";
  for (int i=0; i<N_unique_fl; i++) clog<<s_oo[i]<<",";
  clog<<"]"<<endl;
  clog<<"s_n1=[";
  for (int i=0; i<N_unique_fl; i++) clog<<s_n1[i]<<",";
  clog<<"]"<<endl;
  
  // Create a small logarithimic mesh of 300 points
  mesh1D ioms;
  CreateLogMesh(1, nom_small, 1./beta, iom, ioms);
  
  // Interpolates between low energy QMC and high energy expanded self-energy
  for (int ifl=0; ifl<N_unique_fl; ifl++){
    if (!nonzerou[ifl]) continue;
    // First finding the value of self-energy at the end of sampling.
    dcomplex sb=0;               // value of the self-energy at the end of QMC sampling 
    for (int im=nom-aom; im<nom; im++) sb += Sigma[ifl][im];
    sb/=aom;
    double ob = iom[nom-(aom+1)/2];  // frequency at the end of QMC sampling

    // imag( s_n1[ifl]/(i*w-c) ) = imag( sb )
    double c2 = -s_n1[ifl]*ob/sb.imag()-sqr(ob);
    if (true){
      double c = sqrt(c2);
      // imaginary part
      for (int im=nom; im<iom.size(); im++) Sigma[ifl][im].imag() = -s_n1[ifl]*iom[im]/(sqr(iom[im])+c2);   // the rest is from Hubbard I
      // Real part
      for (int im=nom; im<iom.size(); im++) Sigma[ifl][im].real() = s_oo[ifl] + sqr(ob/iom[im])*(sb.real()-s_oo[ifl]);
    }else{
      // Imaginary part
      int is=nom;
      for (; is<iom.size()-1; is++){ // finding frequency at which the linear interpolation and the tail have similar enough derivative
	// S <=> s_oo[ifl]+s_n1[ifl]/(I*iom[is]);
	double Swi = -s_n1[ifl]/(iom[is]);
	double dSwi = s_n1[ifl]/(iom[is]*iom[is]);
	
	double df0 = (Swi-sb.imag())/(iom[is]-ob);
	double df1 = dSwi;
	if (fabs(df0-df1)<sderiv) break;
      }
      double se = -s_n1[ifl]/(iom[is]);
      double oe = iom[is];
      
      // imaginary part
      for (int im=nom; im<is; im++) Sigma[ifl][im].imag() = sb.imag() + (se-sb.imag())*(iom[im]-ob)/(oe-ob); // linear interpolation of the imaginary part
      for (int im=is; im<iom.size(); im++) Sigma[ifl][im].imag() = -s_n1[ifl]/(iom[im]);   // the rest is from Hubbard I
      // Real part
      for (int im=nom; im<iom.size(); im++) Sigma[ifl][im].real() = s_oo[ifl] + sqr(ob/iom[im])*(sb.real()-s_oo[ifl]);
    }
  }
  /*
  Gf.resize(N_unique_fl,iom.size());
  Gf_1.resize(N_unique_fl,iom.size());
  Gf=0;
  Gf_1=0;
  for (int im=0; im<iom.size(); im++){
    dcomplex iomega(0,iom[im]);
    for (int fl=0; fl<N_unique_fl; fl++)
      if (nonzerou[fl]) Gf_1(fl,im) = (iomega+mu)*Id[fl]-epsk[fl]-Sigma(fl,im)-Delta(fl,im);
    for (int ifl=0; ifl<N_ifl; ifl++)
      if (nonzero[ifl]) Inverse_Gf(ifl_dim[ifl], N_unique_fl, bfl_index[ifl], sign[ifl], conjg[ifl], im, Gf_1, Gf);
    for (int fl=0; fl<N_unique_fl; fl++) Gf(fl,im) /= fl_deg[fl];
  }
  
  ofstream tout("g_qmc.dat");
  for (int io=0; io<iom.size(); io++){
    tout<<iom[io]<<" ";
    for (int ifl=0; ifl<N_unique_fl; ifl++) tout<<Gf[ifl][io]<<" ";
    tout<<endl;
  }
  */
}

void ClusterData::Set_Old_Data_for_HB1(ostream& clog)
{
  //clog<<"WARN: Could not find data for HB1, hence using states from simulation. This is probably OK, but make sure no states were cut in the simulation. "<<endl;
  cnsize = nsize;
  cmsize_ = msize_;
  cmax_size = max_size;
  cNs=Ns;
  
  cF_i.resize(N_flavors,nsize+1);
  cF_M.resize(N_flavors,nsize+1);
  for (int i=1; i<=nsize; i++)
    for (int ib=0; ib<N_flavors; ib++){
      cF_i(ib,i)=F_i(2*ib,i);
      cF_M(ib,i)=F_M(2*ib,i);
    }
  cEne.resize(nsize+1);
  cindx.resize(nsize+1);
  for (int i=1; i<=nsize; i++){
    cEne[i].resize(msize(i));
    for (int m=0; m<msize(i); m++){
      cEne[i][m] = Enes[Eind(i,m)];
    }
    cindx[i]=i;
  }
}

void ClusterData::Read_for_HB1(ifstream& gout, double mu, double U, ostream& clog)
{
  cnsize=0;
  if (!QHB1) return;
  
  if (gout.eof()){
    Set_Old_Data_for_HB1(clog);
    return;
  }
    
  string str;
  clog<<"Continuing with HB1 after "<<Osize<<" operators read!"<<endl;
  int wNm, wN_ifl;
  gout>>wNm>>cnsize>>wN_ifl>>cmax_size;

  if (gout.eof()){
    Set_Old_Data_for_HB1(clog);
    return;
  }
 
  if (wNm!=Nm) {cerr<<"Nm="<<wNm<<" for HB1 and original Nm="<<Nm<<" are not the same"<<endl;exit(1);}
  if (wN_ifl!=N_ifl) {cerr<<"Nm="<<wN_ifl<<" for HB1 and original Nm="<<N_ifl<<" are not the same"<<endl;}

  //  cout<<"wNm="<<wNm<<" cnsize="<<cnsize<<" wN_ifl="<<wN_ifl<<endl;
  
  getline(gout,str); 
  getline(gout,str); // # ind   N   K   Jz size
  
  
  cNs.resize(cnsize+1);
  cmsize_.resize(cnsize+1);
  cindx.resize(cnsize+1);
  cF_i.resize(N_flavors,cnsize+1);
  cEne.resize(cnsize+1);
  cF_M.resize(N_flavors,cnsize+1);
  cF_i=0;
  
  for (int i=1; i<=cnsize; i++){
    int it1, in, ik, isize, tindx;
    double dsz;
    gout>>it1>>tindx>>in>>ik>>dsz>>isize;
    if (it1!=i){cerr<<"Something wrong in parsing cix file. it1!=i!"<<endl;}
    cNs[i] = in;
    cmsize_[i] = isize;
    cindx[i] = tindx;
    for (int ib=0; ib<N_flavors; ib++){
      int ist;
      gout>>ist;
      if (ist!=0){
	cF_i(ib,i)=ist;
      }
    }
    cEne[i].resize(isize);
    for (int is=0; is<isize; is++) gout>>cEne[i][is];
    double spin;
    for (int is=0; is<isize; is++) gout>>spin;
    gout.ignore(1000,'\n');
  }
  getline(gout,str); // # matrix elements
  
  for (int i=1; i<=cnsize; i++){
    for (int ib=0; ib<N_flavors; ib++){
      int it, jt, size1, size2;
      gout>>it>>jt>>size1>>size2;
      if (it!=i || jt!=cF_i(ib,i)) cerr<<"Something wrong reading cix file"<<endl;
      cF_M(ib,i).resize(size2,size1); // constructor changes i->ii
      //      int ii = cF_i(ib,i);
      for (int i1=0; i1<size1; i1++){
	for (int i2=0; i2<size2; i2++){
	  double m;
	  gout>>m;
	  cF_M(ib,i)(i2,i1)=m; // constructor only
	}
      }
      gout.ignore(1000,'\n');
    }
  }
  for (int i=1; i<=cnsize; i++){
    double dE = -cNs[i]*mu + 0.5*cNs[i]*(cNs[i]-1)*U;
    for (size_t m=0; m<cEne[i].size(); m++) cEne[i][m] += dE;
  }
}

void ClusterData::HB1(double beta, double mu, const mesh1D& iom, const function2D<double>& AProb,
		      int nom, int aom, const function2D<dcomplex>& Delta, function2D<dcomplex>& Sigma,
		      const function1D<bool>& nonzero, const function1D<bool>& nonzerou, ostream& clog,
		      double sderiv=0.1, int nom_small=300)
{
  // Interpolate probabilities to the large atomic base
  vector<function1D<double> > Prob(cnsize);
  for (int j=0; j<cnsize; j++){
    int isize = cmsize_[j+1];
    Prob[j].resize(isize);
    Prob[j]=0;
    int indx = cindx[j+1];
    if (indx==0) continue;
    for (int m=0; m<msize(indx); m++){
      Prob[j][m] = AProb[indx-1][m];
    }
  }
  // Create a small logarithimic mesh of 30 points
  mesh1D ioms;
  CreateLogMesh(1, nom_small, 1./beta, iom, ioms);

  // Computes atomic Green's function
  function2D<dcomplex> Gh(N_unique_fl,ioms.size());
  Gh=0;
  function1D<dcomplex> gh(ioms.size()); 


  cout<<"cnsize="<<cnsize<<" cEne.size="<<cEne.size()<<endl;

  for (int ifl=0; ifl<N_ifl; ifl++){
    if (! nonzero[ifl]) continue; // Some GF are projected out. Should be ignored
    for (int i1=0; i1<ifl_dim[ifl]; i1++){
      for (int i2=0; i2<ifl_dim[ifl]; i2++){
	int ia = fl_from_ifl(ifl,i1);
	int ib = fl_from_ifl(ifl,i2);
	int tfl = tfl_index[ifl][i1][i2];
	int bfl = bfl_index[ifl][tfl];
	int sgn = sign[ifl][tfl];
	int dcmp = conjg[ifl][tfl];
	gh=0;
	for (int i=1; i<=cnsize; i++){
	  int j = cF_i(ia,i);
	  int jb = cF_i(ib,i);
	  //int j_check = F_i(2*ia,i);
	  //int jb_check = F_i(2*ib,i);
	  if (j!=jb || j==0) continue;
	  for (int im=0; im<cmsize_[i]; im++){
	    double Pm = Prob[i-1][im];
	    for (int jm=0; jm<cmsize_[j]; jm++){
	      double Pn = Prob[j-1][jm];
	      double mm = cF_M(ia,i)(jm,im)*cF_M(ib,i)(jm,im)*(Pn+Pm);
	      double dE = cEne[j][jm]-cEne[i][im];
	      dcomplex ci = dcomplex(0,1);
	      double dE2 = sqr(dE);
	      double mdE = mm*dE;

	      double mm2 = cF_M(ia,i)(jm,im)*cF_M(ib,i)(jm,im);
	      //double mm2_check = F_M(2*ia,i)(jm,im)*F_M(2*ib,i)(jm,im);
	      //if (mm2 != mm2_check){
	      //cout<<" mm2 != mm2_check"<<endl;
	      //cout<<"ifl="<<ifl<<" i1="<<i1<<" i2="<<i2<<" i="<<i<<" im="<<im<<" Pm="<<Pm<<" mm="<<mm<<" dE="<<dE<<endl;
	      //}
	      for (int io=0; io<ioms.size(); io++){
		double ome = ioms[io];
		double den = 1/(ome*ome+dE2);
		gh[io].real() -= mdE*den;// gh = mm/(iom-dE)
		gh[io].imag() -= mm*ome*den;
	      }
	    }
	  }
	}
	for (int io=0; io<ioms.size(); io++){
	  dcomplex f = sgn*gh[io];
	  if (dcmp) f = f.conj();
	  Gh[bfl][io] += f;
	}
      }
    }
  }
  for (int ifl=0; ifl<N_unique_fl; ifl++) Gh[ifl] *= (1./fl_deg[ifl]);
  //  clog<<endl;
  
  // Computes Hubbard I self-energy
  function2D<dcomplex> Gf(N_unique_fl,ioms.size());
  function2D<dcomplex> Gf_1(N_unique_fl,ioms.size());
  Gf_1=0;
  vector<spline1D<dcomplex> > Sigh(N_unique_fl);
  for (size_t fl=0; fl<Sigh.size(); fl++){
    Sigh[fl].resize(ioms.size());
    for (int im=0; im<ioms.size(); im++) Sigh[fl][im]=0;
  }
  
  
  for (int im=0; im<ioms.size(); im++){
    dcomplex iomega(0,ioms[im]);
    //for (int fl=0; fl<N_unique_fl; fl++) Gf(fl,im) = Gh(fl, im);
    for (int ifl=0; ifl<N_ifl; ifl++)
      if (nonzero[ifl]) Inverse_Gf(ifl_dim[ifl], N_unique_fl, bfl_index[ifl], sign[ifl], conjg[ifl], im, Gh, Gf_1);
    for (int fl=0; fl<N_unique_fl; fl++) Gf_1(fl,im) /= fl_deg[fl];
    for (int fl=0; fl<N_unique_fl; fl++)
      if (nonzerou[fl]) Sigh[fl][im] = (iomega+mu)*Id[fl]-epsk[fl]-Gf_1(fl,im);
  }

  // brisi!
  /*
  Gf=0;
  Gf_1=0;
  ofstream hout("g_hb1.dat");
  ofstream sout("s_hb1.dat");
  ofstream rout("g_hb0.dat");
  for (int im=0; im<ioms.size(); im++){
    dcomplex iomega(0,ioms[im]);
    for (int fl=0; fl<N_unique_fl; fl++) Gf_1(fl,im) = (iomega+mu)*Id[fl]-epsk[fl]-Sigh[fl][im]-Delta[fl](iom.Interp(ioms[im]));
    for (int ifl=0; ifl<N_ifl; ifl++) Inverse_Gf(ifl_dim[ifl], N_unique_fl, bfl_index[ifl], sign[ifl], conjg[ifl], im, Gf_1, Gf);
    for (int fl=0; fl<N_unique_fl; fl++) Gf(fl,im) /= fl_deg[fl];
    hout<<ioms[im]<<" ";
    for (int fl=0; fl<N_unique_fl; fl++) hout<<Gf[fl][im]<<" ";
    hout<<endl;
    sout<<ioms[im]<<" ";
    for (int fl=0; fl<N_unique_fl; fl++) sout<<Sigh[fl][im]<<" ";
    sout<<endl;
    rout<<ioms[im]<<" ";
    for (int fl=0; fl<N_unique_fl; fl++) rout<<Gh[fl][im]<<" ";
    rout<<endl;
  }
  */
  //brisi!
  
  for (size_t fl=0; fl<Sigh.size(); fl++){
    int n = ioms.size()-1;
    dcomplex df0 = (Sigh[fl][1]-Sigh[fl][0])/(ioms[1]-ioms[0]);
    dcomplex dfn = (Sigh[fl][n]-Sigh[fl][n-1])/(ioms[n]-ioms[n-1]);
    Sigh[fl].splineIt(ioms, df0, dfn);
  }

  // Interpolates between low energy QMC and high energy HBI self-energy
  for (int ifl=0; ifl<N_unique_fl; ifl++){
    if (!nonzerou[ifl]) continue;
    dcomplex sb=0;
    for (int im=nom-aom; im<nom; im++) sb += Sigma[ifl][im];
    sb/=aom;
    double ob = iom[nom-(aom+1)/2];

    if (sb.imag()<0){// should be for physical system
      // trying to correct tail  s_tail -> (1/s_tail + cnst/wn)^{-1}
      double sim_tail = Sigh[ifl](ioms.Interp(ob)).imag();
      double cnst = (1/sb.imag()-1/sim_tail);
      for (int im=nom; im<iom.size(); im++){
	double s_tail = Sigh[ifl](ioms.Interp(iom[im])).imag();
	Sigma[ifl][im].imag() = s_tail/(1+s_tail*cnst);
      }
    }else{
      // Imaginary part
      int is=nom;
      for (; is<iom.size()-1; is++){
	intpar p0 = ioms.Interp(iom[is]);
	intpar p1 = ioms.Interp(iom[is-1]);
	double df0 = (Sigh[ifl](p0)-sb).imag()/(iom[is]-ob);
	double df1 = (Sigh[ifl](p0)-Sigh[ifl](p1)).imag()/(iom[is]-iom[is-1]);
	if (fabs(df0-df1)<sderiv) break;
      }
      double se = Sigh[ifl](ioms.Interp(iom[is])).imag();
      double oe = iom[is];

      //cout<<"nom="<<nom<<" Found that is="<<is<<" hence Sig(is)="<<se<<" with om(is)="<<oe<<" and Sig(nom)="<<sb<<" with nom="<<ob<<endl;

      for (int im=nom; im<is; im++){
	Sigma[ifl][im].imag() = sb.imag() + (se-sb.imag())*(iom[im]-ob)/(oe-ob);
	//cout<<"just set "<<iom[im]<<" "<<Sigma[ifl][im].imag()<<endl;
      }
      for (int im=is; im<iom.size(); im++){
	Sigma[ifl][im].imag() = Sigh[ifl](ioms.Interp(iom[im])).imag();
	//cout<<"just set "<<iom[im]<<" "<<Sigma[ifl][im].imag()<<endl;
      }
    }
    // Real part
    double sinf = Sigh[ifl].last().real();
    for (int im=nom; im<iom.size(); im++) Sigma[ifl][im].real() = sinf + sqr(ob/iom[im])*(sb.real()-sinf);
  }

  Gf.resize(N_unique_fl,iom.size());
  Gf_1.resize(N_unique_fl,iom.size());
  Gf=0;
  Gf_1=0;
  for (int im=0; im<iom.size(); im++){
    dcomplex iomega(0,iom[im]);
    for (int fl=0; fl<N_unique_fl; fl++)
      if (nonzerou[fl]) Gf_1(fl,im) = (iomega+mu)*Id[fl]-epsk[fl]-Sigma(fl,im)-Delta(fl,im);
    for (int ifl=0; ifl<N_ifl; ifl++)
      if (nonzero[ifl]) Inverse_Gf(ifl_dim[ifl], N_unique_fl, bfl_index[ifl], sign[ifl], conjg[ifl], im, Gf_1, Gf);
    for (int fl=0; fl<N_unique_fl; fl++) Gf(fl,im) /= fl_deg[fl];
  }
  
  ofstream tout("gs_qmc.dat");
  for (int io=0; io<iom.size(); io++){
    tout<<iom[io]<<" ";
    for (int ifl=0; ifl<N_unique_fl; ifl++) tout<<Gf[ifl][io]<<" "<<Sigma[ifl][io]<<" ";
    tout<<endl;
  }
}
