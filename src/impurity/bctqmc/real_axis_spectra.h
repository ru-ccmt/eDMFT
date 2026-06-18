#include "tinyTraceinv.h"

auto DetermineLambda0(const Array<double,1>& om, const Array<complex<double>,2>& Sm, const Cix& cix, int followPeak_fixed, double T, ostream& clog, double peak_pos=0){
  struct retVals{
    double lambda0;
    int followPeak;
  };
  retVals r;
  r.lambda0=100.; r.followPeak=-1;
  Array<double,1> work(3*cix.Nmaxsize); // for exact diagonalization
  double peak_posT = T*peak_pos;
  intpar ipr = Interp(peak_posT,om);
  int i0=ipr.i, N=om.size();
  int ib_start=0, ib_end=cix.Nblock; // normally we need to go through all sectors
  if (followPeak_fixed>0){           // when followPeak_fixed is set, we just stay in one sector
    ib_start=followPeak_fixed-1;
    ib_end = followPeak_fixed;
  }
  for (int ib=ib_start; ib<ib_end; ++ib){
    int N = cix.Block[ib].size(); // block size
    double lmbda=0;
    //
    //   lambda0 = cix.Ham[ib]+(Sm.real(0.0)-Sm.imag(0.0)*0.5*i) ~ cix.Ham[ib]+Sm.real(0.0)
    //
    if (N==1){ // no need to perform exact diagonalization, because the block is one dimensional
      double E = cix.Ham[ib](0,0);
      int ii = abs(cix.Sgm_symmetry[ib](0,0))-1; // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
      double S0 = real(Sm(ii,i0-2)+Sm(ii,i0-1)+Sm(ii,i0)+Sm(ii,i0+1)+Sm(ii,i0+2))/5.; // self-energy around zero
      lmbda = E + S0 - peak_posT; // strange, on imag axis it seems we use -S0, but probably we have different convention for sign of Sigma.
      if (fabs(S0)>1000){
	clog<<"ERROR S0="<<S0<<endl;
	ofstream out("1_debug.dat"); out.precision(12);
	out<<"# ii="<<ii<<" i="<<i0<<" S[i-2]="<<Sm(ii,i0-2)<<" S[i-1]="<<Sm(ii,i0-1)<<" S[i]="<<Sm(ii,i0)<<" S[i+1]="<<Sm(ii,i0+1)<<" S[i+2]="<<Sm(ii,i0+2)<<endl;
	for (int j=0; j<om.size(); ++j){
	  out<<om(j)<<" "<<Sm(ii,j).real()<<" "<<Sm(ii,j).imag()<<endl;
	}
	out.close();
	exit(1);
      }
      clog<<"lambda["<<ib+1<<"]="<<lmbda<<" Energy["<<ib+1<<"]="<<E<<" S0="<<S0<<endl;
    }else{
      Array<double,2> HS(N,N); HS=0;
      for (int i=0; i<N; ++i){
	for (int j=0; j<N; ++j){
	  int _ii_ = cix.Sgm_symmetry[ib](i,j);
	  int ii = abs(_ii_)-1;
	  if (ii>=0){
	    double S0 = real(Sm(ii,i0-2)+Sm(ii,i0-1)+Sm(ii,i0)+Sm(ii,i0+1)+Sm(ii,i0+2))/5.; // self-energy around zero
	    HS(i,j) = S0*sign(_ii_);
	  }
	}
      }
      HS += cix.Ham[ib];
      Array<double,1> E(N);  // stores eigenvalues
      xsyev(HS.extent(0), HS.data(), HS.extent(1), E.data(), work.data(), work.size(), "N");
      lmbda = min(E)-peak_posT;
      clog<<"lambda["<<ib+1<<"]="<<lmbda<<" E+S0["<<ib+1<<"]="<<E<<endl;
    }
    if (lmbda<r.lambda0){
      r.lambda0=lmbda;
      r.followPeak=ib+1;
    }
  }
  clog<<"lmbda0="<<r.lambda0<<" followPeak="<<r.followPeak<<endl;
  return r;
}
auto DeterminSpectralFunctions(Array<complex<double>,2>& Gm, Array<complex<double>,2>& Sm, Array<double,1>& om, 
			       const Cix& cix, double T, bool AdaptMesh, json& p, MPIcls& mpi, double peak_pos=0,
			       double lambda0_fixed=0, int followPeak_fixed=-1)
{
  // looking for lambda such that \widetilde{G} has maximum at zero frequency.
  // Sufficient condition is that the derivative of 1/\widetilde{G} is zero at zero frequency.
  // One gets a quadratic equation for lambda and thus two roots. Then one chooses the root that maximizes \widetilde{G}.
  // If no root exists, than we take lambda that minimizes linear coeficient in the expansion of 1/\widetilde{G}.
  // The latter equation is linear and one always gets unique solution.
  struct retVals{
    double lambda0;
    int followPeak;
    Array<double,1> few;
  };
  clog.precision(14);
  bool debug = false;
  retVals r;
  ostream& clog = *mpi.out;
  if (lambda0_fixed==0){
    auto [lambda0,followPeak] = DetermineLambda0(om,Sm,cix,followPeak_fixed,T,clog,peak_pos);
    r.lambda0 = lambda0;
    r.followPeak = followPeak;
  }else{
    r.lambda0 = lambda0_fixed;
    clog<<"using previously-fixed lambda0="<<lambda0_fixed<<endl;
    r.followPeak = 0;
  }
  double small=p["Nw_small"], alpha=p["Nw_alpha"], La=p["Nw_La"], Lb=p["Nw_Lb"];
  int Nw_1=p["Nw_1"], Nw_2=p["Nw_2"], Nw_max=p["Nw_max"], Nw_itt=p["Nw_itt"];
  //bool fastFilesystem=false;//p["fastFilesystem"];

  //cout<<"small="<<small<<" alpha="<<alpha<<" La="<<La<<" Lb="<<Lb<<" Nw_1="<<Nw_1<<" Nw_2="<<Nw_2<<" Nw_max="<<Nw_max<<" Nw_itt="<<Nw_itt<<endl;
  r.few.resize(om.size());
  for (int i=0; i<r.few.size(); ++i) r.few(i)=fermi(om(i)/T);
  Array<double,1> dh = GiveDh(om);
  /*
  if (mpi.rank == mpi.master) {
    std::fprintf(stderr, "[root] AdaptMesh=%d\n", AdaptMesh);
  } else {
    std::fprintf(stderr, "[rank %d pre] AdaptMesh=%d\n", mpi.rank, AdaptMesh);
  }
  fflush(stderr);
  MPI_Barrier(MPI_COMM_WORLD);
  */
  if (AdaptMesh){
    if (mpi.rank==mpi.master){
      Array<double,2> App(cix.Nblock,om.size()); App=0;
      // mesh is constructed with the use of fw0 = \sum_p G_p(om)*dom/Integrate(G(om)*dom)
      // so that new points are such that each interval contributes the same weight
      for (int ib=0; ib<cix.Nblock; ++ib){
	int N = cix.Block[ib].size(); // block size
	if (N==1){ // no need to perform exact diagonalization, because the block is one dimensional
	  int ii = abs(cix.Sgm_symmetry[ib](0,0))-1; // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
	  double mune = -cix.Ham[ib](0,0) + r.lambda0;
	  for (int iw=0; iw<om.size(); iw++){
	    double St = Sm(ii,iw).imag();
	    complex<double> Sc(Sm(ii,iw).real(), -St*(1-r.few(iw)));
	    complex<double> Gc = 1.0/(om(iw)+mune-Sc);
	    double Gt = St*real(Gc*conj(Gc));
	    App(ib,iw) = Gt/M_PI;
	    Gm(ii,iw) = complex<double>(Gc.real(), Gt);	    
	  }
	}else{
	  //  Ham = H^+ * E * H
	  //  E = Ht^+ * Ham * Ht
	  // Sdiag(in) = Ht(i,in)*Sm(i,j)*Ht(j,in)
	  Eigen::MatrixXcd Gci(N,N), Gc(N,N);
	  Eigen::MatrixXd St(N,N), Gt(N,N);
	  for (int iw=0; iw<om.size(); ++iw){
	    Gci.setZero();
	    St.setZero();
	    Gci += (om(iw)+r.lambda0) * Eigen::MatrixXcd::Identity(N,N);
	    for (int i=0; i<N; ++i){
	      for (int j=0; j<N; ++j){
		Gci(i,j) -= cix.Ham[ib](i,j);
		int _ii_ = cix.Sgm_symmetry[ib](i,j);
		int ii = abs(_ii_)-1;
		if (ii>=0){
		  double ss = sign(_ii_); //// HERE WAS BUG Jan 6. 2026
		  int iit = abs(cix.Sgm_symmetry[ib](j,i))-1;
		  complex<double> Sc = Sm(ii,iw), Sct = Sm(iit,iw); // Sigma and Sigma^T, with proper sign to all terms, including St!
		  // The bug was in sign of St, which did not have ss in the evaluation. Scx amd Gci were correct.
		  St(i,j) = ss*imag(0.5*(Sc-conj(Sct))); // in matrix form we want to have St = (Sig-Sig^+)/(2*i)
		  complex<double> Scx(Sc.real(), -Sc.imag()*(1-r.few(iw)));
		  Gci(i,j) -= ss*Scx;
		}
	      }
	    }
	    //cout<<"ib="<<ib<<" iw="<<iw<<" Gci="<<Gci<<endl;
	    Gc = Gci.inverse();
	    // Whether Gt is affected by sign change of a few coponents in St is an interesting question.
	    Gt = (Gc * St * Gc.adjoint()).real();
	    for (int i=0; i<N; ++i) App(ib,iw) += Gt(i,i)/M_PI;
	    // This is not needed, but for comparison we write into Gm what we currently have
	    for (int i=0; i<N; ++i){
	      for (int j=0; j<N; ++j){
		int _ii_ = cix.Sgm_symmetry[ib](i,j); // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
		if (_ii_>0 and bl::all(cix.Sgm_indx[_ii_-1][0] == TinyVector<int,3>(ib,i,j))){
		  // We write each entry only once by the first occurence, which is cix.Sgm_indx[ii][0].
		  // we also can safely assume that the first occurence has positive ii index and positive ib index
		  int ii = _ii_-1;
		  //if (DEBUG) clog<<"setting Gm["<<ii<<"]+= G("<<i<<","<<j<<") with deg="<<deg[ii]<<endl;
		  Gm(ii,iw) = complex<double>(Gc(i,j).real(), Gt(i,j));
		}
	      }
	    }
	  }
	}
      }
      int Nps = App.extent(0);
      Array<double,1> integral0(Nps);
      for (int i=0; i<Nps; ++i){
	integral0(i) = sum(App(i,Range::all())*(1-r.few(Range::all()))*dh(Range::all())); // Integrate[Tr(Im(Gc))]
	//integral0(i) = sum(App(i,Range::all())*dh(Range::all()));
      }
      
      // mesh is constructed with the use of fw0 = \sum_p G_p(om)*dom/Integrate(G(om)*dom)
      // so that new points are such that each interval contributes the same weight
      bl::Array<double,1> fw2 = BuildFw2FromApp(App, dh, alpha, 1e-10);

      if (debug){
	ofstream dbg1("cdebug.1"); dbg1.precision(14);
	for (int i=0; i<om.size(); ++i){
	  // fw2 = fw1^alpha;  fw1 = smooth(fw0); fw0(w) = \sum_i Gt(w)(i,i)*dh(w)/Integrate[Gt(w)(i,i)]
	  dbg1<<om(i)<<" "<<fw2(i)<<" ";//<<fw0(i)<<" "<<fw1(i)<<" ";
	  for (int ip=0; ip<Nps; ip++) dbg1<<App(ip,i)<<" ";
	  dbg1<<endl;
	}
      }
      // Build an adaptive frequency mesh xn from a global weight fw2(ω) on an existing grid om.
      Array<double,1> xn = BuildAdaptiveMeshFromFw2(om, fw2, App, La, Lb, Nw_1, Nw_2, small, Nw_max, Nw_itt, clog, debug);
      
      Array<complex<double>,2> Sm_copy(Sm.shape()); Sm_copy = Sm;

      if (debug){
	ofstream dbg2("cdebug.4"); dbg2.precision(14);
	for (int i=0; i<xn.size(); ++i)
	  dbg2<<i<<" "<<xn(i)<<endl;
	// Printing old result
	ofstream outo("Speco.dat"); outo.precision(15);
	for (int i=0; i<om.size(); ++i){
	  outo << om(i)<<"  ";
	  //for (int idim=0; idim<Nps; ++idim) outo<<App(idim,i)<<" ";
	  for (int ii=0; ii<cix.Npseudo_unique; ++ii) outo<<Gm(ii,i).imag()/M_PI<<" ";
	  for (int ii=0; ii<cix.Npseudo_unique; ++ii) outo<<Gm(ii,i).real()<<" ";
	  outo<<endl;
	}
	ofstream outp("Sigo.dat"); outp.precision(15);
	for (int i=0; i<om.size(); ++i){
	  outp << om(i)<<"  ";
	  for (int ii=0; ii<cix.Npseudo_unique; ++ii) outp<<Sm_copy(ii,i).real()<<" "<<Sm_copy(ii,i).imag()<<" ";
	  outp<<endl;
	}
      }
      
      // Now recomputing fx/fn on new mesh xn
      //
      Sm.resize(cix.Npseudo_unique,xn.size());
      Gm.resize(cix.Npseudo_unique,xn.size());
      r.few.resize(xn.size());
      for (int i=0; i<r.few.size(); ++i) r.few(i)=fermi(xn(i)/T);
      //
      for (int ib=0; ib<cix.Nblock; ++ib){
	int N = cix.Block[ib].size(); // block size
	if (N==1){
	  int pos=0;
	  for (int iw=0; iw<xn.size(); ++iw){
	    intpar ip=InterpLeft(xn(iw),pos,om);
	    int ii = abs(cix.Sgm_symmetry[ib](0,0))-1;
	    if (ii>=0){
	      complex<double> Sc  = linInterp(Sm_copy,ii,ip); // Sigma(x(om))
	      if (Sc.imag()<0) Sc = complex<double>(Sc.real(),1e-10);
	      complex<double> Scx(Sc.real(), -Sc.imag()*(1-r.few(iw)));
	      complex<double> Gc = 1.0/(xn(iw) + r.lambda0 - cix.Ham[ib](0,0) - Scx);
	      double Gt = Sc.imag()*real(Gc*conj(Gc));
	      if (Gt<0) Gt=1e-100;
	      Gm(ii,iw) = complex<double>(Gc.real(), Gt);
	      Sm(ii,iw) = Sc;
	    }
	  }
	}else{
	  Eigen::MatrixXcd Gci(N,N), Gc(N,N);
	  Eigen::MatrixXd St(N,N), Gt(N,N);
	  int pos=0;
	  for (int iw=0; iw<xn.size(); ++iw){
	    intpar ip=InterpLeft(xn(iw),pos,om);
	    Gci.setZero();
	    St.setZero();
	    Gci += (xn(iw)+r.lambda0) * Eigen::MatrixXcd::Identity(N,N);
	    for (int i=0; i<N; ++i){
	      for (int j=0; j<N; ++j){
		Gci(i,j) -= cix.Ham[ib](i,j);
		int _ii_ = cix.Sgm_symmetry[ib](i,j);
		int ii = abs(_ii_)-1;
		if (ii>=0){
		  double ss = sign(_ii_); //// HERE WAS BUG Jan 6. 2026
		  int iit = abs(cix.Sgm_symmetry[ib](j,i))-1;
		  complex<double> Sc  = linInterp(Sm_copy,ii,ip); // Sigma(x(om))
		  //if (i==j and Sc.imag()<0) Sc = complex<double>(Sc.real(),1e-10); // for causality
		  Sm(ii,iw) = Sc;
		  complex<double> Sct = Sc;
		  if (iit!=ii){
		    Sct = linInterp(Sm_copy,iit,ip)*ss; // Sigma^T(x(om))
		    //if (i==j and Sct.imag()<0) Sct = complex<double>(Sct.real(),1e-10);
		  }
		  St(i,j) = ss*imag(0.5*(Sc-conj(Sct))); //  HERE WA A BUG  in matrix form we want to have St = (Sig-Sig^+)/(2*i)
		  complex<double> Scx(Sc.real(), -Sc.imag()*(1-r.few(iw)));
		  Gci(i,j) -= ss*Scx;
		}
	      }
	    }
	    Gc = Gci.inverse();
	    Gt = (Gc * St * Gc.adjoint()).real();
	    // Now we write it back
	    for (int i=0; i<N; ++i){
	      for (int j=0; j<N; ++j){
		int _ii_ = cix.Sgm_symmetry[ib](i,j); // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
		if (_ii_>0 and bl::all(cix.Sgm_indx[_ii_-1][0] == TinyVector<int,3>(ib,i,j))){
		  // We write each entry only once by the first occurence, which is cix.Sgm_indx[ii][0].
		  // we also can safely assume that the first occurence has positive ii index and positive ib index
		  int ii = _ii_-1;
		  //if (DEBUG) clog<<"setting Gm["<<ii<<"]+= G("<<i<<","<<j<<") with deg="<<deg[ii]<<endl;
		  Gm(ii,iw) = complex<double>(Gc(i,j).real(), Gt(i,j));
		}
	      }
	    }
	  }
	}
      }
      om.resize(xn.size());
      om = xn;

      // Check quality of the integral
      dh.resize(xn.size());
      dh = GiveDh(xn);
      Array<complex<double>,1> gdh(xn.size());
      for (int ib=0; ib<cix.Nblock; ++ib){
	int N = cix.Block[ib].size(); // block size
	double intg=0;
	for (int i=0; i<N; ++i){
	  int ii = abs(cix.Sgm_symmetry[ib](i,i))-1; // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
	  if (ii>=0){
	    //gdh = Gm(ii,Range::all())*dh(Range::all());
	    gdh = Gm(ii,Range::all())*(1-r.few(Range::all()))*dh(Range::all());
	    intg  += sum(gdh).imag()/M_PI;
	  }
	}
	clog<<"# Final integral["<<ib+1<<"]="<<intg/N <<" and original "<<integral0(ib)/N << " with diff="<<(intg-integral0(ib))/N<<endl;
      }
      // Printing final results
      if (debug){
	ofstream outs("Sign.dat"); outs.precision(15);
	for (int i=0; i<xn.size(); ++i){
	  outs << xn(i)<<"  ";
	  for (int ii=0; ii<cix.Npseudo_unique; ++ii) outs<<Sm(ii,i).real()<<" "<<Sm(ii,i).imag()<<" ";
	  outs<<endl;
	}
	ofstream outn("Specn.dat"); outn.precision(15);
	for (int i=0; i<xn.size(); ++i){
	  outn << xn(i)<<"  ";
	  for (int ii=0; ii<cix.Npseudo_unique; ++ii) outn<<Gm(ii,i).imag()/M_PI<<" ";
	  for (int ii=0; ii<cix.Npseudo_unique; ++ii) outn<<Gm(ii,i).real()<<" ";
	  outn<<endl;
	}
      }
    }
    mpi.Bcast_New_mesh(om,r.few,Sm,Gm);
  }else{
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size(); // block size
      if (N==1){ // no need to perform exact diagonalization, because the block is one dimensional
	int ii = abs(cix.Sgm_symmetry[ib](0,0))-1; // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
	double mune = -cix.Ham[ib](0,0) + r.lambda0;
	for (int iw=0; iw<om.size(); iw++){
	  double St = Sm(ii,iw).imag();
	  complex<double> Sc(Sm(ii,iw).real(), -St*(1-r.few(iw)));
	  complex<double> Gc = 1.0/(om(iw)+mune-Sc);
	  double Gt = St*real(Gc*conj(Gc));
	  Gm(ii,iw) = complex<double>(Gc.real(), Gt);	    
	}
      }else{
	Eigen::MatrixXcd Gci(N,N), Gc(N,N);
	Eigen::MatrixXd St(N,N), Gt(N,N);
	for (int iw=0; iw<om.size(); ++iw){
	  Gci.setZero();
	  St.setZero();
	  Gci += (om(iw)+r.lambda0) * Eigen::MatrixXcd::Identity(N,N);
	  for (int i=0; i<N; ++i){
	    for (int j=0; j<N; ++j){
	      Gci(i,j) -= cix.Ham[ib](i,j);
	      int _ii_ = cix.Sgm_symmetry[ib](i,j);
	      int ii = abs(_ii_)-1;
	      if (ii>=0){
		double ss = sign(_ii_); // HERE WAS BUG Jan 6. 2026
		int iit = abs(cix.Sgm_symmetry[ib](j,i))-1;
		complex<double> Sc = Sm(ii,iw), Sct = Sm(iit,iw); // Sigma and Sigma^T
		St(i,j) = ss*imag(0.5*(Sc-conj(Sct))); // in matrix form we want to have St = (Sig-Sig^+)/(2*i)
		complex<double> Scx(Sc.real(), -Sc.imag()*(1-r.few(iw)));
		Gci(i,j) -= ss*Scx;
	      }
	    }
	  }
	  Gc = Gci.inverse();
	  Gt = (Gc * St * Gc.adjoint()).real();
	  for (int i=0; i<N; ++i){
	    for (int j=0; j<N; ++j){
	      int _ii_ = cix.Sgm_symmetry[ib](i,j); // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
	      if (_ii_>0 and bl::all(cix.Sgm_indx[_ii_-1][0] == TinyVector<int,3>(ib,i,j))){
		// We write each entry only once by the first occurence, which is cix.Sgm_indx[ii][0].
		// we also can safely assume that the first occurence has positive ii index and positive ib index
		int ii = _ii_-1;
		Gm(ii,iw) = complex<double>(Gc(i,j).real(), Gt(i,j));
	      }
	    }
	  }
	  /*
          if (iw==100 and ib==40){
	    cout << "lmbda0="<<r.lambda0<< endl;
	    cout << "Ham=" << cix.Ham[ib] << endl;
            cout << "St["<<ib+1<<"]="<< St <<endl;
            cout << "Gc["<<ib+1<<"]="<< Gc <<endl;
            cout << "Gt["<<ib+1<<"]="<< Gt <<endl;
          }
	  */
	}
      }
    }
  }
  
  return r;
}

auto DeterminProbabilities(const Array<double,1>& om, const Array<complex<double>,2>& Gm, const Cix& cix,
			   const Array<double,1>& few, double T, double mu, double lambda0, int followPeak, int itt, MPIcls& mpi, double small=5e-7)
{
  struct retVals{
    double nd, Q, lambda0, Fimp, Epot, TrSigmaG;
    int followPeak;
    Array<double,1> density_matrix, Prob;
  };
  retVals r;
  r.followPeak = followPeak;
  ostream& clog = *mpi.out;
  Array<double,1> Probability(cix.Natst), nalpha(cix.Nbath);
  Probability=0; nalpha = 0;

  Array<double,1> dh = GiveDh(om);
  
  double sumQ = 0, sum_nd=0, dQ=0;
  Array<double,1> dQm(cix.Npseudo_unique);
  for (int ii=0; ii<cix.Npseudo_unique; ++ii){
    double dQ=0;
    for (int iw=0; iw<om.size(); iw++) dQ += Gm(ii,iw).imag()*few(iw)*dh(iw);
    dQm(ii) = dQ/M_PI;
  }
  r.Epot=0;
  for (int ib=0; ib<cix.Nblock; ++ib){
    int N = cix.Block[ib].size();
    double sumQ1=0;
    for (int i=0; i<N; ++i){
      int i_pseudo = cix.Block[ib](i);
      int ii = cix.Sgm_symmetry[ib](i,i)-1; // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
      // diagonal components should not have negative ii. This is only allowed for off-diagonal
      if (ii>=0){
	double Nocc = cix.iN(i_pseudo-1);
	dQ = dQm(ii); 
	sumQ += dQ;
	sumQ1 += dQ;
	Probability(i_pseudo-1) += dQ;
	sum_nd += dQ*Nocc;
      }
      for (int ibth=0; ibth<cix.Nbath; ++ibth){
	int j_pseudo = cix.Fi(2*ibth+1,i_pseudo);      // <j_pseudo|F_{ibth}|i_pseudo>
	if (j_pseudo!=0){
	  int k_pseudo = cix.Fi(2*ibth,abs(j_pseudo)); // <k_pseudo|F^+_{ibth}|j_pseudo><j_pseudo|F_{ibth}|i_pseudo>
	  if (k_pseudo!=0){
	    int k = cix.inBlock(abs(k_pseudo));
	    int _ii_ = cix.Sgm_symmetry[ib](k,i);
	    int ii = abs(_ii_)-1;
	    if (ii>=0){
	      nalpha(ibth) += sign(_ii_*j_pseudo*k_pseudo)*dQm(ii); // dQ[k_pseudo,i_pseudo] * <k_pseudo|F^+_{ibth}|j_pseudo><j_pseudo|F_{ibth}|i_pseudo>
	      //cout<<"ibth="<<ibth<<" i_pseudo="<<i_pseudo<<" j_pseudo="<<j_pseudo<<" k_pseudo="<<k_pseudo<<" ii="<<ii<<" (i,k)=("<<i<<","<<k<<")"<<" dnalpha="<<sign(_ii_*j_pseudo*k_pseudo)*dQm(ii)<<endl;
	    }
	  }
	}
      }
      for (int j=0; j<N; ++j){
	int _ii_ = cix.Sgm_symmetry[ib](i,j);
	int ii = abs(_ii_)-1;
	if (ii>=0) r.Epot += sign(_ii_)*dQm(ii)*cix.Ham[ib](j,i);
      }
    }
  }
  r.Prob.resize(cix.Npseudo_unique);
  for (int ii=0; ii<cix.Npseudo_unique; ++ii) r.Prob(ii) = dQm(ii)/sumQ;
  
  Probability *= 1/sumQ;
  r.nd = sum_nd/sumQ;
  nalpha *= 1/sumQ;
  r.Q = sumQ;
  r.lambda0 = lambda0;
  r.Epot *= 1/sumQ;
  r.Epot += r.nd*mu; // this is because r.Epot = Tr(H^{atom}_m P_m), but H^{{atom} was previously modified to H^{atom}-mu*N
  r.Fimp = lambda0-log(r.Q)*T;
  
  double Epot_singlep=0;
  for (int b=0; b<cix.Nbath; b++){
    int iib = cix.bath_sigind0(b);
    int ib_sg = 1;
    if (iib<0){ iib=-iib; ib_sg=-1;}
    Epot_singlep += (cix.Ebath(iib)+mu)*ib_sg*nalpha(b);  // previously we subtracted Ebath-=mu, now we add back, like in Epot
  }
  r.TrSigmaG = r.Epot - Epot_singlep;

  r.density_matrix.resize(cix.Nbath_unique_unique); r.density_matrix=0;
  Array<int,1> deg(cix.Nbath_unique_unique); deg=0;
  for (int iorb=0; iorb<cix.Nbath; ++iorb){
    int iib0 = cix.bath_sigind0(iorb);
    double ib_sg = (iib0 < 0 ? -1 : 1);
    int iib   = std::abs(iib0);
    r.density_matrix(iib) += nalpha(iorb)*ib_sg;
    deg(iib) += 1;
  }
  r.density_matrix /= deg;

  if (mpi.rank==mpi.master){
    clog<<"   "<<itt<<") nd= "<<r.nd<<endl;
    clog<<"   "<<itt<<") Q = "<<r.Q<<endl;
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size();
      double Ptot=0;
      for (int i=0; i<N; ++i) Ptot += Probability(cix.Block[ib](i)-1);
      if (Ptot>small){
	int nmax=N;
	for (; nmax>0; --nmax){
	  if (abs(Probability(cix.Block[ib](nmax-1)-1)) > 0.2*small) break;
	}
	int Nocc = cix.iN(cix.Block[ib](0)-1);
	//clog<<"   "<<itt<<") Prob["<<ib+1<<",N="<<Nocc<<"]=[";
	fmt::print(clog,"   {:d}) Prob[{:d},N={:d}]={:10.7f} [",itt,ib+1,Nocc,Ptot);
	for (int i=0; i<nmax; ++i)
	  fmt::print(clog,"{:10.7f},",Probability(cix.Block[ib](i)-1));
	clog<<"]"<<endl;
      }
    }
    clog<<"   "<<itt<<") Fimp = "<<r.Fimp<<endl;
    clog<<"   "<<itt<<") Epot = "<<r.Epot<<endl;
    clog<<"   "<<itt<<") TrSigmaG = "<<r.TrSigmaG<<endl;
    clog<<"   "<<itt<<") Epot_singlep = "<<Epot_singlep<<endl;
    //clog<<"   occupancies: "<<endl<<"   ";
    clog.setf(ios::fixed);
    clog<<"   "<<itt<<") bath_density_matrix="<<r.density_matrix<<endl; //nalpha<<endl;
    clog.unsetf(ios::fixed);
  }
  return r;
}

Array<double,2> DeterminPhiR(const Array<complex<double>,2>& Gm, const Array<double,1>& xn, const Array<double,1>& few, const Array<double,1>& dh, const Array<double,1>& logo, double beta)
{
  // This quantity is needed to evaluate Phi on the real axis. We want Phi = 1/Norder Tr(\Sigma_{pseudo} G_{pseudo}) and in time we just need to evaluate time integral.
  // In frequency is harder. We first show that
  // Tr(Sigma_p * G_p) = -e^{-beta*lambda} \int dx/pi e^{-beta*x} Im(\Sigma_p(x) G_p(x)) = -e^{-beta*lambda} \int dx/pi e^{-beta*x} (\Sigma'_p(x) G''_p(x))+\Sigma''_p(x) G'_p(x))
  // Taking into account the Kramars-Kronig property \Sigma'_p(x) = -1/pi \int \Sigma''_p(y)/(x-y)dy (notice no static part to Sigma')
  // we can show that
  // Tr(Sigma_p * G_p) = -e^{-beta*lambda} \int dx/pi [ e^{-beta*x} \Sigma''_p(x) G'_p(x) - \Sigma''_p(x) \widetilde{G}'_p(x) ]
  //
  // where \widetilde{G}'_p(x) = -1/pi \int e^{-beta y} G''_p(y)/(x-y)dy=-1/pi \int f(y) \widetilde{G}_p(y)/(x-y)dy
  // hence we can compute 
  // Tr(Sigma_p * G_p) = -e^{-beta*lambda} \int dx/pi [f(x)G'_p(x)-f(-x)\widetilde{G}'_p(x)]*widetilde{\Sigma}_p(x)
  // hence we need PhiR_p(x) == -[f(x)G'_p(x)-f(-x)\widetilde{G}'_p(x)]/pi and Phi= 1/Q*\int dx PhiR_p(x)*widetilde{\Sigma}_p(x)
  //
  // Alternative derivation:
  //   Phi = -1/<Q> * \int( PhiR_p(w) * \widetilde{\Sigma}_p(w), {w,-\inf,\inf})
  //   This comes from the fact that
  //   Phi = -1/<Q> * Int[ e^{-\beta x} * Im(G_p(x)\Sigma_p(x)) ,{x,-\inf,\inf}]/pi
  //   which we transform to
  //   PhiR(w) = -( Re(G_p(w))*fe(w) + fe(-w)*Int[ widetilde{G}_p(x)*f(x)/(w-x)*(-1/pi),{x,-inf,inf}] )/pi
  //   PhiR(w) = -( Re(G_p(w))*fe(w) + fe(-w)*KK_w(f(y)\widetilde{G}_p(y)) )/pi
  //
  // Note that this element-wise propery works in the absence of broken symmetry. In general we might need to work with matrices also here.
  int  Npseudo_unique=Gm.extent(0);
  Array<double,2> PhiR(Npseudo_unique,xn.size());
  Array<double,1> mGm(xn.size()), mGr(xn.size());
  Array<double,1> iGm(xn.size());
  Array<double,1> pGr(xn.size());
  for (int ip=0; ip<Npseudo_unique; ++ip){
    iGm(Range::all()) = imag(Gm(ip,Range::all())); // this is -Im(G)/f(-om) = -Gt(om)
    pGr(Range::all()) = real(Gm(ip,Range::all())); // this is  Re(G)
    mGm = iGm*few;                                 // this is -Gt(om)*f(om)
    mGr = KramarsKronig(mGm, xn, dh, logo);        // this is -Re(Gt) == -KK( Gt(om)*f(om))
    PhiR(ip,Range::all()) = -(pGr(Range::all())*few(Range::all()) + mGr(Range::all())*(1.0-few(Range::all())))/M_PI;
  }
  /*
  ofstream out("PhiR.debug"); out.precision(12);
  for (int iw=0; iw<xn.size(); ++iw){
    out<<xn(iw)<<" ";
    for (int ip=0; ip<Npseudo_unique; ++ip){
      //double ph1 = -Gm(ip,iw).real()*few(iw)/M_PI;
      //double ph2 = -mGr_temp(ip,iw)*(1-few(iw))/M_PI;
      //out<<PhiR(ip,iw)<<" "<<ph1<<" "<<ph2<<" ";
      out<<PhiR(ip,iw)<<" ";
    }
    out<<endl;
  }
  */
  /*
  Array<complex<double>,1> dSm(xn.size());
  int N = xn.size();
  double dF=0;
  for (int ip=0; ip<Npseudo_unique; ++ip){
    dSm(0) = (Sm(1)-Sm(0))/(xn(1)-xn(0));
    for (int iw=1; iw<xn.size()-1; ++iw) dSm(iw) = 0.5*(Sm(iw+1)-Sm(iw))/(xn(iw+1)-xn(iw)) + 0.5*(Sm(iw)-Sm(iw-1))/(xn(iw)-xn(iw-1));
    dSm(N-1) = (Sm(N-1)-Sm(N-2))/(xn(N-1)-xn(N-2));
    double df=0;
    for (int iw=0; iw<xn.size(); ++iw){
      double cz = Gm(ip,iw).imag()*(1-dSm(iw).real())-Gm(ip,iw).real()*(dSm(iw).imag()+Sm(iw).imag()*few(iw));
      df += cz*few(iw)*dh(iw)/M_PI;
    }
    dF += df;
  }
  double F0 = log(dF)/beta;
  cout<<"F_0="<<F0<<endl;
  */
  return PhiR;
}


      
double DifferenceS(const Array<double,1>& xn, const Array<double,2>& Smt_new, const Array<double,1>& xno, const Array<double,2>& Smt)
{
  if (xn.size()==xno.size() and sum(abs(xn-xno))<1e-5){
    return sum(abs(Smt_new-Smt))/Smt_new.size();
  }else{
    double diff=0;
    int pos=0;
    for (int i=0; i<xn.size(); ++i){
      intpar ipx = InterpLeft(xn(i), pos, xno);
      for (int ip=0; ip<Smt.extent(0); ++ip){
	double Smt_old = linInterp(Smt, ip, ipx);
	diff += fabs(Smt_new(ip,i)-Smt_old);
      }
    }
    diff /= Smt_new.size();
    return diff;
  }
}

#include "Eigen/Dense"
using Eigen::MatrixXcd;  // complex<double> matrix

void ComputePhysicalSelfEnergyOnMesh(Array<complex<double>,2>& Sigm, Array<complex<double>,2>& SigD, Array<complex<double>,2>& Gc,
				     const Array<double,1>& om, const Array<complex<double>,2>& Gimp, const Array<complex<double>,2>& Fimp1, const Array<complex<double>,2>& Fimp2,
				     const Array<complex<double>,2>& iDeltaw,
				     const BathBlocks& bb, const Cix& cix, ostream& clog/**mpi.out*/, double small, bool debugPrint=false)
{

  vector<bool> component_set(cix.Nbath_unique_unique,false);
  for (int ib=0; ib<(int)bb.blocks.size(); ++ib){
    int N = bb.blocks[ib].size();
    if (N==1){ // no matrix needed
      for (int iorb=0; iorb<cix.Nbath; ++iorb){
	if (bb.orb_ind[iorb][0]==ib){
	  int iib0 = cix.bath_sigind0(iorb);
	  double ib_sg = (iib0 < 0 ? -1 : 1);
	  int iib   = std::abs(iib0);
	  if (not component_set[iib]){
	    double Eimp = cix.Ebath(iib);
	    Sigm(iib,Range::all())  = Fimp1(iib,Range::all())/Gimp(iib,Range::all());
	    SigD(iib,Range::all()) = -1.0/Gimp(iib,Range::all())*ib_sg-iDeltaw(iib,Range::all())*ib_sg-Eimp*ib_sg+om(Range::all());
	    Gc  (iib,Range::all()) =  1.0/(om(Range::all())-Eimp*ib_sg-Sigm(iib,Range::all())*ib_sg-iDeltaw(iib,Range::all())*ib_sg);
	    component_set[iib]=true;
	    //
	    for (int iw=0; iw<om.size(); ++iw){
	      if (Sigm(iib,iw).imag()>0){ // causality violation
		//std::cout << "WARNING Sg["<<iorb<<",om="<<om(iw)<<"] is not causal "<< Sigm(iib,iw).imag() << std::endl;
		Sigm(iib,iw) = std::complex<double>(Sigm(iib,iw).real(),-small);
	      }
	    }
	  }
	}
      }
    }else{
      bool SC = bb.is_block_sc[ib];
      MatrixXcd Gm(N, N), Fm(N, N), Dm(N, N), om_Eimp(N, N); // Eigen matrices for linear algebra
      ofstream out1, out2, out3, out4, out5, out6;
      if (debugPrint){
	out1.open("Sm.dat");
	out2.open("Gmc.dat");
	out3.open("Dm.dat");
	out4.open("Gm.dat");
	out5.open("Fm.dat");
	out6.open("om_Eimp.dat");
      }
      for (int iw=0; iw<om.size(); ++iw){
	intpar ip;
	if (SC) ip = Interp(-om(iw), om);
	Gm.setZero();
	Fm.setZero();
	Dm.setZero();
	om_Eimp.setZero();
	// Diagonal: ω * I
	for (int i=0; i<N; ++i) om_Eimp(i,i) = om(iw);
	// Fill Gm, Fm, Dm, om_Eimp from bath data
	for (int iorb=0; iorb<cix.Nbath; ++iorb){
	  if (bb.orb_ind[iorb][0]==ib){
	    int j=bb.orb_ind[iorb][1], k=bb.orb_ind[iorb][2];
	    int iib0 = cix.bath_sigind0(iorb);
	    double ib_sg = (iib0 < 0 ? -1 : 1);
	    int iib   = std::abs(iib0);
	    if (SC && j>=N/2 && k>=N/2){
	      if (iib0>=0){
		// SC/Nambu “conjugated” sector, namely, G_{22}(om+i*delta) = -G^H_{11}(-om+i*delta)
		Gm(k,j) = -std::conj(linInterp(Gimp, iib, ip));
		Dm(k,j) = -std::conj(linInterp(iDeltaw, iib, ip));
		om_Eimp(k,j) += cix.Ebath(iib);
	      }else{
		// SC/Nambu “conjugated” sector, namely, G_{22}(om+i*delta) = -G^H_{11}(-om+i*delta)
		Gm(k,j) = std::conj(linInterp(Gimp, iib, ip));
		Dm(k,j) = std::conj(linInterp(iDeltaw, iib, ip));
		om_Eimp(k,j) -= cix.Ebath(iib);
	      }
	    }else{
	      if (iib0>=0){
		// “normal” state or “normal” sector
		Gm(j,k) = Gimp(iib,iw);
		Dm(j,k) = iDeltaw(iib,iw);
		om_Eimp(j,k) -= cix.Ebath(iib);
	      }else{
		// “normal” state or “normal” sector
		Gm(j,k) = -Gimp(iib,iw);
		Dm(j,k) = -iDeltaw(iib,iw);
		om_Eimp(j,k) += cix.Ebath(iib);
	      }
	    }
	    if (SC){
	      if (k>=N/2){
		if (j>=N/2) // F_{22} has an extra minus sign, and takes down, which is stored in F2
		  Fm(k,j) = ib_sg*std::conj(linInterp(Fimp2, iib, ip));
		else        // F_{12} has also an extra minus sign, and it destroys cooper pair, hence it is stored in F2
		  Fm(j,k) = -ib_sg*Fimp2(iib,iw);
	      }else{
		Fm(j,k) = ib_sg*Fimp1(iib,iw);
	      }
	    }else{
	      //Fm(j,k) = ib_sg*0.5*(Fimp1(iib,iw)+Fimp2(iib,iw));
	      Fm(j,k) = ib_sg*Fimp1(iib,iw);
	    }
	  }
	}
	MatrixXcd Gminv = Gm.inverse();          // Gm^{-1}
	MatrixXcd Sm    = Gminv * Fm;            // Sm = Gm^{-1} * Fm
	MatrixXcd SmD   = om_Eimp - Gminv - Dm;  // Dyson: Σ = (ω - E_imp) - G^{-1} - Δ
	if (SC){ // we need to symmetrize anomalous components
	  const int N2 = N/2;
	  MatrixXcd S12 = Sm.block(0, N2, N2, N2);
	  MatrixXcd S21 = Sm.block(N2, 0, N2, N2);
	  MatrixXcd A = 0.5*(S12 + S21);
	  Sm.block(0, N2, N2, N2) = A;
	  Sm.block(N2, 0, N2, N2) = A;
	}
	{
	  complex<double> I(0.0,1.0);
	  Eigen::MatrixXcd A = (0.5*I)* (Sm - Sm.adjoint());
	  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);
	  if (es.info() != Eigen::Success){
	    clog<<"WARNING Eigenvalues of self-energy could not be computed!"<<endl;
	  }else{
	    bool causality_violated = false;
	    Eigen::VectorXd evals = es.eigenvalues();
	    for (int i = 0; i < evals.size(); ++i){
	      if (evals(i)<0){
		causality_violated = true;
		clog<<"WARNING correcting causality S[om="<<om(iw)<<",ib="<<ib<<"] eval="<<evals(i)<<endl;
	      }
	      evals(i) = std::max(evals(i), small);
	    }
	    if (causality_violated){
	      const Eigen::MatrixXcd& V = es.eigenvectors();
	      Eigen::MatrixXcd Apsd = V * evals.asDiagonal() * V.adjoint();
	      // Clean up Hermitian symmetry again
	      Apsd = 0.5 * (Apsd + Apsd.adjoint());
	      Eigen::MatrixXcd H = complex<double>(0.5,0) * (Sm + Sm.adjoint());
	      Sm = H - I*Apsd;
	    }
	  }
	}
	MatrixXcd Ginv   = om_Eimp - Sm - Dm;    // G^{-1} = (ω - E_imp - Σ - Δ)
	MatrixXcd Gmc  = Ginv.inverse();         // G = (ω - E_imp - Σ - Δ)^{-1}
	// Storin Σ back into Blitz arrays
	vector<bool> component_set(cix.Nbath_unique_unique,false);
	for (int iorb=0; iorb<cix.Nbath; ++iorb){
	  if (bb.orb_ind[iorb][0]==ib){
	    int iib0 = cix.bath_sigind0(iorb);
	    double ib_sg = (iib0 < 0 ? -1 : 1);
	    int iib   = std::abs(iib0);
	    if (not component_set[iib]){
	      int j=bb.orb_ind[iorb][1], k=bb.orb_ind[iorb][2];
	      Sigm(iib,iw) = Sm(j,k)*ib_sg;  // just find the first occurence of this component.
	      SigD(iib,iw) = SmD(j,k)*ib_sg;
	      Gc(iib,iw) = Gmc(j,k)*ib_sg;
	      component_set[iib] = true;
	      // There might be a bit of a problem when G_{dn} is different than G_{up} in SC state, i.e.
	      // both SC and magnetism are present. In this case we would need to write into iw that
	      // correspond to -omega, which is not readily available.
	      // We will take care of that specific case when needed.
	    }
	  }
	}
	if (debugPrint){
	  out1<<om(iw)<<" ";
	  out2<<om(iw)<<" ";
	  out3<<om(iw)<<" ";
	  out4<<om(iw)<<" ";
	  out5<<om(iw)<<" ";
	  out6<<om(iw)<<" ";
	  for (int i=0; i<N; ++i)
	    for (int j=0; j<N; ++j){
	      out1 << Sm(i,j).real() << " " << Sm(i,j).imag()<< " ";
	      out2 << Gmc(i,j).real()<< " " << Gmc(i,j).imag()<< " ";
	      out3 << Dm(i,j).real() << " " << Dm(i,j).imag()<< " ";
	      out4 << Gm(i,j).real() << " " << Gm(i,j).imag()<< " ";
	      out5 << Fm(i,j).real() << " " << Fm(i,j).imag()<< " ";
	      out6 << om_Eimp(i,j).real()<< " " << om_Eimp(i,j).imag()<< " ";
	    }
	  out1<<endl;
	  out2<<endl;
	  out3<<endl;
	  out4<<endl;
	  out5<<endl;
	  out6<<endl;
	}
      }
    }
  }
}
