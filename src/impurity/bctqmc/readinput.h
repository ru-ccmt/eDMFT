#include <iostream>
#include <fstream>
#include "basics.h"
#include "interpolate.h"
#include "tinyMatmul.h"
#include "dynamic_mesh.h"

using namespace std;
using namespace blitz;

auto ReadPseudoSig00(string Sfilename, int Npseudo, MPIcls* mpi)
{
  struct retVal{
    Array<double,1> xn;
    Array<complex<double>,2> Sm;
    int Sinitialized;
  };
  retVal r;
  r.Sinitialized=1;
  //if (Sm.extent(0)!=Npseudo || Sm.extent(1)!=tau.size()) Sm.resize(Npseudo,tau.size()); // pseudoparticle self-energies
  if (mpi->rank==mpi->master){
    ifstream fsig0(Sfilename.c_str());
    string line;
    getline(fsig0,line); // comment
    if (!fsig0.good()){
      r.Sinitialized=0;
      goto endread;
    };
    deque<double> wom;
    deque<Array<double,1>> wSg;
    while(fsig0){
      double w;
      fsig0>>w;
      if (!fsig0 or fsig0.eof()) break;
      Array<double,1> Sg(Npseudo*2);
      for (int j=0; j<Sg.size(); j++) fsig0 >> Sg(j);
      getline(fsig0,line);
      wom.push_back(w);
      wSg.push_back(Sg);
    }
    (*mpi->out)<<"# Read "<<Sfilename<<" with "<<wom.size()<<" entries"<<endl;
    if (wom.size()==0){
      (*mpi->out)<<"ERROR : No file "<<Sfilename<<endl;
      r.Sinitialized=0;
    }
    r.xn.resize(wom.size());
    r.Sm.resize(Npseudo,wom.size());
    for (int i=0; i<wom.size(); ++i){
      r.xn(i) = wom[i];
      for (int ip=0; ip<Npseudo; ++ip)
	r.Sm(ip,i) = complex<double>(wSg[i](2*ip),wSg[i](2*ip+1));
    }
 }
 endread: mpi->bcast_ints(&r.Sinitialized,1);
  if (r.Sinitialized==0){
    (*mpi->out)<<"WARNING: File "<<Sfilename<<" could not be read"<<endl;
    return r;
  }else{
    mpi->bcast_Sm(r.xn,r.Sm,Npseudo);
    return r;
  }
}


auto ReadDelta_mpi(string fDelta, int Nbath, MPIcls& mpi)
{
  struct retVal{
    Array<double,1> om;
    Array<complex<double>,2> Deltaw;
  };
  //(*mpi->out) << "Before reading Delta"<<endl;
  retVal r;
  if (mpi.rank==mpi.master){
    ifstream fDelta_(fDelta.c_str());
    if (!fDelta_) {
      throw std::runtime_error("ReadDelta_mpi: cannot open file " + fDelta);
    }
    deque<double> wom;
    deque<vector<double>> wDelta;
    string line;
    int line_no = 0;
    while (std::getline(fDelta_, line)) {
      ++line_no;
      // skip empty or comment lines
      {
	string tmp = line;
	// trim leading spaces
	tmp.erase(0, tmp.find_first_not_of(" \t\r\n"));
	if (tmp.empty() || tmp[0] == '#') continue;
      }
      std::istringstream iss(line);
      double w;
      if (!(iss >> w)) {
	std::cerr << "ReadDelta_mpi: could not read frequency on line "
		  << line_no << " in file " << fDelta << '\n';
	throw std::runtime_error("ReadDelta_mpi: parse error");
      }

      // read the rest of the numbers on this line
      std::vector<double> cols;
      double val;
      while (iss >> val) {
	cols.push_back(val);
      }
      const int expected_cols = 2*Nbath;
      if (static_cast<int>(cols.size()) < expected_cols) {
	std::cerr << "ReadDelta_mpi: too few columns on line " << line_no
		  << " in file " << fDelta << ". Expected "
		  << expected_cols << " values after frequency, got "
		  << cols.size() << '\n';
	throw std::runtime_error("ReadDelta_mpi: insufficient columns");
      } else if (static_cast<int>(cols.size()) > expected_cols) {
	std::cerr << "ReadDelta_mpi: WARNING: extra columns on line " << line_no
		  << " in file " << fDelta << ". Expected "
		  << expected_cols << ", got "
		  << cols.size() << " — extra values will be ignored.\n";
      }
      // Store exactly 2*Nbath values
      //Array<double,1> wD(2*Nbath);
      //for (int i = 0; i < expected_cols; ++i) {
      //	wD(i) = cols[i];
      //}
      wom.push_back(w);
      wDelta.push_back(cols);
      /*
      double w;
      fDelta_>>w;
      Array<double,1> wD(2*Nbath);
      for (int i=0; i<wD.size(); ++i) fDelta_>>wD(i);
      if (!fDelta_ or fDelta_.eof()) break;
      getline(fDelta_,line);
      wom.push_back(w);
      wDelta.push_back(wD);
      */
    }
    (*mpi.out)<<"# Read "<<fDelta<<" with "<<wom.size()<<" entries"<<endl;
    r.om.resize(wom.size());
    for (int i=0; i<wom.size(); i++) r.om(i) = wom[i];
    r.Deltaw.resize(Nbath,wom.size());
    r.Deltaw=0;
    for (int i=0; i<wom.size(); i++){
      for (int ib=0; ib<std::min<int>(wDelta[i].size(),Nbath); ++ib)
	r.Deltaw(ib,i) = complex<double>(wDelta[i][2*ib],wDelta[i][2*ib+1]);
    }
  }
  mpi.BcastDeltaw(r.om,r.Deltaw,Nbath,fDelta);
  return r;
}

/*
int ReadDiagsFromFile(int norder, map<string,double>& stored2, MPIcls& mpi)
{
  Array<int,3> diags;
  Array<int,2> pseudos;
  int stored_size=0;
  if (mpi.rank==mpi.master){
    // In serial, we could just create has-table from file, but in parallel we need to first
    // store data in continous array, which are pseudos and diags.
    string filename = "diags."s+to_string(norder);
    ifstream inp(filename);
    deque<pair<Array<int,2>,Array<int,1>>> stored;
    Array<int,2> diag(norder,2);
    Array<int,1> pseudo(2*norder);
    int iin=0;
    while (inp and !inp.eof()){
      inp>>iin;
      if (inp.eof()) break;
      for (int j=0; j<norder; j++) inp>>diag(j,0)>>diag(j,1);
      for (int j=0; j<2*norder; j++) inp>>pseudo(j);
      stored.push_back(make_pair(diag.copy(),pseudo.copy()));
      //string s0 = diag_p_2String(diag,pseudo);
      //stored2[s0]=1;
    }
    (*mpi.out)<<"From file "<<filename<<" read "<<iin<<" diagrams."<<endl;
    diags.resize(stored.size(),norder,2);
    pseudos.resize(stored.size(),2*norder);
    for (int i=0; i<stored.size(); ++i){
      diags(i,Range::all(),Range::all()) = stored[i].first;
      pseudos(i,Range::all()) = stored[i].second;
    }
    stored_size=stored.size();
  }
  mpi.bcast_ints(&stored_size,1);
  if (stored_size==0) return 0;
  mpi.bcast_diags(diags,pseudos,stored_size,norder);
  Array<int,2> diag(norder,2);
  Array<int,1> pseudo(2*norder);
  for (int i=0; i<stored_size; ++i){
    diag=diags(i,Range::all(),Range::all());
    pseudo=pseudos(i,Range::all());
    string s0 = diag_p_2String(diag,pseudo);
    stored2[s0]=1;
  }
  return stored_size;
}
*/


////////////////////////////
// For imaginary axis data
////////////////////////////
void InterpolateIfNeeded(Array<double,1>& iom, Array<complex<double>,2>& Deltaw, double beta, ostream& out)
{
  if (fabs(iom(0)-M_PI/beta) > 1e-6){// temperature is different than the input.
    out<<"WARN: The input Delta seems to have different temperature, hence interpolating"<<" iom(0)="<<iom(0)<<" expected="<<M_PI/beta<<endl;
    int nom = static_cast<int>((iom(iom.size()-1)*beta/M_PI-1)/2.+0.5)+1; // last point
    Array<double,1> om(nom);
    Array<complex<double>,2> Delta(Deltaw.extent(0),nom);
    for (int i=0; i<nom; i++) om(i) = (2*i+1.)*M_PI/beta; // the new mesh corresponding to this temperature
    Spline1D<double> gspliner(iom.size());  // we will spline the results
    Spline1D<double> gsplinei(iom.size());  // we will spline the results
    for (int ib=0; ib<Deltaw.extent(0); ++ib){
      // Seeting up the spline
      for (int j=0; j<iom.size(); j++){
	gspliner[j] = Deltaw(ib,j).real(); // filling up the spline
	gsplinei[j] = Deltaw(ib,j).imag(); // filling up the spline
      }
      //double dfN = 0.0;
      //complex<double> df0 = (Deltaw(ib,1)-Deltaw(ib,0))/(iom(1)-iom(0));
      //gspliner.splineIt( iom, df0.real(), dfN);              // Splines the laplace data, knowing the derivatives at both ends
      //gsplinei.splineIt( iom, df0.imag(), dfN);              // Splines the laplace data, knowing the derivatives at both ends
      gspliner.splineIt( iom );
      gsplinei.splineIt( iom );
      for (int j=0; j<nom; ++j){
	intpar p = Interp(om(j),iom);
	Delta(ib,j) = complex<double>( gspliner(p) , gsplinei(p) );
      }
    }
    iom.resize(nom); iom = om;
    Deltaw.resize(Deltaw.extent(0),nom); Deltaw = Delta;
  }
}

struct HighFrequencyParams {
    double a;
    double eps;
};
//pair<double,double>
HighFrequencyParams FindHighFrequencyExp(int Nf, const Array<double,1>& omi, const Array<complex<double>,1>& fi)
{// Looking for the form : a/(iom-eps). We will determine a and eps and store into ah.first=a and ah.second=eps
  double S=0, Sy=0, Sz=0;
  for (int j=omi.size()-Nf; j<omi.size(); j++){
    double y = fi(j).imag()*omi(j);         // 1/omega    term
    double z = fi(j).real()*omi(j)*omi(j);  // 1/omega^2  term
    Sy += y;
    Sz += z;
    S++;
  }
  double ax = Sy/S;
  double bx = Sz/S;
  
  //pair<double,double> ah;
  //ah.first = -ax;    // a
  //ah.second = bx/ax; // eps
  HighFrequencyParams ah;
  ah.a = -ax;
  ah.eps=0;
  if (ax!=0) ah.eps = bx/ax;
  return ah;
}
void InverseFourier(const Array<complex<double>,1>& Giom, const Array<double,1>& iom, Array<double,1>& Gtau, const Array<double,1>& tau, const /*pair<double,double>*/HighFrequencyParams& ah, double beta)
{
  Array<complex<double>,1> dG(iom.size());
  for (int n=0; n<iom.size(); n++){
    complex<double> g_infty = ah.a/(complex<double>(0,iom(n))-ah.eps);
    dG(n) = Giom(n)-g_infty;
    //cout<<iom(n)<<" "<<dG(n).real()<<" "<<dG(n).imag()<<endl;
  }
  
  for (int t=0; t<tau.size(); t++){
    double tau_ = tau(t);
    double dsum=0;
    for (int n=0; n<iom.size(); n++)
      dsum += cos(iom(n)*tau_)*dG(n).real()+sin(iom(n)*tau_)*dG(n).imag();
    
    Gtau(t) = 2*dsum/beta;
    if (ah.eps>0)
      Gtau(t) -= ah.a*exp(-ah.eps*tau(t))/(1.+exp(-ah.eps*beta));
    else
      Gtau(t) -= ah.a*exp(ah.eps*(beta-tau(t)))/(1.+exp(ah.eps*beta));
  }
  Gtau *= -1.;// Here we want to make it positive in the interval [0,beta]
}
void FourierOfDelta(const Array<double,1>& iom, const Array<complex<double>,1>& Deltaw,
		    Array<double,1>& tau, Array<double,1>& Deltat, int Nhigh=40)
{
  double beta = tau(tau.extent(0)-1);
  int Nhighf = min(Nhigh, static_cast<int>(iom.size()*0.4)); // correction 21.4.2008 : Shuld not compute high frequency from low frequency points
  //pair<double,double> ah = FindHighFrequencyExp(Nhighf, iom, Deltaw);
  HighFrequencyParams ah = FindHighFrequencyExp(Nhighf, iom, Deltaw);
  //cout<<"beta="<<beta<<" Nhighf="<<Nhighf<<" ah="<<ah.a<<","<<ah.eps<<endl;
  InverseFourier(Deltaw, iom, Deltat, tau, ah, beta); 
}

void SetTranspose(Array<double,2>& Ht, const Array<double,2>& H){
  Ht.resize(H.extent(1),H.extent(0));
  for (int i=0; i<H.extent(1); ++i)
    Ht(i,Range::all()) = H(Range::all(),i);
}

void CreateSigma000(Array<complex<double>,2>& Sm,
		    const Array<double,1>& om, const Array<double,2>& Deltai, const Array<double,1>& feo, const Cix& cix, double beta, bool fastFilesystem, MPIcls& mpi){
  // If inteligence guess for the pseudo-particles self-energy is not found,
  // it creates a guess using atomic type of approximation.
  //
  // Sigma_b(omega) += |<f|c_s^\dagger|b>|^2 f(Ef-omega)f(-Ef)/f(-omega) * Delta(Ef-omega)
  // Sigma_f(omega) += |<f|c_s^\dagger|b>|^2 f(Eb-omega)f(-Eb)/f(-omega) * Delta(omega-E_b)
  bool debug=false;
  ostream& clog = *mpi.out;
  if (false){
    Array<double,2> G0(cix.Npseudo_unique,om.size()); // atomic-like G is stored.
    G0 = 0; // 0 it means it has not been set yet
    Array<double,1> work(3*cix.Nmaxsize); // for exact diagonalization
    Array<double,2> Ht, expE;
    Array<double,1> HH, E;
    Array<double,1> gw(om.size());
    Array<double,1> Eb_min(cix.Nblock); Eb_min=0;
    //vector<Array<double,3>> tmpG0;
    //if (saveG0) tmpG0.resize(cix.Nblock);
    //vector<Array<double,2>> H_temp(cix.Nblock);
    //vector<Array<double,1>> E_temp(cix.Nblock);
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size(); // block size
      bool need_to_compute=false;   // check if at least one G0 components is not yet computed.
      for (int i=0; i<N*N; ++i){
	int _ii_ = cix.Sgm_symmetry[ib](i/N,i%N); // use symmetry to not compute the same componente multiple times.
	int ii = abs(_ii_)-1;
	if (_ii_!=0 and G0(ii,0)==0){
	  need_to_compute=true;
	  break;
	}
      }
      //cout<<"ib="<<ib<<" need_to_compute="<<need_to_compute<<endl;
      if (need_to_compute){
	if (N==1){ // no need to perform exact diagonalization, because the block is one dimensional
	  double E = cix.Ham[ib](0,0)-cix.Emin;
	  int _ii_ = cix.Sgm_symmetry[ib](0,0); // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
	  int ii = abs(_ii_)-1;
	  double fm = fermi(-E*beta); 
	  for (int i=0; i<om.size(); ++i){
	    double ff = (om(i)>0) ? fermi((E-om(i))*beta)*fm/(1-feo(i)) : fermi((om(i)-E)*beta)*(1-fm)/feo(i);// f(Ef-w)*f(-Ef)/f(-w)
	    G0(ii,i) = sign(_ii_)*M_PI*ff; // note that all exponents need to be negative, hence subtracting Emin
	  }
	  Eb_min(ib) = E;
	  // this is normally handled by lambda to adjust the common origin from which we meassure atomic states.
	}else{
	  Array<double,2> H(cix.Ham[ib].copy()); // copy constructor. We don't want to modify cix.Ham. It will contain eigenvectors.
	  E.resize(N);                    // store eigenvalues
	  xsyev(H.extent(0), H.data(), H.extent(1), E.data(), work.data(), work.size(), "V");
	  // Note that diagonalization gives: Ham = H^+ * E * H
	  // In our case
	  // G = H^+ exp(-E*t) * H => G(i,j) = H^H(i,k) exp(-E(k)*t) H(k,j) = Ht^*(i,k) exp(-E(k)*t) Ht(j,k)
	  expE.resize(N,om.size());  // compute expE == exp(-(E-Emin)*tau)
	  for (int k=0; k<N; ++k){
	    double Ek = E(k)-cix.Emin;
	    double fm = fermi(-Ek*beta); 
	    for (int i=0; i<om.size(); ++i)
	      expE(k,i) = (om(i)>0) ? fermi((Ek-om(i))*beta)*fm/(1-feo(i)) : fermi((om(i)-Ek)*beta)*(1-fm)/feo(i);// f(Ef-w)*f(-Ef)/f(-w)
	  }
	  expE *= M_PI;
	  Eb_min(ib) = min(E)-cix.Emin;
	  //if (saveG0){
	  //  tmpG0[ib].resize(N,N,tau.size());
	  //  tmpG0[ib]=0;
	  //}
	  //Array<double,1> EE(E.copy()); EE= EE-cix.Emin;
	  //cout<<"ib="<<ib+1<<" E-Emin="<<EE<<" expE(-1)="<<expE(Range::all(),expE.ubound(1))<<endl;
	  SetTranspose(Ht,H); // It is more efficient to compute the product with Ht than H.
	  HH.resize(N);       // will contain H^+ * H
	  for (int i=0; i<N; ++i){
	    for (int j=0; j<N; ++j){
	      int _ii_ = cix.Sgm_symmetry[ib](i,j);
	      int ii = abs(_ii_)-1;
	      //if (saveG0){
	      //  for (int k=0; k<N; ++k) HH(k) = Ht(i,k)*Ht(j,k); // if Ham is complex, we would need to conjugate first Ht(i,k) 
	      //  gw = 0;
	      //  for (int k=0; k<N; ++k) gw(Range::all()) += HH(k)*expE(k,Range::all());
	      //  tmpG0[ib](i,j,Range::all()) = gw(Range::all());
	      //}
	      if (_ii_!=0 and G0(ii,0)==0){
		// this entry is still needed
		for (int k=0; k<N; ++k) HH(k) = Ht(i,k)*Ht(j,k); // if Ham is complex, we would need to conjugate first Ht(i,k) 
		//// if we want to check that HH*E gives original Ham.
		//double eh=0;
		//for (int k=0; k<N; ++k) eh += HH(k)*E(k);
		//
		gw = 0;
		for (int k=0; k<N; ++k) gw(Range::all()) += HH(k)*expE(k,Range::all());
		G0(ii,Range::all()) = sign(_ii_)*gw; // setting G0
	      }
	    }
	  }
	  //H_temp[ib].resize(N,N); H_temp[ib]=H;
	  //E_temp[ib].resize(N); E_temp[ib]=E;
	}
      }else{
	int _ii_ = cix.Sgm_symmetry[ib](0,0);
	int ii = abs(_ii_)-1;
	int ib_equiv = cix.Sgm_indx[ii][0][0];
	Eb_min(ib) = Eb_min(ib_equiv);
      }
    }
    //cout<<"Eb_min="<<Eb_min<<endl;

    if (mpi.rank==mpi.master and debug){
      ofstream out("PPG0_start"); out.precision(12);
      for (int iw=0; iw<om.size(); ++iw){
	out<<om(iw)<<" ";
	for (int ip=0; ip<cix.Npseudo_unique; ++ip) out<< G0(ip,iw)<<" ";
	out<<endl;
      }
    }

    Sm.resize(cix.Npseudo_unique,om.size()); // resizing output self-energy
    complex<double> I(0,1);
    Sm=0;
    for (int ip=0; ip<cix.Npseudo_unique; ++ip){ // we go over unique self-energy components, rather than all possible components
      TinyVector<int,3> bij = cix.Sgm_indx[ip][0];// gives block, and 2D index (i,j) inside block
      if (bij[0]<0) cout<<"ERROR ib="<<bij[0]<<" should be positive. Check algorithm"<<endl;
      int i1 = cix.Block[bij[0]](bij[1]);         // these are indices of atomic states started from 1...2^(orb*spin)
      int i2 = cix.Block[bij[0]](bij[2]);         // number of this atomic state as found in cix file
      for (int iorb=0; iorb<cix.Nbath; ++iorb){   // over all baths for back propagating NCA diagram
	int j1 = cix.Fi(2*iorb, i1);              // j1 and j2 are atomic states after applying |j1>=psi^+_{iorb}|i1>
	int j2 = cix.Fi(2*iorb, i2);              // |j2>=psi^+_{iorb}|i2>
	int s = cix.bath_sigind(iorb);            // which unique bath
	if (j1!=0 and j2!=0){                     // is j1=0 or j2=0 diagram vanishes.
	  int ib = cix.iBlck(abs(j1));            // ib is block in which j1 and j2 are located. Could be obtained by cix.Fbc(2*iorb,,bij[0])
	  int _ii_ = cix.Sgm_symmetry[ib-1](cix.inBlock(abs(j1)),cix.inBlock(abs(j2))); // indec of G[ii]==G(j1,j2)
	  int ii = abs(_ii_)-1;
	  if (ii>=0){
	    complex<double> tsgn = sign(_ii_)*sign(j1)*sign(j2); tsgn *= I;
	    double E0 = Eb_min(ib-1);
	    int pos = om.size()-1;
	    for (int iw=0; iw<om.size(); ++iw){
	      intpar ipw = InterpRight(E0-om(iw),pos,om);
	      Sm(ip,iw) += tsgn*G0(ii,iw)*linInterp(Deltai,s,ipw);
	    }
	  }
	  if (debug) cout<<"  back: ("<<i1<<","<<i2<<") iorb="<<iorb<<" s="<<sign(j1)*sign(j2)<<" ("<<j1<<","<<j2<<")="<<ii<<" ib="<<ib<<endl;
	}
	// Now continue with forward propagating hybridization
	j1 = cix.Fi(2*iorb+1, i1);
	j2 = cix.Fi(2*iorb+1, i2);
	if (j1!=0 and j2!=0){
	  int ib = cix.iBlck(abs(j1));
	  int _ii_ = cix.Sgm_symmetry[ib-1](cix.inBlock(abs(j1)),cix.inBlock(abs(j2))); // G(j1,j2);
	  int ii = abs(_ii_)-1;
	  if (ii>=0){
	    complex<double> tsgn = sign(_ii_)*sign(j1)*sign(j2); tsgn *= I;
	    double E0 = Eb_min(ib-1);
	    int pos = 0;
	    for (int iw=0; iw<om.size(); ++iw){
	      intpar ipw = InterpLeft(om(iw)-E0,pos,om);
	      Sm(ip,iw) += tsgn*G0(ii,iw)*linInterp(Deltai,s,ipw);
	    }
	  }
	  if (debug) cout<<"  forw: ("<<i1<<","<<i2<<") iorb="<<iorb<<" s="<<sign(j1)*sign(j2)<<" ("<<j1<<","<<j2<<")="<<ii<<" ib="<<ib<<endl;
	}
      }
      const double dsmall = 1e-1, L2=om(0)*om(0);
      for (int iw=0; iw<om.size(); iw++){
	Sm(ip,iw) += I*dsmall*exp(-om(iw)*om(iw)/(0.05*L2));
	Sm(ip,iw) *= exp(-om(iw)*om(iw)/(0.2*L2));
      }
    }
  }else{
    //cout<<"Starting to create Sigma000"<<endl;
    //cout<<"cix.Nblock="<<cix.Nblock<<" Npseudo_unique="<<cix.Npseudo_unique<<" om="<<om.size()<<endl;
    Sm.resize(cix.Npseudo_unique,om.size()); // resizing output self-energy
    //cout<<"resized="<<Sm.shape()<<endl;
    Sm=0;
    // This is simplest algorithm, which sets Im(Sm) as gaussian 
    const double L2=om(0)*om(0);
    const complex<double> I(0.0,1.0);
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size(); // block size
      for (int i=0; i<N; ++i){
	int _ii_ = cix.Sgm_symmetry[ib](i,i);
	int ii = abs(_ii_)-1;
	if (ii>=0){
	  for (int iw=0; iw<om.size(); iw++)
	    Sm(ii,iw) = 0.3*exp(-om(iw)*om(iw)/(0.2*L2))*I;
	}
      }
    }
  }
  
  Array<double,1> dh(om.size()); dh = GiveDh(om);
  Array<double,1> logo(om.size()); logo = CalcLogoOnMesh(om);
  Array<double,1> Smi(om.size()), Smr(om.size());
  
  for (int l=0; l<cix.Npseudo_unique; l++){
    for (int i=0; i<om.size(); i++) Smi(i) = -Sm(l,i).imag()*(1-feo(i));
    Smr = KramarsKronig(Smi, om, dh, logo);
    for (int i=0; i<om.size(); i++) Sm(l,i) = complex<double>(Smr(i), Sm(l,i).imag());
  }
  
  if (mpi.rank==mpi.master){//DEBUG and fastFilesystem){
    ofstream out("PPSig_start"); out.precision(12);
    for (int iw=0; iw<om.size(); ++iw){
      out<<om(iw)<<" ";
      for (int ip=0; ip<cix.Npseudo_unique; ++ip) out<< Sm(ip,iw).real()<<" "<<Sm(ip,iw).imag()<<"  ";
      out<<endl;
    }
  }
  if (mpi.rank==mpi.master) clog<<"WARNING: Creating a guess for Starting Pseudo self-energy in imaginary time"<<endl;
}

void CreateSigma000(Array<double,2>& Sigt,
		    const Array<double,1>& tau, const Spline1DMany& Deltas, //const vector<Spline1D<double>> Deltas,
		    const Cix& cix, double beta, bool fastFilesystem, MPIcls& mpi){
  // This is for imaginary axis calculation
  const bool DEBUG=false;
  const bool saveG0=true;
  ostream& clog = *mpi.out;
  // If previous pseudo-particles self-energy is not found, it creates a guess using atomic type NCA approximation.
  // We assume that G0(pp,tau) = exp(-E_pp*tau) and evaluate NCA diagrams with it.
  // In the case of one band, we would have
  // Sigma_b(omega) += |<f|c_s^\dagger|b>|^2 exp(-(Ef-E_gs)*tau) * Delta(-tau)
  // Sigma_f(omega) += |<f|c_s^\dagger|b>|^2 exp(-(Eb-E_gs)*tau) * Delta( tau)
  Array<double,2> G0(cix.Npseudo_unique,tau.size()); // atomic-like G is stored.
  G0 = 0; // 0 it means it has not been set yet
  Array<double,1> work(3*cix.Nmaxsize); // for exact diagonalization
  Array<double,2> Ht, expE;
  Array<double,1> HH, E;
  Array<double,1> gt(tau.size());
  vector<Array<double,3>> tmpG0;
  if (saveG0) tmpG0.resize(cix.Nblock);
  //vector<Array<double,2>> H_temp(cix.Nblock);
  //vector<Array<double,1>> E_temp(cix.Nblock);
  for (int ib=0; ib<cix.Nblock; ++ib){
    int N = cix.Block[ib].size(); // block size
    bool need_to_compute=false;   // check if at least one G0 components is not yet computed.
    for (int i=0; i<N*N; ++i){
      int _ii_ = cix.Sgm_symmetry[ib](i/N,i%N); // use symmetry to not compute the same componente multiple times.
      int ii = abs(_ii_)-1;
      if (_ii_!=0 and G0(ii,0)==0){
	need_to_compute=true;
	break;
      }
    }
    //cout<<"ib="<<ib<<" need_to_compute="<<need_to_compute<<endl;
    if (need_to_compute){
      if (N==1){ // no need to perform exact diagonalization, because the block is one dimensional
	double E = cix.Ham[ib](0,0);
	int _ii_ = cix.Sgm_symmetry[ib](0,0); // note that cix.Sgm_symmetry starts with 1, hence need to subtract 1.
	int ii = abs(_ii_)-1;
	for (int it=0; it<tau.size(); ++it) G0(ii,it) = sign(_ii_)*exp(-(E-cix.Emin)*tau(it)); // note that all exponents need to be negative, hence subtracting Emin
	// this is normally handled by lambda to adjust the common origin from which we meassure atomic states.
      }else{
	Array<double,2> H(cix.Ham[ib].copy()); // copy constructor. We don't want to modify cix.Ham. It will contain eigenvectors.
	E.resize(N);                    // store eigenvalues
	xsyev(H.extent(0), H.data(), H.extent(1), E.data(), work.data(), work.size(), "V");
	// Note that diagonalization gives: Ham = H^+ * E * H
	// In our case
	// G = H^+ exp(-E*t) * H => G(i,j) = H^H(i,k) exp(-E(k)*t) H(k,j) = Ht^*(i,k) exp(-E(k)*t) Ht(j,k)
	expE.resize(N,tau.size());  // compute expE == exp(-(E-Emin)*tau)
	for (int k=0; k<N; ++k){
	  if (E(k)<cix.Emin) cout<<"WARNING cix.Emin="<<cix.Emin<< " should be smaller than E(k)="<<E(k)<<endl;
	  for (int it=0; it<tau.size(); ++it)
	    expE(k,it) = exp(-(E(k)-cix.Emin)*tau(it));
	}
	if (saveG0){
	  tmpG0[ib].resize(N,N,tau.size());
	  tmpG0[ib]=0;
	}
	Array<double,1> EE(E.copy()); EE= EE-cix.Emin;
	//cout<<"ib="<<ib+1<<" E-Emin="<<EE<<" expE(-1)="<<expE(Range::all(),expE.ubound(1))<<endl;
	SetTranspose(Ht,H); // It is more efficient to compute the product with Ht than H.
	HH.resize(N);       // will contain H^+ * H
	for (int i=0; i<N; ++i){
	  for (int j=0; j<N; ++j){
	    int _ii_ = cix.Sgm_symmetry[ib](i,j);
	    int ii = abs(_ii_)-1;

	    if (saveG0){
	      for (int k=0; k<N; ++k) HH(k) = Ht(i,k)*Ht(j,k); // if Ham is complex, we would need to conjugate first Ht(i,k) 
	      gt = 0;
	      for (int k=0; k<N; ++k) gt(Range::all()) += HH(k)*expE(k,Range::all());
	      tmpG0[ib](i,j,Range::all()) = gt(Range::all());
	    }
	    
	    if (_ii_!=0 and G0(ii,0)==0){
	      // this entry is still needed
	      for (int k=0; k<N; ++k) HH(k) = Ht(i,k)*Ht(j,k); // if Ham is complex, we would need to conjugate first Ht(i,k) 
	      //// if we want to check that HH*E gives original Ham.
	      //double eh=0;
	      //for (int k=0; k<N; ++k) eh += HH(k)*E(k);
	      //
	      gt = 0;
	      for (int k=0; k<N; ++k) gt(Range::all()) += HH(k)*expE(k,Range::all());
	      G0(ii,Range::all()) = sign(_ii_)*gt; // setting G0

	    }
	  }
	}
	//H_temp[ib].resize(N,N); H_temp[ib]=H;
	//E_temp[ib].resize(N); E_temp[ib]=E;
      }
    }
  }

  if (mpi.rank==mpi.master and saveG0){
    for (int ib=0; ib<cix.Nblock; ++ib){
      int N = cix.Block[ib].size();
      if (N<2 or tmpG0[ib].extent(0)==0) continue;
      ofstream out("checkPGG0."+to_string(ib+1));
      for (int it=0; it<tau.size(); ++it){
	out<<tau(it)<<" ";
	for (int i=0; i<N; ++i){
	  for (int j=0; j<N; ++j){
	    out<<tmpG0[ib](i,j,it)<<" ";
	  }
	}
	out<<endl;
      }
    }
  }
  
  // using spline to evaluate Delta(beta-tau)==mDelta and Delta(tau)==pDelta
  Array<double,2> mDelta(Deltas.size(),tau.size()), pDelta(Deltas.size(),tau.size());
  for (int s=0; s<Deltas.size(); ++s){
    for (int it=0; it<tau.size(); ++it) pDelta(s,it) = Deltas(it,s);//Deltas[s][it];
    for (int it=0; it<tau.size(); ++it) mDelta(s,it) = Deltas[s](Interp(beta-tau(it),tau));
  }
  if (mpi.rank==mpi.master){//DEBUG){
    ofstream test("checkDelta.dat");
    for (int it=0; it<tau.size(); ++it){
      test<<tau(it)<<" ";
      for (int s=0; s<Deltas.size(); ++s)
	test<<mDelta(s,it)<<" ";
      //test<<mDelta(it,s)<<" ";
      test<<"  ";
      for (int s=0; s<Deltas.size(); ++s)
	test<<pDelta(s,it)<<" ";
      //test<<pDelta(it,s)<<" ";
      test<<endl;
    }
  }
  Sigt.resize(cix.Npseudo_unique,tau.size()); // resizing output self-energy
  Sigt=0;
  for (int ip=0; ip<cix.Npseudo_unique; ++ip){ // we go over unique self-energy components, rather than all possible components
    TinyVector<int,3> bij = cix.Sgm_indx[ip][0];// gives block, and 2D index (i,j) inside block
    if (bij[0]<0) cout<<"ERROR ib="<<bij[0]<<" should be positive. Check algorithm"<<endl;
    int i1 = cix.Block[bij[0]](bij[1]);         // these are indices of atomic states started from 1...2^(orb*spin)
    int i2 = cix.Block[bij[0]](bij[2]);         // number of this atomic state as found in cix file
    for (int iorb=0; iorb<cix.Nbath; ++iorb){   // over all baths for back propagating NCA diagram
      int j1 = cix.Fi(2*iorb, i1);              // j1 and j2 are atomic states after applying |j1>=psi^+_{iorb}|i1>
      int j2 = cix.Fi(2*iorb, i2);              // |j2>=psi^+_{iorb}|i2>
      int s = cix.bath_sigind(iorb);            // which unique bath
      if (j1!=0 and j2!=0){                     // is j1=0 or j2=0 diagram vanishes.
	int ib = cix.iBlck(abs(j1));            // ib is block in which j1 and j2 are located. Could be obtained by cix.Fbc(2*iorb,,bij[0])
	int _ii_ = cix.Sgm_symmetry[ib-1](cix.inBlock(abs(j1)),cix.inBlock(abs(j2))); // indec of G[ii]==G(j1,j2)
	int ii = abs(_ii_)-1;
	if (ii>=0)
	  Sigt(ip,Range::all()) += sign(_ii_)*G0(ii,Range::all())*mDelta(s,Range::all())*sign(j1)*sign(j2);
	if (DEBUG) cout<<"  back: ("<<i1<<","<<i2<<") iorb="<<iorb<<" s="<<sign(j1)*sign(j2)<<" ("<<j1<<","<<j2<<")="<<ii<<" ib="<<ib<<endl;
      }
      // Now continue with forward propagating hybridization
      j1 = cix.Fi(2*iorb+1, i1);
      j2 = cix.Fi(2*iorb+1, i2);
      if (j1!=0 and j2!=0){
	int ib = cix.iBlck(abs(j1));
	int _ii_ = cix.Sgm_symmetry[ib-1](cix.inBlock(abs(j1)),cix.inBlock(abs(j2))); // G(j1,j2);
	int ii = abs(_ii_)-1;
	if (ii>=0)
	  Sigt(ip,Range::all()) += sign(_ii_)*G0(ii,Range::all())*pDelta(s,Range::all())*sign(j1)*sign(j2);
	if (DEBUG) cout<<"  forw: ("<<i1<<","<<i2<<") iorb="<<iorb<<" s="<<sign(j1)*sign(j2)<<" ("<<j1<<","<<j2<<")="<<ii<<" ib="<<ib<<endl;
      }
    }
  }
  if (mpi.rank==mpi.master){//DEBUG and fastFilesystem){
    ofstream out("PPSig_start"); out.precision(12);
    for (int it=0; it<tau.size(); ++it){
      out<<tau(it)<<" ";
      for (int ip=0; ip<cix.Npseudo_unique; ++ip) out<< Sigt(ip,it)<<" ";
      out<<endl;
    }
  }
  if (mpi.rank==mpi.master) clog<<"WARNING: Creating a guess for Starting Pseudo self-energy in imaginary time"<<endl;
}
