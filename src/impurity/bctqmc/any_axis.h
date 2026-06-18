class AnyAxis{
public:
  double beta;
  AnyAxis(double _beta_) : beta(_beta_){};
  virtual ~AnyAxis() = default;

  // Unified access to “the frequency grid” for printing / templated code
  virtual const Array<double,1>& freq() const = 0;

  void Print_Physical1(const string& filename, const Array<complex<double>,2>& Gc) const
  {
    const auto& w = freq();
    ofstream sout(filename); sout.precision(12);
    for (int i=0; i<w.size(); ++i){
      sout << w(i) << " ";
      for (int ip=0; ip<Gc.extent(0); ++ip) sout<<Gc(ip,i).real()<<" "<<Gc(ip,i).imag()<<" ";
      sout<<endl;
    }
  }
  void Print_Physical2(const string& filename, const Array<complex<double>,2>& Gc,
		       const Array<double,1>& density_matrix0, double nd, double mu, double Fimp, double Epot, double TrSigmaG, double Ekin, double nf,
		       const Array<double,1>& density_matrix, const Array<int,1>& deg_bath) const
  {
    const auto& w = freq();
    ofstream sout(filename); sout.precision(12);
    // Internally we keep reduced, symmetry-averaged bath occupancies.  For the
    // user-facing header we print the physical occupancies per reduced channel,
    // i.e. multiplied by the number of equivalent bath orbitals represented by
    // that channel.
    //
    // The pseudo/density-matrix estimate (nd, density_matrix0) is numerically
    // more stable than the Matsubara-summed occupancy (nf, density_matrix), so
    // we expose the stable quantity under the legacy nf/mom header names used
    // by the Python post-processing, while still preserving the frequency-sum
    // estimate under nd/density_matrix.
    sout << "# nf="<<nd<<" mu="<<mu<<" T="<<1.0/beta<<" TrSigmaG="<<TrSigmaG<<" Epot="<<Epot<<" Ekin="<<Ekin
	 << " mom=["<<density_matrix0(0)*deg_bath(0);
    for (int a=1; a<density_matrix0.extent(0); ++a) sout<<","<<density_matrix0(a)*deg_bath(a);
    sout<<"] Fimp="<<Fimp<<" nd="<<nf<<" density_matrix=["<<density_matrix(0)*deg_bath(0);
    for (int a=1; a<density_matrix.extent(0); ++a) sout<<","<<density_matrix(a)*deg_bath(a);
    sout<<"] "<<endl;
    
    for (int i=0; i<w.size(); ++i){
      sout << w(i) << " ";
      for (int ip=0; ip<Gc.extent(0); ++ip) sout<<Gc(ip,i).real()<<" "<<Gc(ip,i).imag()<<" ";
      sout<<endl;
    }
  }
};

void ReadArray(istream& fsig0, Array<double,1>& Sg){
  for (int j=0; j<Sg.size(); j++) fsig0 >> Sg(j);
}
void ReadArray(istream& fsig0, Array<complex<double>,1>& Sg){
  for (int j=0; j<Sg.size(); j++){
    double r, i;
    fsig0 >> r >> i;
    Sg(j) = complex<double>(r,i);
  }
}

template<class T>
auto _ReadPseudoSig00_(string Sfilename, int Npseudo, MPIcls& mpi)
{// T==complex<double> or T==double for real or imaginary axis
  struct retVal{
    Array<double,1> xn;
    Array<T,2> Sm;
    int Sinitialized;
  };
  retVal r;
  r.Sinitialized=1;
  if (mpi.rank==mpi.master){
    ifstream fsig0(Sfilename.c_str());
    string line;
    getline(fsig0,line); // possible comment
    if (!fsig0.good()){
      r.Sinitialized=0;
      goto endread;
    };
    deque<double> wom;
    deque<Array<T,1>> wSg;
    {// checking if line == comment
      std::size_t found = line.find("#");
      if (found!=std::string::npos) line.erase(found);
      line = trim(line);
      if (not line.empty()){ // line was not a comment
	std::stringstream inp2;
	inp2 << line << std::endl;
	double w;
	inp2 >> w;
	Array<T,1> Sg(Npseudo);
	ReadArray(inp2, Sg);
	wom.push_back(w);
	wSg.push_back(Sg);
      }
    }
    while(fsig0){
      double w;
      fsig0>>w;
      if (!fsig0 or fsig0.eof()) break;
      Array<T,1> Sg(Npseudo);
      ReadArray(fsig0, Sg);
      //for (int j=0; j<Sg.size(); j++) fsig0 >> Sg(j);
      getline(fsig0,line);
      wom.push_back(w);
      wSg.push_back(Sg);
    }
    (*mpi.out)<<"# Read "<<Sfilename<<" with "<<wom.size()<<" entries"<<endl;
    if (wom.size()==0){
      (*mpi.out)<<"ERROR : No file "<<Sfilename<<endl;
      r.Sinitialized=0;
    }
    r.xn.resize(wom.size());
    r.Sm.resize(Npseudo,wom.size());
    for (int i=0; i<wom.size(); ++i){
      r.xn(i) = wom[i];
      for (int ip=0; ip<Npseudo; ++ip)
	r.Sm(ip,i) = wSg[i](ip);
    }
  }
 endread: mpi.bcast_ints(&r.Sinitialized,1);
  if (r.Sinitialized==0){
    (*mpi.out)<<"WARNING: File "<<Sfilename<<" could not be read"<<endl;
    return r;
  }else{
    mpi.bcast_Sm(r.xn,r.Sm,Npseudo);
    return r;
  }
}
void Prints(ostream& out, const complex<double>& Gm){
  out<<Gm.real()<<" "<<Gm.imag()<<" ";
}
void Prints(ostream& out, const double& Gm){
  out<<Gm<<" ";
}

/*
void ReadStatusFile_any(ifstream& inp, bool& read_succeeded, Array<double,1>& freq_times1, Array<int,2>& cdiag, Array<int,1>& cpseudo, Array<int,1>& cDlt_type,
			MPIcls& mpi, double ctf_beta, int Nsize, int Norder, int N_min_order, int N_max_order)
{
  read_succeeded=true;
  Array<int,2> diag(Norder,2);
  Array<int,1> pseudo(2*Norder), Dlt_type(Norder);
  Array<double,1> freq_times(Nsize);

  for (int i=0; i<Norder; ++i) inp >> diag(i,0) >> diag(i,1);
  for (int i=0; i<2*Norder; ++i) inp >> pseudo(i);
  for (int i=0; i<Norder; ++i) inp >> Dlt_type(i);
  for (int i=0; i<Nsize; ++i){
    inp >> freq_times(i);
    if (fabs(freq_times(i))>ctf_beta){
      (*mpi.out)<<"Note: status file contains frequency/time outside cutoff="<<ctf_beta<<", hence aborting reading" << endl;
      read_succeeded=false;
      return;
    }
  }
  read_succeeded = read_succeeded and inp and not inp.eof();
  if (read_succeeded){
    cdiag.resize(Norder,2); cdiag=diag;
    cpseudo.resize(2*Norder); cpseudo=pseudo;
    cDlt_type.resize(Norder); cDlt_type=Dlt_type;
    freq_times1(Range(0,Nsize-1))=freq_times(Range(0,Nsize-1));
  }
}
*/
