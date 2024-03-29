// @Copyright 2007 Kristjan Haule
// 
template <int b>
int sigP(double Ps)
{
  return Ps>b ? 1 : -1;
}

class NOperators;

class TBath{
public:
  int ifl,bind;
  TBath(int ifl_=-1, int bind_=0): ifl(ifl_), bind(bind_) {}
};

class Njj{
public:
  function2D<double> M;
  int istate, i_m_state;
  int ifl1, ifl2, ii;
};


class Uentry{
public:
  int a, b, c, d;
  double u;
  Uentry(int a_=0, int b_=0, int c_=0, int d_=0, double u_=0) : a(a_), b(b_), c(c_), d(d_), u(u_){};
};

class ClusterData{
public:
  function2D<int> F_i;                  //yes
private:
  function2D<function2D<double> > F_M;  //yes
public:
  //function2D<double> Ene;               //yes
  function1D<double> Enes;              // Small energy array, which has only existing states
  function2D<int> Eind;                 // Index to the Energy array, which will also be used for exponents
private:
  function2D<double> Ene0;              //yes
  function2D<double> Spin;              //yes
  function1D<int> Ks;                   //yes
  function1D<int> msize_;               //yes
  function2D<double> Patom_;            //no
  function1D<double> Patom;             //no
  vector<function2D<double> > Ma;       //no
  vector<int> gflip_index;              //yes
public:
  int Nm, nsize, N_ifl, N_flavors, Nvfl;//yes
  int max_size, Osize, DOsize;          //yes
  function1D<double> Sz;                //yes
  function1D<int> ifl_dim;              //yes
  function1D<int> ifl_from_fl;          //yes
  function1D<int> bfl_from_fl;          //yes
  function2D<int> fl_from_ifl;          //yes
  vector<function2D<int> > tfl_index;   //yes
  vector<function2D<int> > vfl_index;   //yes
  vector<deque<int> > v2fl_index;       //yes
  vector<TBath> vfli_index;             //yes
  vector<deque<int> > sign;             //yes
  vector<deque<int> > conjg;            //yes
  vector<deque<int> > bfl_index;        //yes
  Number Zatom;                         //no
  function1D<function2D<function2D<double> > > HF_M; //yes
  function2D<int> DF_i, DF_inv;         // (DOsize,nsize+1)
  function2D<function2D<double> > DF_M; // (DOsize,nsize+1,size_new,size_old);
  function1D<double> natom;             //no
  int N_unique_fl;                      //yes
  function1D<int> fl_deg;               //yes
  function1D<int> Ns;                   //yes
  function1D<double> epsk;              //yes
  string RealSigma, ImagSigma;          //yes
  map<int,double> Id;                   //yes
  function2D<int> gflip_fl;             //yes
  function1D<pair<int,int> > gflip_ifl; //yes
  function2D<bool> ifl_equivalent;      //yes
  function1D<int> gs_ind;               //no
  bool QHB1;                            //yes
  int cnsize, cmax_size;                //no
  function1D<int> cNs;                  //no
  function1D<int> cmsize_;              //no
  function2D<int> cF_i;                 //no
  function1D<int> cindx;                //no
  vector<vector<double> > cEne;         //no
  function2D<function2D<double> > cF_M; //no
  //function2D<function2D<double> > Nj;   //yes
  //function1D<function2D<double> > Uc;   //yes
  deque<Uentry>  UC;
  //vector<deque<Njj> > Njjs;             // no
  function2D<function2D<double> > Njm;    //
  function2D<char> Njm_c;  // stores information if Njm(istate,fl2) is zero or idenity
  //function1D<char> Njm_z;  // stores information if all Njm(istate,:) are zero or idenity
  function2D<double> Njm_r; // ratio
  //function2D<double> Njs;  // for segment
  //deque<pair<int,int> > fl_fl;
  //
public:
  void ParsInput(ifstream& gout, ostream& clog);
  void EvaluateMa(double beta, double mu, double U);
  int size() const {return nsize;}
  int msize(int i) const {return msize_[i];}
  int praState(int i) const {return i+1;}
  int Fi(int iop, int ist) const {return F_i(iop,ist);}
  const funProxy<int>& Fi(int iop) const {return F_i[iop];}
  const function2D<double>& FM(int op, int st) const { return F_M(op,st);}
  const funProxy<function2D<double> >& FM(int op) const {return F_M[op];}
  double patom(int i) const {return Patom[i];}
  int N_baths(int ifl) const {return bfl_index[ifl].size();}
  void RecomputeCix(function2D<double>& AProb, double asign, double treshold, function1D<NState>& praStates);
  friend class NOperators;
  void Read_for_HB1(ifstream& gout, double mu, double U, ostream& clog);
  void HB1(double beta, double mu, const mesh1D& iom, const function2D<double>& AProb, int nom, int aom, const function2D<dcomplex>& Delta, function2D<dcomplex>& Sigma, const function1D<bool>& nonzero, const function1D<bool>& nonzerou, ostream& clog, double sderiv, int nom_small);
  void HB3(double beta, double mu, double U, const mesh1D& iom, const function2D<double>& AProb, int nom, int aom, const function2D<dcomplex>& Delta, function2D<dcomplex>& Sigma, const function1D<bool>& nonzero, const function1D<bool>& nonzerou, ostream& clog, double sderiv, int nom_small);

  void Check(function1D<NState>& praStates);
  double P_atom(int i, int m) const {return Patom_(i+1,m);}
  int BcastClusterData(int my_rank, int Master);
  void Compute_Nmatrices();
  void Compute_Nmatrices2(double U, bool Print);
  void Set_Old_Data_for_HB1(ostream& clog);
  void Construct_Ucf(function4D<double>& Ucf, double U) const;
  //void Get_HB2_Uc(function2D<function2D<double> >& Uc, double U);
};

void ClusterData::ParsInput(ifstream& gout, ostream& yout)
{//*****************************************************************************//
 // Many arrays defined in the following routine are explained in the example of 
 // 2x2 cluster DMFT for superconductings state. In the latter case, the following
 // input is read from the input 'cix' file:
 //
 // ifl  ifl_dim   bfl_index     gflip_index 
 //  0     1	    0                0
 //  1	   1       -0*               1
 //  2	   2        1  3  3 -1*      2
 //  3	   2        1 -3 -3 -1*      3
 //  4     1        2                4
 //  5     1       -2*               5
 //
 // The following arrays and variables are subsequently defined in this routine:
 // N_ifl = 6                          : number of blocks of block-diagonal hybridization
 // N_unique_fl = 4                    : number of columns read from the file
 // N_flavors  = 8                     : size of the actual matrix of hybridization N_flavors x N_flavors
 // ifl_dim[N_ifl] = [1,1,2,2,1,1]     : dimension of each block in block-diagonal hybridization
 // bfl_index[N_ifl][nbfl] = [[ 0],[ 0],[ 1, 3, 3, 1],[ 1, 3, 3, 1],[ 2], [2]]  : stores information how to construct matrix of hybridizations
 // sign[N_ifl][nbfl]      = [[ 1],[ 1],[ 1, 1, 1,-1],[ 1,-1,-1,-1],[ 1],[-1]]  : stores information how to construct matrix of hybridizations
 // conjg[N_ifl][nbfl]     = [[ 0],[ 1],[ 0, 0, 0, 1],[ 0, 0, 0, 1],[ 0],[ 1]]  : stores information how to construct matrix of hybridizations
 // gflip_index[N_ifl]     = [0,1,2,3,4,5]                                      : instructions for global flips
 //
 //
 // tfl_index[N_ifl][b1,b2] = (ifl=0 : [0]                        : index array for matrix of hybridizations
 //                           (ifl=1 : [0]
 //                           (ifl=2 : [[0,1],[2,3]]
 //                           (ifl=3 : [[0,1],[2,3]]
 //                           (ifl=4 : [0]
 //                           (ifl=5 : [0]
 //
 // vfl_index[N_ifl][b1,b2] = (ifl=0 : [[0]]
 //                           (ifl=1 : [[1]]
 //                           (ifl=2 : [[2,3],[4,5]]
 //                           (ifl=3 : [[6,7],[8,9]]
 //                           (ifl=4 : [[10]]
 //                           (ifl=5 : [[11]]
 //
 // v2fl_index[N_ifl][index] =  (ifl=0 : [0]
 //                             (ifl=1 : [1]
 //                             (ifl=2 : [2,3,4,5]
 //                             (ifl=3 : [6,7,8,9]
 //                             (ifl=4 : [10]
 //                             (ifl=5 : [11]
 //
 // Nvfl = 12
 // vfli_index[Nvfl] = [(0,0),
 //                     (1,0),
 //                     (2,0),
 //                     (2,1),
 //                     (2,2),
 //                     (2,3),
 //                     (3,0),
 //                     (3,1),
 //                     (3,2),
 //                     (3,3),
 //                     (4,0),
 //                     (5,0)]
 //
 //
 // fl_deg[N_unique_fl] = [2,4,2,4]
 //
 // ifl_from_fl[N_flavors] = (fl=0 : 0      bfl_from_fl[N_flavors] = (fl=0 : 0
 //                          (fl=1 : 1                               (fl=1 : 0
 //                          (fl=2 : 2                               (fl=2 : 0
 //                          (fl=3 : 2                               (fl=3 : 1
 //                          (fl=4 : 3                               (fl=4 : 0
 //                          (fl=5 : 3                               (fl=5 : 1
 //                          (fl=6 : 4                               (fl=6 : 0
 //                          (fl=7 : 5                               (fl=7 : 0
 //
 // fl_from_ifl[N_ifl][nbfl] =     (ifl=0 : [0]
 //                                (ifl=1 : [1]
 //                                (ifl=2 : [2, 3]
 //                                (ifl=3 : [4, 5]
 //                                (ifl=4 : [6]
 //                                (ifl=5 : [7]
 //
 // ifl_equivalent[N_ifl][N_ifl] = ( 1 1 0 0 0 0
 //                                ( 1 1 0 0 0 0
 //                                ( 0 0 1 1 0 0
 //                                ( 0 0 1 1 0 0
 //                                ( 0 0 0 0 1 1
 //                                ( 0 0 0 0 1 1
 //  
 //**************************************************************************************//
 //**************************************************************************************//
 // The second example if given for the f-shell in tetragonal crystal field where mixing
 // between Jz=5/2 and Jz=-3/2 occurs.
 //**************************************************************************************//
 // J   Jz  ifl   ifl_dim  bfl_index  gflip_index   fl_from_ifl
 // 5/2 -5/2, 3/2   2       0,7,7,1      0            0  4
 // 5/2  5/2,-3/2   2       0,7,7,1      0            5  1
 // 5/2  -1/2       1         2          1            2
 // 5/2   1/2       1         2          1            3
 // 7/2  -7/2       1         3          2            6
 // 7/2  -5/2       1         4          3            7
 // 7/2  -3/2       1         5          4            8
 // 7/2  -1/2       1         6          5            9
 // 7/2   1/2       1         6          5           10
 // 7/2   3/2       1         5          4           11
 // 7/2   5/2       1         4          3           12
 // 7/2   7/2       1         3          2           13
 // 
 // Note that hybridization is given in block form in different order than the electron operators.
 // The order of orbitals in the hybridization matrix is (-5/2,3/2,5/2,-3/2,-1/2,1/2,....)
 // while the order of electron operators is             (-5/2,-3/2,-1/2,1/2,3/5,5/2,....)
 // therefore fl_from_ifl index is not equal to identity.
 // 
 // fl_from_ifl= [[0,4],[5,1],[2],[3],[6],[7],[8],[9],[10],[11],[12],[13]]
 //               -5/2 -3/2 -1/2 1/2 3/2 5/2 
 // ifl_from_fl= [  0,   1,   2,  3,  0,  1,  4,  5,  6,  7,  8,  9,  10, 11]
 // bfl_from_fl= [  0,   1,   0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0]
 //
 // fl_deg = [2,2,2,2,2,2,2,4]
 // N_unique_fl = 8
 // N_ifl = 12
 // N_flavors = 14
 // tfl_index= [[[0,1],[2,3]],[[0,1],[2,3]],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0]]
 // vfl_index= [[[0,1],[2,3]],[[4,5],[6,7]],[8],[9],[10],[11],[12],[13],[14],[15],[16],[17]]
 // v2fl_index=[[0,1,2,3],[4,5,6,7],[8],[9],[10],[11],[12],[13],[14],[15],[16],[17]]
 // vfli_index=[(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3),(2,0),(3,0),(4,0),(5,0),(6,0),(7,0),(8,0),(9,0),(10,0),(11,0)]
 //*****************************************************************************//
  gout.ignore(1000,'\n');
  gout.ignore(1000,'\n');
  gout>>Nm>>nsize>>N_ifl>>max_size;

  if (!gout || Nm<0 || nsize<0 || N_ifl<0 || max_size<0){cerr<<"Something wrong reading cix file!"<<endl; exit(1);}
  gout.ignore(1000,'\n');
  gout.ignore(1000,'\n');

  ifl_dim.resize(N_ifl);
  sign.resize(N_ifl);
  conjg.resize(N_ifl);
  tfl_index.resize(N_ifl);
  vfl_index.resize(N_ifl);
  bfl_index.resize(N_ifl);
  gflip_index.resize(N_ifl);
  v2fl_index.resize(N_ifl);
  
  N_flavors=0;
  int max_dim=0;
  Nvfl=0;  
  for (int i=0; i<N_ifl; i++){
    int it;
    gout>>it>>ifl_dim[i];
    if (it!=i) {cerr<<"Something wrong in reading symmetries of baths in cix file!"<<endl; exit(1);}
    tfl_index[i].resize(ifl_dim[i],ifl_dim[i]);
    vfl_index[i].resize(ifl_dim[i],ifl_dim[i]);
    int ii=0;
    string str;
    for (int i1=0; i1<ifl_dim[i]; i1++){
      for (int i2=0; i2<ifl_dim[i]; i2++){
	int sgn=1; int cgn=0;
	gout>>str;
	size_t im = str.find("-");
	if (im != string::npos){
	  sgn=-1;
	  str.replace(im,1,"");
	}
	size_t ic = str.find("*");
	if (ic != string::npos){
	  cgn=1;
	  str.replace(ic,1,"");
	}
	int ind = atoi(str.c_str());
	bfl_index[i].push_back(ind);
	sign[i].push_back(sgn);
	conjg[i].push_back(cgn);
	Id[ind] = (i1==i2) ? 1 : 0;
	tfl_index[i][i1][i2] = ii++;
	
	vfl_index[i][i1][i2] = Nvfl;
	v2fl_index[i].push_back(Nvfl);
	Nvfl++;
      }
    }
    if (ifl_dim[i]>max_dim) max_dim = ifl_dim[i];
    N_flavors += ifl_dim[i];
    {
      int is;
      gout>>is;
      gflip_index[i]=is;
    }
    gout.ignore(1000,'\n');
  }
  gout.ignore(1000,'\n');


  N_unique_fl = 0;
  for (int ifl=0; ifl<N_ifl; ifl++)
    for (size_t b=0; b<bfl_index[ifl].size(); b++)
      if (bfl_index[ifl][b]>N_unique_fl) N_unique_fl = bfl_index[ifl][b];
  N_unique_fl++;
  
  fl_deg.resize(N_unique_fl);
  fl_deg=0;
  for (int ifl=0; ifl<N_ifl; ifl++)
    for (size_t b=0; b<bfl_index[ifl].size(); b++) fl_deg[bfl_index[ifl][b]]++;


  //// START NEW CHANGE 4.10.2008
  if (!gout){cerr<<"Something wrong in reading cix file! Exiting 0! "<<endl; exit(1);}
  string str0;
  getline(gout,str0);  // reading word "FL_FROM_IFL" or epsk[N_unique_fl] 
  // The order of psi_i operators is not necessary simply connected with the index in hybridization Delta.
  // The index fl_from_ifl connects electron operator psi_i with the index in hybridization matrix.
  fl_from_ifl.resize(N_ifl,max_dim);
  if (str0.find("FL_FROM_IFL",0)!=string::npos){
    deque<int> check_fl_from_ifl;
    int n_flavors=0;
    for (int ifl=0; ifl<N_ifl; ifl++){
      for (int ib=0; ib<ifl_dim[ifl]; ib++){
	gout>>n_flavors;
	if (!gout){cerr<<"Something wrong in reading IFL_FROM_FL in cix file! Exiting! "<<endl; exit(1);}
	fl_from_ifl(ifl,ib) = n_flavors;
	check_fl_from_ifl.push_back(n_flavors);
      }
      gout.ignore(1000,'\n');
    }
    sort(check_fl_from_ifl.begin(), check_fl_from_ifl.end());
    for (int i=0; i<check_fl_from_ifl.size(); i++){
      if (i!=check_fl_from_ifl[i]) {cerr<<"FL_FROM_IFL is wrong! Must be index array 0...N-1!"<<endl; exit(1);}
    }
    gout.ignore(1000,'\n');
    getline(gout,str0);
  }else{
    int n_flavors=0;
    for (int ifl=0; ifl<N_ifl; ifl++){
      for (int ib=0; ib<ifl_dim[ifl]; ib++){
	fl_from_ifl(ifl,ib) = n_flavors;
	n_flavors++;
      }
    }
  }
  if (!gout){cerr<<"Something wrong in reading cix file! Exiting! "<<endl; exit(1);}
  gout.ignore(1000,'\n'); // comment

  // reading from stringstream because we have epsk in string form
  std::istringstream sstr0(str0);
  epsk.resize(N_unique_fl);

  epsk=0;
  for (int ik=0; ik<N_unique_fl; ik++) sstr0>>epsk[ik];

  ifl_from_fl.resize(N_flavors);
  bfl_from_fl.resize(N_flavors);
  for (int ifl=0; ifl<N_ifl; ifl++){
    for (int jb=0; jb<ifl_dim[ifl]; jb++){
      int tfl = fl_from_ifl(ifl,jb);
      ifl_from_fl[tfl] = ifl;
      bfl_from_fl[tfl] = jb;
    }
  }

  //// END NEW CHANGE 4.10.2008
  
  vfli_index.resize(Nvfl);
  for (int i=0; i<N_ifl; i++){
    int ii=0;
    for (int i1=0; i1<ifl_dim[i]; i1++){
      for (int i2=0; i2<ifl_dim[i]; i2++){
	int jj = vfl_index[i][i1][i2];
	vfli_index[jj] = TBath(i,ii);
	ii++;
      }
    }
  }

  map<int,deque<int> > gfl_tmp;
  for (int ifl=0; ifl<N_ifl; ifl++){
    int ii = gflip_index[ifl];
    gfl_tmp[ii].push_back(ifl);
  }
  int sz=0;
  for (map<int,deque<int> >::iterator ia=gfl_tmp.begin(); ia!=gfl_tmp.end(); ia++){
    deque<int>& d = ia->second;
    sz += d.size()*(d.size()-1)/2;
  }
  
  gflip_fl.resize(sz+1,N_flavors);
  gflip_ifl.resize(sz);
  for (int iu=0; iu<sz; iu++){
    for (int j=0; j<gflip_fl.size_Nd(); j++) gflip_fl[iu][j] = j;
  }

  int iu=0;
  for (map<int,deque<int> >::iterator ia=gfl_tmp.begin(); ia!=gfl_tmp.end(); ia++){
    deque<int>& d = ia->second;
    for (size_t j1=0; j1<d.size(); j1++){
      for (size_t j2=j1+1; j2<d.size(); j2++){
	int ifa = d[j1];
	int ifb = d[j2];
	if (ifl_dim[ifa]!=ifl_dim[ifb]){cerr<<"Dimensions of similar baths have to be the same!"<<endl; exit(1);}
	gflip_ifl[iu] = make_pair(ifa,ifb);
	for (int ib=0; ib<ifl_dim[ifa]; ib++){
	  int fa = fl_from_ifl[ifa][ib];
	  int fb = fl_from_ifl[ifb][ib];
	  gflip_fl[iu][fa] = fb;
	  gflip_fl[iu][fb] = fa;
	  gflip_fl[sz][fa] = fb;
	  gflip_fl[sz][fb] = fa;
	}
	iu++;
      }
    }
  }
  ifl_equivalent.resize(N_ifl,N_ifl);
  for (int ifl=0; ifl<N_ifl; ifl++){
    for (int jfl=0; jfl<N_ifl; jfl++){
      ifl_equivalent[ifl][jfl] = bfl_index[ifl]==bfl_index[jfl];
    }
  }

  F_i.resize(2*N_flavors,nsize+1);
  F_M.resize(2*N_flavors,nsize+1);
  //Ene.resize(nsize+1,max_size);
  Ene0.resize(nsize+1,max_size);
  Spin.resize(nsize+1,max_size);
  Ns.resize(nsize+1);
  Ks.resize(nsize+1);
  Sz.resize(nsize+1);
  msize_.resize(nsize+1);
  
  msize_=0;
  F_i=0;
  
  int tsize=0;
  for (int i=1; i<=nsize; i++){
    int it1, in, ik, isize;   double dsz;
    gout>>it1>>in>>ik>>dsz>>isize;
    if (it1!=i){cerr<<"Something wrong in parsing cix file!"<<endl;}
    Ns[i] = in;
    Ks[i] = ik;
    Sz[i] = dsz;
    msize_[i] = isize;
    for (int ib=0; ib<N_flavors; ib++){
      int ist;
      gout>>ist;
      if (ist!=0){
	F_i(2*ib,i)=ist;   // first c^+
	F_i(2*ib+1,ist)=i; // next c
      }
    }
    for (int is=0; is<isize; is++) gout>>Ene0(i,is);
    Ene0[i].resize(isize);
    for (int is=0; is<isize; is++) gout>>Spin(i,is);
    Spin[i].resize(isize);
    gout.ignore(1000,'\n');
    
    tsize += isize;
  }

  // We will store energies in more compressed--efficient way for later computing.
  Eind.resize(nsize+1,max_size);
  Enes.resize(tsize);
  int iind=0;
  for (int i=1; i<=nsize; i++){
    int m=0;
    Eind[i].resize(msize_[i]);
    for (int j=0; j<msize_[i]; j++){
      Eind(i,j) = iind;
      Enes[iind] = Ene0(i,j);
      iind++;
    }
  }

  gout.ignore(1000,'\n');
  //Ene0 = Ene;
  
  for (int i=1; i<=nsize; i++){
    for (int ib=0; ib<N_flavors; ib++){
      int it, jt, size1, size2;
      gout>>it>>jt>>size1>>size2;
      if (it!=i || jt!=F_i(2*ib,i)) cerr<<"Something wrong reading cix file"<<endl;
      F_M(2*ib,i).resize(size2,size1); // constructor changes i->ii
      int ii = F_i(2*ib,i);
      F_M(2*ib+1,ii).resize(size1,size2);// destructor changes ii->i
      for (int i1=0; i1<size1; i1++){
	for (int i2=0; i2<size2; i2++){
	  double m;
	  gout>>m;
	  F_M(2*ib,i)(i2,i1)=m;    // constructor
	  F_M(2*ib+1,ii)(i1,i2)=m; // destructor
	}
      }
      gout.ignore(1000,'\n');
    }
  }

  string str;
  QHB1=false;
  getline(gout,str);
  //cout<<"str on HB="<<str<<endl;
  common::QHB2=false;
  if (str.find("HB1",0)!=string::npos){
    clog<<"Have HB1"<<endl;
    QHB1=true;
  }else if(str.find("HB2",0)!=string::npos){
    clog<<"Have HB2"<<endl;
    QHB1=true;
    common::QHB2=true;
  }else{
    // high frequency expansion, the first few moments of self-energy
    if (gout){
      getline(gout,str);
      RealSigma = str;
      getline(gout,str);
      ImagSigma = str;
    }
  }
  bool read_comment=true;
  if (common::QHB2){
    getline(gout,str); // # Uc = U[ifl,j1,j2,ifl]-U[ifl,j1,ifl,j2]:
    //cout<<"str on Uc="<<str<<endl;
    int a, b, c, d;
    double u;
    for (int i=0; i<100000; i++){
      getline(gout,str);
      stringstream ss(str); //covert input to a stream for conversion to int
      ss >> a >> b >> c >> d >> u;
      //cout<<"a,b,c,d="<<a<<" "<<b<<" "<<c<<" "<<d<<" u="<<u<<endl;
      if (!ss) break;
      UC.push_back(Uentry(a,b,c,d,u));
    }
    read_comment=false;
  }
  if (read_comment) getline(gout,str); // # number of operators needed

  //gout>>Osize;
  //gout.ignore(1000,'\n');
  string line;
  getline(gout,line);
  istringstream lines(line);
  lines >> Osize;
  DOsize=0;
  int _DOsize_=0;
  lines >> _DOsize_;
  if (lines){
    DOsize = _DOsize_;
  }

  clog<<"Number of operators needed "<<Osize<<"  "<<DOsize<<endl;
  getline(gout,str);
  HF_M.resize(Osize);
  
  for (int op=0; op<Osize; op++){// Operator for the high frequency moments, called Op in the string above
    if (!gout.good()) cerr<<"Something wrong in reading cix file high frequency expansion"<<endl;
    HF_M[op].resize(N_unique_fl,nsize+1);
    
    streampos pos1 = gout.tellg();
    int in=0;
    int im=0;
    int it;
    for (int ib=0; ib<N_unique_fl; ib++){
      gout>>it;
      gout.ignore(1000,'\n');
      if (it==1) in++;
      if (Id[ib]==1) im++;
    }
    bool off_diagonal=true;
    if (in==im){
      off_diagonal=false;
    }else if (in==N_unique_fl){
      off_diagonal=true;
    }else{
      cerr<<"Something wrong reading cix: Number of entries for Operator is not correct"<<endl;
    }
    clog<<"off_diagonal="<<off_diagonal<<endl;
    
    gout.seekg(pos1);
    for (int i=1; i<=nsize; i++){
      for (int ib=0; ib<N_unique_fl; ib++){
	int it=i;
	int size1 = msize_[i];
	int size2 = msize_[i];
	if (Id[ib]==1 || off_diagonal){  // Only diagonal baths need this. Otherwise just set it to zero!
	  gout>>it>>size1>>size2;
	  if (it!=i) cerr<<"Something wrong reading cix file reading Operators: i="<< i<<" it="<<it<<" size="<<size1<<endl;
	  if (size1!=msize_[i] || size2!=msize_[i]) cerr<<"Sizes of HFM are wrong! "<<i<<" "<<size1<<" "<<size2<<" "<<msize_[i]<<endl;
	}
	HF_M[op](ib,i).resize(size2,size1); // constructor changes i->ii
	for (int i1=0; i1<size1; i1++){
	  for (int i2=0; i2<size2; i2++){
	    double m=0.0;
	    if (Id[ib]==1 || off_diagonal) gout>>m;
	    HF_M[op](ib,i)(i2,i1)=m;
	  }
	}
	if (Id[ib]==1 || off_diagonal) gout.ignore(1000,'\n');
      }
    }	
    getline(gout,str);
  }

  //function2D<int> DF_i(DOsize,nsize+1);
  //function2D<function2D<double> > DF_M(DOsize,nsize+1,size2,size1);
  if (DOsize>0){
    DF_i.resize(DOsize,nsize+1);
    DF_inv.resize(DOsize,nsize+1);
    DF_i=0;
    DF_inv=0;
    DF_M.resize(DOsize,nsize+1);
    for (int op=0; op<DOsize; op++){ // Operators to be usef for dynamic susceptibility
      if (!gout.good()) cerr<<"Something wrong in reading cix file high frequency expansion"<<endl;
      for (int i=1; i<=nsize; i++){
	int it, jt, size1, size2;
	gout>>it>>jt;
	if (it!=i){
	  cerr<<"ERROR in cix file, dynamic operators, i="<<i<<" while it="<<it<<endl;
	  exit(1);
	}
	if (jt<0 || jt>nsize){
	  cerr<<"ERROR in cix file, dynamic operators, jt not in range. It is "<<jt<<" and should be between 1 and "<<nsize<<endl;
	  exit(1);
	}
	if (jt==0) {
	  DF_i(op,i)=0;
	  gout.ignore(1000,'\n');
	  continue;
	}
	gout>>size1>>size2;
	if (size1!=msize_[i]){
	  cerr<<"ERROR in cix file, dynamic operators, size1["<<i<<"]="<<size1<<" msize="<<msize_[i]<<endl;
	  exit(1);
	}
	if (size2!=msize_[jt]){
	  cerr<<"ERROR in cix file, dynamic operators, size2["<<i<<"]="<<size2<<" msize="<<msize_[jt]<<endl;
	  exit(1);
	}
	DF_i(op,i)=jt;
	DF_inv(op,jt)=i;
	DF_M(op,i).resize(size2,size1);
	for (int i1=0; i1<size1; i1++){
	  for (int i2=0; i2<size2; i2++){
	    double m;
	    gout>>m;
	    DF_M(op,i)(i2,i1)=m;
	  }
	}
	gout.ignore(1000,'\n');
      }
    }
    /*
    for (int i=1; i<=nsize; i++){
      if (DF_i(0,i)>0){
	cout<<i<<" "<<DF_i(0,i)<<" ";
	for (int i1=0; i1<DF_M(0,i).size_Nd(); i1++){
	  for (int i2=0; i2<DF_M(0,i).size_N(); i2++){
	    cout<<" "<<DF_M(0,i)(i2,i1)<<" ";
	  }
	}
	cout<<endl;
      }
    }
    */
  }
}

void ClusterData::EvaluateMa(double beta, double mu, double U)
{
  Ma.resize(N_flavors);
  for (size_t i=0; i<Ma.size(); i++) Ma[i].resize(size()+1,max_size);
  
  function2D<double> D(max_size,max_size);
  for (int i=1; i<=size(); i++){
    for (int op=0; op<N_flavors; op++){
      int inew = F_i(2*op+1,i);
      if (inew==0) continue;
      D.MProduct(F_M(2*op,inew),F_M(2*op+1,i));
      for (int m=0; m<D.size_N(); m++)	Ma[op](i,m) = D(m,m);
    }
  }

  for (int i=1; i<=nsize; i++){
    double dE = -Ns[i]*mu + 0.5*Ns[i]*(Ns[i]-1)*U;
    //for (int m=0; m<Ene0[i].size(); m++){
    for (int m=0; m<msize(i); m++)
      //Ene(i,m) += dE;
      Enes[Eind(i,m)] += dE;  // correcting for F0==U and the chemical potential
  }

  
  gs_ind.resize(nsize+1);
  for (int i=1; i<=nsize; i++){
    int imin = 0;
    /*
    for (int m=1; m<Ene[i].size(); m++)
      if (Ene(i,m)<Ene(i,imin)) imin = m;
    */
    for (int m=1; m<msize(i); m++)
      if (Enes[Eind(i,m)]<Enes[Eind(i,imin)]) imin = m;
    gs_ind[i] = imin;
    //    clog<<setw(3)<<i<<" "<<setw(4)<<gs_ind[i]<<endl;
  }
 
  Zatom=0;
  for (int i=1; i<=nsize; i++)
    for (int m=0; m<msize(i); m++)
      Zatom += Number(1,-Enes[Eind(i,m)]*beta);
  
  Patom.resize(size());
  Patom_.resize(size()+1,max_size);
  for (int i=1; i<=size(); i++){
    double sum=0;
    for (int m=0; m<msize(i); m++){
      Patom_(i,m) = divide(Number(1,-Enes[Eind(i,m)]*beta),Zatom);
      sum += Patom_(i,m);
    }
    Patom[i-1] = sum;
  }

  natom.resize(N_flavors);
  for (int b=0; b<N_flavors; b++){
    double sum=0;
    for (int i=1; i<=size(); i++)
      for (int m=0; m<msize(i); m++)
	sum += Patom_(i,m)*Ma[b](i,m);
    natom[b] = sum;
  }

}

bool CheckStream(istream& inputf, int& n, int& m)
{
  istream input(inputf.rdbuf());
  input.seekg(0,ios::beg);

  n=0;
  string str;
  bool begincomment = false;
  getline(input,str);
  if (!input.good()) {
    cerr << "ERROR: Wrong file format for input stream or no data!" << endl;
    return false;
  }
  if (str.find('#')<string::npos){
    begincomment = true;
  }else n++;
  
  getline(input,str); n++;
  stringstream oneline(str);
  m=0; double t;
  while (oneline){oneline>>t; m++;}
  m--;
  while (input){ getline(input,str); n++;}
  n--;

//   clog << " Number of entries: "<< n <<endl;
//   clog << " Number of columns: "<< m <<endl;

  inputf.seekg(0,ios::beg);
  if (begincomment) getline(inputf, str);
  if (!inputf){ cerr<<"Reopening didn't suceeded!"<<endl; return false;}
  return true;
}

bool ReadData(istream& input, mesh1D& om, function2D<dcomplex>& fi, int mdata, int baths)
{
  clog<<"Reading data with "<<om.size()<<" entries"<<endl;
  if (mdata<2*baths+1) {cerr<<"Not enough columns in input Delta file!"<<endl; return false;}
  vector<double> data(mdata);
  int i=-1;
  while (input && ++i<om.size()){
    for (int j=0; j<mdata; j++) input>>data[j];
    input.ignore(500,'\n');
    for (int j=0; j<baths; j++) fi[j][i] = dcomplex(data[2*j+1],data[2*j+2]);
    om[i] = data[0];
  }
  input.clear();              // forget we hit the end of file
  input.seekg(0, ios::beg);   // move to the start of the file
  string str;
  getline(input, str);        // read comment
  return true;
}

vector<pair<double,double> > FindHighFrequencyExp_original(int Nf, const mesh1D& omi, const function2D<dcomplex>& fi)
{
  vector<pair<double,double> > ah(fi.size_N());
  for (int b=0; b<fi.size_N(); b++){
    double S=0, Sx=0, Sy=0, Sxx=0, Sxy=0;
    double Sz=0, Sxz=0;
    for (int j=omi.size()-Nf; j<omi.size(); j++){
      double y = fi[b][j].imag()*omi[j];         // 1/omega    term
      double z = fi[b][j].real()*omi[j]*omi[j];  // 1/omega^2  term
      double x = omi[j];
      Sy += y;
      Sz += z;
      Sx += 1/(x*x);
      Sxx += 1/(x*x*x*x);
      Sxy += y/(x*x);
      Sxz += z/(x*x);
      S++;
    }
    double dd = S*Sxx-Sx*Sx;
    double ax = (Sxx*Sy-Sx*Sxy)/dd;
    double bx = (Sxx*Sz-Sx*Sxz)/dd;
    //    double bx = (S*Sxy-Sx*Sy)/dd;
    //    clog<<a<<" "<<bx<<endl;
    ah[b].first = -ax;
    ah[b].second = bx/ax;
  }
  return ah;
}

vector<pair<double,double> > FindHighFrequencyExp(int Nf, const mesh1D& omi, const function2D<dcomplex>& fi)
{// Looking for the form : a/(iom-eps). We will determine a and eps and store into ah.first=a and ah.second=eps
  //  cout.precision(16);
  vector<pair<double,double> > ah(fi.size_N());
  for (int b=0; b<fi.size_N(); b++){
    double S=0, Sy=0, Sz=0;
    for (int j=omi.size()-Nf; j<omi.size(); j++){
      double y = fi[b][j].imag()*omi[j];         // 1/omega    term
      double z = fi[b][j].real()*omi[j]*omi[j];  // 1/omega^2  term
      Sy += y;
      Sz += z;
      S++;
    }
    double ax = Sy/S;
    double bx = Sz/S;
    //cout<<"b="<<b<<" got ax="<<ax<<" got bx="<<bx<<endl;
    ah[b].first = -ax;    // a
    if (fabs(ax)>1e-15) ah[b].second = bx/ax; // eps
    else ah[b].second = 0.0;
  }
  return ah;
}

void CreateLogMesh(int M1, int M2, double T, const mesh1D& omi, mesh1D& oms)
{
  if (M1+M2+1>=omi.size()){
    oms = omi;
    return;
  }
  oms.resize(M1+M2+1);
  double small=0.1;
  for (int i=0; i<M1; i++) oms[i] = omi[i];

  double alpha = log((omi.size()-1.)/M1)/(M2-1.);
  
  int inp = static_cast<int>(0.5*(oms[M1-1]/(M_PI*T)-1)+small);
  int l=0;
  for (int i=0; i<M2; i++){
    int in = static_cast<int>(M1*exp(alpha*i)+0.5);
    if (in!=inp) oms[M1+(l++)] = (2*in+1)*M_PI*T;
    inp = in;
  }
  int last = static_cast<int>(0.5*(omi.last()/(M_PI*T)-1)+small);
  if (inp!=last) oms[M1+(l++)] = (2*last+1)*M_PI*T;
  oms.resize(M1+l);
  oms.SetUp(0.0);
}


void InverseFourier(const functionb<dcomplex>& Giom, const mesh& iom, functionb<double>& Gtau, const mesh& tau, const pair<double,double>& ah, double beta)
{
  function1D<dcomplex> dG(iom.size());
  for (int n=0; n<iom.size(); n++){
    dcomplex g_infty = ah.first/(dcomplex(0,iom[n])-ah.second);
    dG[n] = Giom[n]-g_infty;
  }
  
  for (int t=0; t<tau.size(); t++){
    double tau_ = tau[t];
    double dsum=0;
    for (int n=0; n<iom.size(); n++)
      dsum += cos(iom[n]*tau_)*dG[n].real()+sin(iom[n]*tau_)*dG[n].imag();
    
    Gtau[t] = 2*dsum/beta;
    if (ah.second>0)
      Gtau[t] -= ah.first*exp(-ah.second*tau[t])/(1.+exp(-ah.second*beta));
    else
      Gtau[t] -= ah.first*exp(ah.second*(beta-tau[t]))/(1.+exp(ah.second*beta));
  }
}

bool ReadDelta(int Ns, istream& input, int n, int m, mesh1D& omi, function2D<dcomplex>& Deltaw, vector<pair<double,double> >& ah, double beta, int Nhigh=40)
{
  mesh1D tomi(n);
  function2D<dcomplex> tDeltaw(Ns,tomi.size());
  if (!ReadData(input, tomi, tDeltaw, m, Ns)) return false;
  tomi.SetUp(0);
  int Nhighf = min(Nhigh, static_cast<int>(tomi.size()*0.4)); // correction 21.4.2008 : Shuld not compute high frequency from low frequency points
  ah = FindHighFrequencyExp(Nhighf, tomi, tDeltaw); 
  
  //int imax = static_cast<int>((tomi.last()*beta/M_PI-1)/2. + 0.5);
  int imax = static_cast<int>((tomi.last()*beta/M_PI+1)/2. + 1e-6); // correction, Jun 2015
  omi.resize(imax); // Matsubara mesh for this temperature
  for (int i=0; i<omi.size(); i++)
    omi[i] = (2*i+1)*M_PI/beta;

  // Is the mesh iom from file Delta.inp correct for this temperature?
  double mesh_diff=0;
  for (int i=0; i<min(omi.size(),tomi.size()); i++) mesh_diff += fabs(omi[i]-tomi[i]);
  //cout<<"mesh_diff="<<mesh_diff<<endl;
  
  Deltaw.resize(Ns,omi.size());
  for (int ib=0; ib<Ns; ib++){
    //cout<<"ah="<<ah[ib].first<<","<<ah[ib].second<<endl;
    if (mesh_diff>1.){ // The mesh of Matsubara points is not good enough
      spline1D<double> tDeltaw_r(tomi.size());
      spline1D<double> tDeltaw_i(tomi.size());
      for (int i=0; i<tomi.size(); i++){
	tDeltaw_r[i]=tDeltaw[ib][i].real();
	tDeltaw_i[i]=tDeltaw[ib][i].imag();
      }
      tDeltaw_r.splineIt(tomi,0,0);
      tDeltaw_i.splineIt(tomi,0,0);
      for (int i=0; i<omi.size(); i++)
	if (omi[i]<tomi.last()){
	  intpar pw = tomi.Interp(omi[i]);
	  Deltaw[ib][i] = dcomplex( tDeltaw_r(pw) , tDeltaw_i(pw) );
	} else{
	  Deltaw[ib][i] = ah[ib].first/( dcomplex(0,omi[i]) - ah[ib].second);
	}
    }else{
      for (int i=0; i<omi.size(); i++){
	if (omi[i]<tomi.last()+M_PI/beta) // correction 16.12.2006
	  Deltaw[ib][i] = tDeltaw[ib](tomi.Interp(omi[i]));
	else
	  Deltaw[ib][i] = ah[ib].first/( dcomplex(0,omi[i]) - ah[ib].second);
      }
    }
  }
  return true;
}

void DeltaFourier(int Ntau, double beta, mesh1D& tau, vector<spline1D<double> >& Delta, const mesh1D& omi, const function2D<dcomplex>& Deltaw, const vector<pair<double,double> >& ah)
{
  //tau.MakeEquidistantMesh(Ntau,0,beta);
  GiveDoubleExpMesh(tau,beta,Ntau);
  
  for (size_t ib=0; ib<Delta.size(); ib++) Delta[ib].resize(tau.size())
					     ;
  for (size_t ib=0; ib<Delta.size(); ib++){
    InverseFourier(Deltaw[ib], omi, Delta[ib], tau, ah[ib], beta); // Here needs to be corrected!
    //*** DEBUG
    //ofstream check(NameOfFile("check.dat",ib).c_str());
    //check<<"# ah = "<<ah[ib].first<<" "<<ah[ib].second<<endl; 
    //for (int i=0; i<tau.size(); i++) check<<tau[i]<<" "<<Delta[ib][i]<<endl;
    //***
    int n = tau.size();
    //double df0 = (Delta[ib][1]-Delta[ib][0])/(tau[1]-tau[0]);      // If Delta[ib][itau]>0 set it to something negative
    //double dfn = (Delta[ib][n-1]-Delta[ib][n-2])/(tau[n-1]-tau[n-2]);
    spline1D<double>& Dc = Delta[ib];
    double x1 = 0.5*(tau[1]+tau[0]);
    double df1 = (Dc[1]-Dc[0])/(tau[1]-tau[0]); // derivative in the first midpoint
    double x2 = 0.5*(tau[2]+tau[1]);
    double df2 = (Dc[2]-Dc[1])/(tau[2]-tau[1]); // derivative in the second midpoint
    double df0 = df1 + (df2-df1)*(0.0-x1)/(x2-x1); // extrapolated derivative at 0
    x1 = 0.5*(tau[n-1]+tau[n-2]);
    df1 = (Dc[n-1]-Dc[n-2])/(tau[n-1]-tau[n-2]); // derivative at the last mindpoint
    x2 = 0.5*(tau[n-2]+tau[n-3]);
    df2 = (Dc[n-2]-Dc[n-3])/(tau[n-2]-tau[n-3]); // derivative at the point before the last point
    double dfn = df1 + (df2-df1)*(common::beta-x1)/(x2-x1);  // extrapolated derivative
    
    Delta[ib].splineIt(tau,df0,dfn);
  }
}

// void ClusterData::Fill_fl_fl()
// {
//   for (int ii=0; ii<N_ifl; ii++) fl_fl.push_back(make_pair(ii,ii));
// }

/*
void ClusterData::Compute_Nmatrices()
{
  Njm.resize(size()+1,Nvfl);
  Njm_c.resize(size()+1,Nvfl);
  Njm_z.resize(size()+1);
  Njm_c = 0;
  Njm_z = 0;
  for (int ist=0; ist<size(); ist++){
    int istate = praState(ist);
    for (int ifl=0; ifl<N_ifl; ifl++){
      for (int bfls=0; bfls<ifl_dim[ifl]; bfls++){
	for (int bfle=0; bfle<ifl_dim[ifl]; bfle++){
	  int fls = fl_from_ifl[ifl][bfls];
	  int fle = fl_from_ifl[ifl][bfle];
	  int op_nodag = 2*fle+1;
	  int i_m_state = Fi(op_nodag)[istate];
	  if (i_m_state==0) continue;
	  int op_dagg = 2*fls;
	  int istate_new = Fi(op_dagg)[i_m_state];
	  if (istate_new!=istate) continue;
	  int fl2 = vfl_index[ifl](bfls,bfle);// Note that we are computing  psi^+_{bfls} psi_{bfle}
	  Njm(istate,fl2).resize(msize(istate),msize(istate));
	  function2D<double>& M = Njm(istate,fl2);
	  Multiply(M, FM(op_dagg)[i_m_state],FM(op_nodag)[istate]);
	  // check if it is zero
	  double dsum=0.0;
	  for (int i=0; i<M.size_N(); i++)
	    for (int j=0; j<M.size_Nd(); j++)
	      dsum += fabs(M(i,j));
	  if (dsum>1e-10) Njm_c(istate,fl2)=2; // means Njm is nonzero, but not idenity
	  dsum=0.0;
	  for (int i=0; i<M.size_N(); i++)
	    for (int j=0; j<M.size_Nd(); j++){
	      double c = i==j ? 1.0 : 0.0;
	      dsum += fabs(M(i,j)-c);
	    }
	  if (dsum<1e-10) Njm_c(istate,fl2)=1; // means Njm is equal to idenity
	  if (Njm_c(istate,fl2)==2) Njm_z[istate]=1; // at least one nontrivial N for this istate
	}
      }
    }
  }
  
  //cout<<"fl_fl="<<endl;
  //for (int ii=0; ii<fl_fl.size(); ii++){
  //  cout<<ii<<" "<<fl_fl[ii].first<<" "<<fl_fl[ii].second<<endl;
  //}
  //cout<<endl;
  //for (int istate=1; istate<=size(); istate++){
  //  cout<<"istate="<<istate<<endl;
  //  for (int fl=0; fl<Njs[istate].size(); fl++){
  //    cout<<"  fl="<<fl<<" : "<<Njs(istate,fl)<<endl;
  //  }
  //}
  //for (int istate=1; istate<=size(); istate++){
  //  for (int i=0; i<Njjs[istate].size(); i++){      
  //    cout<<istate<<" istate="<<setw(2)<<Njjs[istate][i].istate<<" (ifl1,ifl2)="<<setw(2)<<Njjs[istate][i].ifl1<<setw(2)<<Njjs[istate][i].ifl2<<" ii(fl1,fl2)="<<Njjs[istate][i].ii<<" i_m_state="<<setw(2)<<Njjs[istate][i].i_m_state<<" sizes="<<Njjs[istate][i].M.size_N()<<","<<Njjs[istate][i].M.size_Nd()<<endl;
  //    for (int j1=0; j1<Njjs[istate][i].M.size_N(); j1++){
  //    	for (int j2=0; j2<Njjs[istate][i].M.size_Nd(); j2++)
  //    	  cout<<Njjs[istate][i].M(j1,j2)<<" ";
  //    }
  //    cout<<endl;      
  //    int ii=Njjs[istate][i].ii;
  //    //cout<<Njjs[istate][i].ifl1<<"="<<fl_fl[ii].first<<"   "<<Njjs[istate][i].ifl2<<"="<<fl_fl[ii].second<<endl;
  //    if (Njjs[istate][i].ifl1!=fl_fl[ii].first || Njjs[istate][i].ifl2!=fl_fl[ii].second) cerr<<"ERROR: ifl1,ifl2="<<Njjs[istate][i].ifl1<<" "<<Njjs[istate][i].ifl2<<" while fl_fl gives ifl1,ifl2="<<fl_fl[ii].first<<" "<<fl_fl[ii].second<<endl;
  //  }
  //}
  //for (int istate=1; istate<=size(); istate++){
  //  for (int fl2=0; fl2<Nvfl; fl2++){
  //    function2D<double>& M = Njm(istate,fl2);
  //    cout<<" istate="<<setw(2)<<istate<<" fl2="<<setw(2)<<fl2<<" Njm_c="<<int(Njm_c(istate,fl2))<<" Njm_z="<<int(Njm_z[istate])<<" size="<<M.size_N()<<","<<M.size_Nd()<<endl;
  //    for (int j1=0; j1<M.size_N(); j1++){
  //    	for (int j2=0; j2<M.size_Nd(); j2++)
  //    	  cout<<M(j1,j2)<<" ";
  //    }
  //    cout<<endl;
  //  }
  //}
}
*/

void ClusterData::RecomputeCix(function2D<double>& AProb, double asign, double treshold, function1D<NState>& praStates)
{
  vector<deque<int> > redundant(nsize);
  function1D<int> new_msize_(msize_);
  //  deque<int> remove;
  for (int ist=0; ist<nsize; ist++){
    int ii = ist+1;
    deque<int> redun;
    for (int m=0; m<msize_[ii]; m++)
      if (fabs(AProb(ist,m))<treshold) redun.push_back(m); // this state should be eliminated
    redundant[ist] = redun;
    //    if (redun.size()==msize_[ist+1])  remove.push_back(ist);
    new_msize_[ii] = msize_[ii] - redun.size();
  }

  for (int ist=0; ist<nsize; ist++){
    int ii = ist+1;
    if (new_msize_[ii]==0) praStates[ist].SetIstate(0);
    
    for (int ib=0; ib<N_flavors; ib++){
      int jj = F_i(2*ib,ii);
      if (jj==0) continue;
      int size1 = msize_[ii];
      int size2 = msize_[jj];
      
      for (int k=redundant[ii-1].size()-1; k>=0; k--){
	int ik = redundant[ii-1][k];
	for (int m1=0; m1<size2; m1++)
	  for (int l=ik+1; l<size1; l++){
	    F_M(2*ib,ii)(m1,l-1) = F_M(2*ib,ii)(m1,l);
	    F_M(2*ib+1,jj)(l-1,m1) = F_M(2*ib+1,jj)(l,m1);
	  }
	size1--;
      }
      for (int k=redundant[jj-1].size()-1; k>=0; k--){
	int ik = redundant[jj-1][k];
	for (int m2=0; m2<size1; m2++)
	  for (int l=ik+1; l<size2; l++){
	    F_M(2*ib,ii)(l-1,m2) = F_M(2*ib,ii)(l,m2);
	    F_M(2*ib+1,jj)(m2,l-1) = F_M(2*ib+1,jj)(m2,l);
	  }
	size2--;
      }
      if (size1==0 || size2==0){
	F_i[2*ib][ii]=0;
	F_i[2*ib+1][jj]=0;
	size1=0;
	size2=0;
      }
      double sum=0;
      for (int m1=0; m1<size2; m1++)
	for (int m2=0; m2<size1; m2++)
	  sum += fabs(F_M(2*ib,ii)(m1,m2));
      if (sum<treshold){
	F_i[2*ib][ii]=0;
	F_i[2*ib+1][jj]=0;	
	size1=0;
	size2=0;
      }
      F_M(2*ib,ii).resize(size2,size1);
      F_M(2*ib+1,jj).resize(size1,size2);
    }
  }
  
  for (int ist=0; ist<nsize; ist++){
    int ii = ist+1;
    for (int k=redundant[ii-1].size()-1; k>=0; k--){
      int ik = redundant[ii-1][k];
      for (int l=ik+1; l<msize_[ii]; l++){
	Ene0[ii][l-1] = Ene0[ii][l];
	Spin[ii][l-1] = Spin[ii][l];
      }	
    }
    Ene0[ii].resize(new_msize_[ii]);
    Spin[ii].resize(new_msize_[ii]);
  }
  
  for (int op=0; op<Osize; op++){
    for (int ist=0; ist<nsize; ist++){
      int ii = ist+1;
      for (int ib=0; ib<N_unique_fl; ib++){
	function2D<double>& hf = HF_M[op](ib,ii);
	int size = msize_[ii];
	for (int k=redundant[ii-1].size()-1; k>=0; k--){
	  int ik = redundant[ii-1][k];
	  for (int m1=0; m1<size; m1++) for (int l=ik+1; l<size; l++) hf(m1,l-1) = hf(m1,l);
	  for (int m1=0; m1<size; m1++) for (int l=ik+1; l<size; l++) hf(l-1,m1) = hf(l,m1);
	}
	hf.resize(new_msize_[ii],new_msize_[ii]);
      }
    }
  }
  msize_ = new_msize_;
  
  ofstream gout("new.cix"); gout.precision(12);
  gout<<"# Updated Cix file for cluster DMFT with CTQMC"<<endl;
  gout<<"# cluster_size, number of states, number of baths, maximum_matrix_size"<<endl;
  gout<<Nm<<" "<<nsize<<" "<<N_ifl<<" "<<max_size<<endl;
  gout<<"# baths, dimension, symmetry"<<endl;
  for (int i=0; i<N_ifl; i++){
    gout<<left<<setw(3)<<i<<" "<<setw(4)<<ifl_dim[i]<<" ";
    for (int b=0; b<sqr(ifl_dim[i]); b++){
      if (sign[i][b]<0) gout<<"-";
      gout<<bfl_index[i][b];
      if (conjg[i][b]) gout<<"*";
      gout<<" ";
    }
    gout<<"   "<<setw(3)<<gflip_index[i]<<endl;
  }
  gout<<"# cluster energies for non-equivalent baths, eps[k]"<<endl;

  for (int ik=0; ik<epsk.size(); ik++) gout<<epsk[ik]<<" ";
  gout<<endl;
  
  gout<<right;
  gout<<" #   N   K   Sz size"<<endl;
  for (int ist=0; ist<nsize; ist++){
    int ii = ist+1;
    gout<<setw(2)<<ii<<" "<<setw(3)<<Ns[ii]<<" "<<setw(3)<<Ks[ii]<<" "<<setw(4)<<right<<Sz[ii]<<" "<<setw(3)<<msize_[ii]<<"    ";
    for (int ib=0; ib<N_flavors; ib++) gout<<setw(3)<<F_i(2*ib,ist+1)<<" ";
    //    gout<<" ";
    
    for (int is=0; is<Ene0[ii].size(); is++) gout<<setw(18)<<Ene0(ii,is)<<" ";
    for (int is=0; is<Spin[ii].size(); is++) gout<<setw(18)<<Spin(ii,is)<<" ";
    gout<<endl;
  }
  gout<<"# matrix elements"<<endl;
  for (int ist=0; ist<nsize; ist++){
    for (int ib=0; ib<N_flavors; ib++){
      gout<<setw(2)<<ist+1<<" "<<setw(3)<<F_i(2*ib,ist+1)<<"  ";
      function2D<double>& fm = F_M(2*ib,ist+1);
      gout<<setw(3)<<fm.size_Nd()<<" "<<setw(3)<<fm.size_N()<<"  ";
      for (int im1=0; im1<fm.size_Nd(); im1++)
	for (int im2=0; im2<fm.size_N(); im2++)
	  gout<<setw(20)<<fm[im2][im1]<<" ";
      gout<<endl;
    }
  }
  gout<<"# high frequency expansion (Real and Imaginary part of Sigma)"<<endl;
  gout<<RealSigma<<endl;
  gout<<ImagSigma<<endl;
  gout<<"# number of operators needed"<<endl;
  
  gout<<Osize<<endl;
  for (int op=0; op<Osize; op++){
    gout<<"# high frequency moment matrix elements, called Op"<<op+1<<" above"<<endl;
    for (int i=1; i<=nsize; i++){
      for (int ib=0; ib<N_unique_fl; ib++){
	int size = msize_[i];
	gout<<setw(2)<<i<<"  "<<setw(3)<<size<<" "<<setw(3)<<size<<"      ";
	function2D<double>& hf = HF_M[op](ib,i);
	for (int i1=0; i1<size; i1++)
	  for (int i2=0; i2<size; i2++)
	    gout<<setw(20)<<hf(i2,i1)<<" ";
	gout<<endl;
      }
    }
  }

  //Ene = Ene0;
  for (int i=1; i<=nsize; i++){
    double dE = -Ns[i]*common::mu + 0.5*Ns[i]*(Ns[i]-1)*common::U;
    for (int m=0; m<Ene0[i].size(); m++) Enes[Eind(i,m)] = Ene0(i,m) + dE;
  }
}

void ClusterData::Construct_Ucf(function4D<double>& Ucf, double U) const
{
  Ucf.resize(N_flavors,N_flavors,N_flavors,N_flavors);
  Ucf = 0.0;
  for (int k=0; k<UC.size(); k++){
    if (UC[k].a>=N_flavors || UC[k].b>=N_flavors || UC[k].c>=N_flavors || UC[k].d>=N_flavors) {cerr<<" In HB2::U component out of range..."<<endl; exit(1);}
    Ucf(UC[k].a,UC[k].b,UC[k].c,UC[k].d) = UC[k].u;
  }
  for (int a=0; a<N_flavors; a++)
    for (int b=0; b<N_flavors; b++)
      Ucf(a,b,b,a) += U;  // F0 contribution
}
/*
void ClusterData::Get_HB2_Uc(function2D<function2D<double> >& Uc, double U)
{
  function4D<double> Ucf(N_flavors,N_flavors,N_flavors,N_flavors);
  Construct_Ucf(Ucf, U);
  Uc.resize(N_ifl,N_ifl);
  for (int ifl=0; ifl<N_ifl; ifl++)
    for (int jfl=0; jfl<N_ifl; jfl++){
      Uc(ifl,jfl).resize(v2fl_index[ifl].size(), v2fl_index[jfl].size());
      Uc(ifl,jfl) = 0.0;
    }

  // Uc[ifl,jfl][(i,b),(j,k)] = U[(ifl,i),(jfl,j),(jfl,k),(ifl,b)] - U[(ifl,i),(jfl,j),(ifl,b),(jfl,k)]
  for (int ifl=0; ifl<N_ifl; ifl++){
    for (int i=0; i<ifl_dim[ifl]; i++){
      for (int b=0; b<ifl_dim[ifl]; b++){
	int ib = tfl_index[ifl](i,b);
	int a = fl_from_ifl[ifl][i];
	int d = fl_from_ifl[ifl][b];
	for (int jfl=0; jfl<N_ifl; jfl++){
	  for (int j=0; j<ifl_dim[jfl]; j++){
	    for (int k=0; k<ifl_dim[jfl]; k++){
	      int jk = tfl_index[jfl](j,k);
	      int b = fl_from_ifl[jfl][j];
	      int c = fl_from_ifl[jfl][k];
	      Uc(ifl,jfl)(ib,jk) = Ucf(a,b,c,d)-Ucf(a,b,d,c);
	    }
	  }
	}
      }
    }
  }
}
*/
void Inverse_Gf(int dim, int N_unique_fl, const deque<int>& bfl_index, const deque<int>& sign, const deque<int>& conjg, int im, const function2D<dcomplex>& Gf, function2D<dcomplex>& Gf_1)
{
  static function2D<dcomplex> g(dim,dim), gi(dim,dim);
  if (dim==1){
    Gf_1[bfl_index[0]][im] += 1/Gf[bfl_index[0]][im];
    return;
  }
  g.resize(dim,dim);
  gi.resize(dim,dim);
  int ib=0;
  for (int i1=0; i1<dim; i1++){
    for (int i2=0; i2<dim; i2++){
      dcomplex gf = Gf(bfl_index[ib],im);
      if (sign[ib]==-1) gf*=-1;
      if (conjg[ib]) gf = gf.conj();	  
      g[i1][i2] = gf;
      ib++;
    }
  }
  gi.Inverse(g);
  ib=0;
  for (int i1=0; i1<dim; i1++){
    for (int i2=0; i2<dim; i2++){
      dcomplex gf = gi[i1][i2];
      if (sign[ib]==-1) gf*=-1;
      if (conjg[ib]) gf = gf.conj();
      Gf_1(bfl_index[ib], im) += gf;
      ib++;
    }
  }
}

void F_times_Inverse_G(functionb<dcomplex>& Sg, const functionb<dcomplex>& Ff, const functionb<dcomplex>& Gf, const ClusterData& cluster)
{//(cluster.ifl_dim[ifl], cluster.N_unique_fl, cluster.bfl_index[ifl], cluster.sign[ifl], cluster.conjg[ifl], im, Gf, Gf_1);
  Sg = 0.0;
  for (int ifl=0; ifl<cluster.N_ifl; ifl++){
    int dim = cluster.ifl_dim[ifl];
    const deque<int>& _bfl_index_ = cluster.bfl_index[ifl];
    if (dim==1){
      int ii = cluster.bfl_index[ifl][0];
      Sg[ii] += Ff[ii]/Gf[ii];
      continue;
    }
    function2D<dcomplex> g(dim,dim), gi(dim,dim), fi(dim,dim), si(dim,dim);
    int ib=0;
    for (int i1=0; i1<dim; i1++){ // creates a matrix of G and F
      for (int i2=0; i2<dim; i2++){
	dcomplex gf = Gf[_bfl_index_[ib]];
	dcomplex ff = Ff[_bfl_index_[ib]];
	if (cluster.sign[ifl][ib]==-1){
	  gf *= -1;
	  ff *= -1;
	}
	if (cluster.conjg[ifl][ib]){
	  gf = gf.conj();
	  ff = ff.conj();
	}
	g[i1][i2] = gf;
	fi[i1][i2] = ff;
	ib++;
      }
    }
    gi.Inverse(g);      // This is just G^{-1}
    si.MProduct(fi,gi); // This gives  F * G^{-1}
    
    ib=0;
    for (int i1=0; i1<dim; i1++){ // Now compresses S into smaller array
      for (int i2=0; i2<dim; i2++){
	dcomplex sf = si[i1][i2];
	if (cluster.sign[ifl][ib]==-1) sf *= -1;
	if (cluster.conjg[ifl][ib]) sf = sf.conj();
	Sg[_bfl_index_[ib]] += sf;
	ib++;
      }
    }
  }
  for (int i=0; i<cluster.N_unique_fl; i++) Sg[i] /= cluster.fl_deg[i];
}

// F_sampled = Fsvd[fl_aj][fl_ik][:]
template<class T>
void Multiply_F_with_U(function2D<T>& Ft, const vector<function2D<T> >& F_sampled, double U, const ClusterData& cluster){
  function4D<double> Uc(cluster.N_flavors,cluster.N_flavors,cluster.N_flavors,cluster.N_flavors);
  cluster.Construct_Ucf(Uc, U);
  int Ft_sizeNd = F_sampled[0].size_Nd();
  Ft.resize(cluster.N_unique_fl,Ft_sizeNd);
  Ft = 0.0;
  function1D<T> tfc(Ft_sizeNd);
  // \psi_ia  \psi^+_ii \psi^+_ij \psi_ik
  for (int ia=0; ia<cluster.N_flavors; ia++){
    int ifla = cluster.ifl_from_fl[ia];
    int bea = cluster.bfl_from_fl[ia];
    for (int ib=0; ib<cluster.N_flavors; ib++){
      int iflb = cluster.ifl_from_fl[ib];
      if (iflb != ifla) continue;
      int bsb = cluster.bfl_from_fl[ib];
      int ab_ind = cluster.tfl_index[ifla](bea,bsb);
      tfc = 0.0;
      for (int ii=0; ii<cluster.N_flavors; ii++){
	int ifli = cluster.ifl_from_fl[ii];
	int bsi = cluster.bfl_from_fl[ii];
	for (int ij=0; ij<cluster.N_flavors; ij++){
	  int iflj = cluster.ifl_from_fl[ij];
	  int bsj = cluster.bfl_from_fl[ij];
	  for (int ik=0; ik<cluster.N_flavors; ik++){
	    int iflk = cluster.ifl_from_fl[ik];
	    int bek = cluster.bfl_from_fl[ik];
	    if (ifla==ifli && iflj==iflk){
	      // corresponds to  \psi_a \psi_i^+
	      int fl_ai = cluster.vfl_index[ifla](bea,bsi);  // combined index for (ia,ii)
	      // corresponds to  \psi_j^+ \psi_k
	      int fl_jk = cluster.vfl_index[iflk](bsj,bek);  // combined index for (ij,ik)
	      double UU = 0.5*(Uc(ii,ij,ik,ib)-Uc(ij,ii,ik,ib));
	      // <m|\psi_a \psi_i^+|n><n|\psi_j^+ \psi_k|m>
	      const funProxy<T>& _Fsvd_ = F_sampled[fl_ai][fl_jk];
	      for (int l=0; l<Ft_sizeNd; l++) tfc[l] += _Fsvd_[l] * UU;
	    }else if (ifla==iflj && ifli==iflk){
	      // corresponds to \psi_a \psi_j^+
	      int fl_aj = cluster.vfl_index[ifla](bea,bsj);  // combined index for (ia,ij)
	      // corresponds to \psi_i^+ \psi_k
	      int fl_ik = cluster.vfl_index[iflk](bsi,bek);  // combined index for (ii,ik)
	      double UU = 0.5*(Uc(ii,ij,ik,ib)-Uc(ij,ii,ik,ib));
	      const funProxy<T>& _Fsvd_ = F_sampled[fl_aj][fl_ik];
	      for (int l=0; l<Ft_sizeNd; l++) tfc[l] -= _Fsvd_[l] * UU;
	    }
	  }
	}
      }
      int i_unique = cluster.bfl_index[ifla][ab_ind];
      Ft[i_unique] += tfc;
    }
  }
}

void ClusterData::Compute_Nmatrices2(double U, bool Print=false){
  function4D<double> Uc(N_flavors,N_flavors,N_flavors,N_flavors);
  Construct_Ucf(Uc, U);
  
  NState state0(common::max_size,common::max_size), state1(common::max_size,common::max_size), state2(common::max_size,common::max_size);
  NState state3(common::max_size,common::max_size), st_fp(common::max_size,common::max_size);
  
  Njm.resize(size()+1,N_flavors);
  Njm_r.resize(size()+1,N_flavors);
  Njm_c.resize(size()+1,N_flavors);
  //Njm_z.resize(size()+1);
  Njm_c = 0;
  //Njm_z = 0;
  Njm_r = 0;
  
  // \psi^+_ii \psi^+_ij \psi_ik | state0>
  for (int ist=0; ist<nsize; ist++){
    state0.SetPraState(ist,*this);
    for (int ia=0; ia<N_flavors; ia++){
      int ifla = ifl_from_fl[ia];
      int bea = bfl_from_fl[ia];
      st_fp.apply(FM(2*ia), Fi(2*ia), state0); // psi^+_{ia}|state0>

      Njm(ist+1,ia).resize(st_fp.M.size_N(),st_fp.M.size_Nd());
      function2D<double>& _njm_ = Njm(ist+1,ia);
      _njm_=0;
      
      for (int ik=0; ik<N_flavors; ik++){
	int iflk = ifl_from_fl[ik];
	int bek = bfl_from_fl[ik];
	state1.apply(FM(2*ik+1), Fi(2*ik+1), state0);// \psi_{ik}|state0>
	if (state1.istate==0) continue;
	for (int ij=0; ij<N_flavors; ij++){
	  int iflj = ifl_from_fl[ij];
	  int bsj = bfl_from_fl[ij];
	  state2.apply(FM(2*ij), Fi(2*ij), state1);// \psi^+_{ij}|state1>
	  if (state2.istate==0) continue;
	  for (int ii=0; ii<N_flavors; ii++){
	    int ifli = ifl_from_fl[ii];
	    int bsi = bfl_from_fl[ii];
	    state3.apply(FM(2*ii), Fi(2*ii), state2); // psi^+_{ii}|state2>
	    if (state3.istate==0) continue;
	    double UU = 0.5*(Uc(ii,ij,ik,ia)-Uc(ij,ii,ik,ia));
	    if (UU==0) continue;
	    // U_{i,j,k,ia} psi_i^+ psi_j^+ psi_k psi_ia |m> is Coulomb repulsion
	    // Here we use
	    // Njm(state0,ia,:,:) = \sum_{i,j,k} U_{i,j,k,ia} psi_i^+ psi_j^+ psi_k |state0>
	    if (state3.istate != st_fp.istate){
	      cout<<"WARNING : mode S is likely not working with this cix-file. Check that G mode gives the same answer!"<<endl;
	      cout<<"WARN : Probably UCoulomb in cix file is wrong : ";
	      cout<<"ist="<<ist+1<<" ia="<<ia<<" k="<<ik<<" j="<<ij<<" i="<<ii<<" state3.istate="<<state3.istate<<" nstate.istate="<<st_fp.istate<<endl;
	      continue;
	    }
	    state3.M *= UU;
	    _njm_ += state3.M;
	  }
	}
      }

      ///// Checking what is the structure of the operator
      ///// It could be:
      //      - zero, in which case Njm_c(ist+1,ia)=0
      //      - proportional to F^+, in which case Njm_c(ist+1,ia)=1 and Njm_r(ist+1,ia) contains proportionality constant
      //      - arbitrary, in which ase Njm_c(ist+1,ia)=2
      //      if at least one of ia orbitals is non-trivial, Njm_z[ist+1]=1. Otherwise Njm_z[ist+1]=0.
      double dsum=0.0;
      for (int i=0; i<_njm_.size_N(); i++)
	for (int j=0; j<_njm_.size_Nd(); j++)
	  dsum += fabs(_njm_(i,j));

      //cout<<"First dsum="<<dsum<<endl;
      
      if (dsum>1e-10){
	double ratio=0;
	double max_fp=0;
	// first finds the ratio between Njm and F^+ by taking the ratio at the largest F^+ component
	for (int i=0; i<_njm_.size_N(); i++){
	  for (int j=0; j<_njm_.size_Nd(); j++){
	    if (fabs(st_fp.M(i,j)) > max_fp){
	      max_fp = fabs(st_fp.M(i,j));
	      ratio = _njm_(i,j)/st_fp.M(i,j);
	      Njm_r(ist+1,ia) = ratio;
	    }
	  }
	}
	Njm_c(ist+1,ia)=1; // means Njm is proportional to F^+ (at least we think so at the moment)
	for (int i=0; i<_njm_.size_N(); i++){
	  for (int j=0; j<_njm_.size_Nd(); j++){
	    if (fabs( _njm_(i,j) - st_fp.M(i,j) * ratio ) > 1e-5 ){
	      Njm_c(ist+1,ia)=2; // means Njm is not proportional to F^+ (we just figured that out)
	      goto loop_out2;
	    }
	  }
	}

	if (Njm_c(ist+1,ia)==1 && Njm_r(ist+1,ia)==0) {
	  cout<<"It should not happen "<<endl;
	}
      loop_out2:
	//if (Njm_c(ist+1,ia)==2) Njm_z[ist+1]=1; // at least one nontrivial N for this istate
	continue;
      }
      
      if (Print){
	cout<<"ist="<<ist+1<<" ia="<<ia<<" : "<<st_fp.M.size_N()<<"x"<<st_fp.M.size_Nd();
	cout<<" c="<<(int)Njm_c(ist+1,ia)<<" r="<<setw(8)<<Njm_r(ist+1,ia);
	// Ratio between this operator and F^+ should be close to : n(ts+0)*U
	// X[i,j] = _njm_[i,j]/st_fp.M[i,j]
	cout<<" M/F=";
	for (int i1=0; i1<st_fp.M.size_N(); i1++){
	  for (int i2=0; i2<st_fp.M.size_Nd(); i2++){
	    double X;
	    if (fabs(st_fp.M(i1,i2)) > 1e-16){ 
	      X = _njm_(i1,i2)/st_fp.M(i1,i2);
	    }else{
	      if (fabs(_njm_(i1,i2))>1e-4) cout<<"WARNING : It seems M/F is not possible to compute. You should likely change mode S to mode G"<<endl;
	      X = 0.0;
	    }
	    cout<<setw(8)<< X <<", ";
	  }
	  cout<<"; ";
	}
	cout<<" F=";
	for (int i1=0; i1<st_fp.M.size_N(); i1++){
	  for (int i2=0; i2<st_fp.M.size_Nd(); i2++){
	    cout<<setw(8)<<st_fp.M(i1,i2)<<", ";
	  }
	  cout<<"; ";
	}
	cout<<endl;
      }
    }
  }
}


#ifdef AS
template <>
int sigP<0>(double Ps){
  return Ps>0 ? 1 : 1;
}
#endif

