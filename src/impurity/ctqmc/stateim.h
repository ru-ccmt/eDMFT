// @Copyright 2007 Kristjan Haule
// 

NState& NState::operator=(const NState& C)
{
  //clog<<"Nstate "<<C.istate<<" is being assigned"<<endl;
  istate = C.istate;
  exponent = C.exponent;
  M.SetEqual(C.M);
  //M = C.M;
  return *this;
}
NState::NState(const NState& C)
{
  //clog<<"Copy constructor for Nstate called"<<endl;
  istate = C.istate;
  exponent = C.exponent;
  M.SetEqual(C.M);
}

inline void NState::SetEqual(const NState& m)
{
  istate = m.istate;
  exponent = m.exponent;
  M.SetEqual(m.M);
}

string NState::TotalSize() const
{
  stringstream str;
  str<<M.fullsize_N()<<" "<<M.fullsize_Nd();
  return str.str();
}

string NState::CurrentSize() const
{
  stringstream str;
  str<<M.size_N()<<" "<<M.size_Nd();
  return str.str();
}

//inline void NState::Evolve(const NState& C0, const ClusterData& cluster, const function2D<double>& exp_)
inline void NState::Evolve(const NState& C0, const ClusterData& cluster, const funProxy<double>& exp_)
{
  //t_evolve.start();
  int N = C0.M.size_N();
  int Nd = C0.M.size_Nd();
  M.resize(N,Nd);
  istate = C0.istate;
  if (N==1 && Nd==1){
    double maxi = fabs(C0.M[0][0]);
    exponent = log(maxi) + exp_[cluster.Eind(istate,0)] + C0.exponent;
    M(0,0) = C0.M(0,0)>0 ? 1. : -1.;
    return;
  }
  //if (exp_[istate].size()!=M.size_N()) cerr<<"Exponents are not of right size!"<<endl;
  if (cluster.msize(istate)!=M.size_N()) cerr<<"exponents and matrices are of wrong size"<<endl;
    
  const double* __restrict__ C0_M = C0.M.MemPt();
  int Nd0 = C0.M.lda();
  //const double* __restrict__ _exp_ = exp_.MemPt();
  //int Nde = exp_.lda();
  const double* __restrict__ _exp_ = &exp_[cluster.Eind(istate,0)];
  double* __restrict__ _M_ = M.MemPt();
  int Ndm = M.lda();

  /*
  cout<<"Evolving "<<istate<<" ";
  for (int i=0; i<N; i++) cout<<_exp_[i]<<" ";
  cout<<endl;
  */
  
  // finds biggest entry in the matrix of the state and extract new exponents from the biggest entry
  double max_exp=-100000;
  for (int i=0; i<N; i++){
    double maxi = fabs(C0_M[i*Nd0+0]); // maxi is the biggest entry in the row
    for (int j=1; j<Nd; j++)	
      if (fabs(C0_M[i*Nd0+j])>maxi) maxi = fabs(C0_M[i*Nd0+j]);
    //double expn = log(maxi) + _exp_[istate*Nde+i];// each row needs to be multiplied with different time evolution
    double expn = log(maxi) + _exp_[i];      // each row needs to be multiplied with different time evolution
    if (expn>max_exp) max_exp = expn;        // because atomic energies are different within the state
  }
  exponent = max_exp + C0.exponent; // this is the maximal exponent which is assigned to the time evolved state
  
  for (int i=0; i<N; i++){ // The entries in the matrix need to be updated accordingly
    //double expo = exp( _exp_[istate*Nde+i]-max_exp );
    double expo = exp( _exp_[i]-max_exp );
    for (int j=0; j<Nd; j++) _M_[i*Ndm+j] = C0_M[i*Nd0+j]*expo;
  }
  //t_evolve.stop();
}

inline void NState::apply(const functionb<function2D<double> >& FM, const functionb<int>& Fi, const NState& C)
{
  //t_apply.start();
  int orig_state = C.istate;
  int new_state = Fi[orig_state];
  istate = new_state;
  exponent = C.exponent;
  //    if (istate!=0) M.SMProduct(FM[orig_state],C.M);
  //M.resize_clear( FM[orig_state].size_N(), C.M.size_Nd());  // new change 2012!

  if (istate==0){
    M.resize(0,0);
    return;
  }
  
  M.resize( FM[orig_state].size_N(), C.M.size_Nd());
  Multiply(M, FM[orig_state],C.M);
  //t_apply.stop();
}


inline bool NState::empty() const
{
  if (istate==0) return true;
  for (int i=0; i<M.size_N(); i++)
    for (int j=0; j<M.size_Nd(); j++)
      if (fabs(M(i,j))>common::minF) return false;
  return true;
}

inline void NState::Print(ostream& out) const
{
  for (int i=0; i<M.size_N(); i++){
    for (int j=0; j<M.size_Nd(); j++){
      out<<setw(10)<<M(i,j)<<" ";
    }
    out<<endl;
  }
}

inline Number NState::TraceProject(const NState& s){
  if (s.istate!=istate) return 0.0;
  if (M.size_N()!=M.size_Nd()){cerr<<"Non quadratic state in TraceProject!"<<endl;}
  double sum=0;
  for (int i=0; i<M.size_N(); i++) sum += M(i,i);
  return Number(sum,exponent+s.exponent);
}

inline Number NState::Project_to_Pra(const NState& s, functionb<Number>& Proj){
  if (s.istate!=istate) {
    for (int i=0; i<M.size_N(); i++) Proj[i]=0;
    return 0.0;
  }
  if (M.size_N()!=M.size_Nd()){cerr<<"Non quadratic state in TraceProject!"<<endl;}
  double sum=0;
  for (int i=0; i<M.size_N(); i++){
    sum += M(i,i);
    Proj[i] = Number(M(i,i),exponent+s.exponent);
  }
  return Number(sum,exponent+s.exponent);
}

inline Number NState::ScalarProduct(const NState& s)
{
  if (s.istate!=istate) return 0.0;
  double sum=0;
  for (int i=0; i<M.size_N(); i++)
    for (int j=0; j<M.size_Nd(); j++)
      sum += M(i,j)*s.M(i,j);
  return Number(sum, exponent+s.exponent);
}

inline void NState::SetPraState(int i, const ClusterData& cluster)
{
  istate = cluster.praState(i);
  int msize = cluster.msize(istate);
  M.resize(msize,msize);
  M = 0;
  for (int l=0; l<msize; l++) M(l,l)=1;
  exponent = 0;
}

std::ostream& operator<< (std::ostream& stream, const NState& st)
{
  stream<<"istate="<<st.istate<<" exponent="<<st.exponent<<" M["<<st.M.size_N()<<","<<st.M.size_Nd()<<"]="<<endl;
  for (int i=0; i<st.M.size_N(); i++){
    for (int j=0; j<st.M.size_Nd(); j++)
      stream<<st.M(i,j)<<" ";
    stream<<endl;
  }
  return stream;
}
