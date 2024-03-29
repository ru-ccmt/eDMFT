// @Copyright 2007 Kristjan Haule
// 

void SortTimes(double t_start1, int fls1, double t_end1, int fle1, double t_start2, int fls2, double t_end2, int fle2, function1D<double>& tm, function1D<int>& op);

class NOperators{
public:
  const ClusterData& cluster;           // this is reference to the information about the cluster atomic states (to the object containing that information)
  int time_order_sign;                  // current sign due to time ordering of the operators
  function1D<double> time;              // stores time of all kinks
  function1D<int> index;                // Kinks are not stored in successive order due to efficiency reasons. Index is an array which sorts them according to their time.
  function1D<int> type;                 // is type 2*fl (for creation) or 2*fl+1 (for destruction) operator, where fl is the index of the operator type (for one band model either spin-up or spin-down)
  deque<int> Empty;                     // empty space in the container
  function2D<double> Exps;              // contains factors : -dt*E_m
  function1D<int> tr_index;             // index for trial moves (not yet accepted moves)
  // Kinks are stored in this class for computing the local trace. They are also stored in the intervals class.
  // However, there are more intervals - as many as the number of independent baths. Sometimes we need to link
  // a kink from this class with the same kink in the intervals class. The arrays below contain the information
  // which can be used to find the link in intervals provided we know the position of the same kink in operators
  function1D<IntervalIndex> p_ifl;      // link between kink in this class and kink in intervals
  int Empty_used;                       
  int maxi_, mini_;
public:
  NOperators(int Nmax, const ClusterData& cluster_);
  
  pair<int,int> Add_Cd_C(int fls, double t_start, int fle, double t_end, const IntervalIndex& ifls, const IntervalIndex& ifle);
  bool Try_Add_Cd_C_(int fls, double& t_start, int fle, double& t_end, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right);
  pair<int,int> Try_Add_Cd_C(int fls, double t_start, int fle, double t_end);
  
  void Remove_C_Cd(int ipe, int ips);
  bool Try_Remove_C_Cd_(int fle, int ie, int fls, int is, int& ipe, int& ips, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right, int istep);
  void Try_Remove_C_Cd(int ipe, int ips);
  
  void Move(int ip_old, int ip_new);
  bool Try_Move_(int opera, double t_old, int i_old, double& t_new, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right);
  void TryMove(int opera, double t_old, int i_old, int& ip_old, int& ipo, double t_new, int& ip_new);
    
  bool Try_Add_2_Cd_2_C_(int fls1, double t_start1, int fle1, double t_end1, int fls2, double t_start2, int fle2, double t_end2, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right);
  pair<int,int> Try_Add_2_Cd_2_C(int fls1, double t_start1, int fle1, double t_end1, int fls2, double t_start2, int fle2, double t_end2);
  bool Try_Remove_2_C_2_Cd_(int fle1, int ie1, int fls1, int is1, int fle2, int ie2, int fls2, int is2, int& ipe1, int& ips1, int& ipe2, int& ips2, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right);
  void Try_Remove_C_Cd(int ipe1, int ips1, int ipe2, int ips2);
  
  double t(int i) const { return time[index[i]];}
  int typ(int i) const {return type[index[i]];}
  const IntervalIndex& iifl(int i) const {return p_ifl[index[i]];}
  
  int typ_transpose(int i) const {int o = type[index[i]]; return (o%2==0) ? o+1 : o-1;}
  int sign() const {return time_order_sign;}
  bool empty() const {return size()==0;}
  double tr_t(int i) const {return time[tr_index[i]];  }
  int tr_typ(int i) const {return type[tr_index[i]];}
  int size() const {return index.size()-Empty.size();}
  bool full() const {return Empty.size()<=6;}
  int max_size() const {return index.size();}
  int mini(){return mini_;}
  int maxi(){return maxi_;}

  const funProxy<double>& exp_(int it) const{ int itt = it<size() ? index[it] : index.size(); return Exps[itt]; }
  const funProxy<double>& exp_last() const { return Exps[index.size()];}
  const funProxy<double>& tr_exp_(int it) const { return Exps[tr_index[it]];}
  const funProxy<double>& tr_exp_last() const { return Exps[index.size()+1];}

  void GlobalFlipChange(const functionb<int>& gflip_fl, int ifa, int ifb);
  void GlobalFlipChange(const functionb<int>& gflip_fl, const function1D<pair<int,int> >& gflip_ifl);
  void ExchangeTwo(int ip_a, int ip_b);

  int FindSuccessive(int opera, int which, int istart) const;
  int FindIndex(int opera, double tc, int lo, int istep) const;
  int FindIndexSlow(int opera, double tc) const;
  int FindExactTime(double tc, int lo) const;
  void Resize(int newNmax);
  void ReleaseSomeTempMemory();
private:
  void ChangeEntry(int isort, double time_, int type_, double dt, int new_place, function1D<int>& x_index, const IntervalIndex& p_ifl_);
  void Add(int op, double t, int& is, int& count, const IntervalIndex& iflx);
  int Remove(int to_remove);
  void TryAdd_basic(int which, int op, double t, int& is);
  void TryRemove_basic(int which, int to_remove);
};


inline NOperators::NOperators(int Nmax, const ClusterData& cluster_) :
  cluster(cluster_), time_order_sign(1), time(Nmax+2), index(Nmax), type(Nmax+2), Empty(Nmax),
  Exps(Nmax+2,cluster.Enes.size()), tr_index(Nmax), p_ifl(Nmax+2)
{
  for (int i=0; i<Nmax; i++) Empty[i]=i;
  /*
  for (int i=0; i<Expn.size(); i++){
    Expn[i].resize(cluster.nsize+1,cluster.max_size);
    for (int j=1; j<Expn[i].size_N(); j++) Expn[i][j].resize(cluster.msize(j));
  }
  for (int i=1; i<Expn[Nmax].size_N(); i++){
    for (int m=0; m<Expn[Nmax][i].size(); m++)
      Expn[Nmax](i,m) = -cluster.Ene(i,m)*common::beta;
  }
  */
  for (int m=0; m<cluster.Enes.size(); m++){
    Exps(Nmax,m) = -cluster.Enes[m]*common::beta;
    //Exps(0,m) = -cluster.Enes[m]*common::beta;   // because we start from the atomic state
  }
}

inline void NOperators::Resize(int newNmax)
{
  int Nmax = index.size();
  time.resize_preserve(newNmax+2);   // stores time of all kinks
  time.resize_virtual (newNmax+2);
  index.resize_preserve(newNmax);    // Kinks are not stored in successive order due to efficiency reasons. Index is an array which sorts them according to their time.
  index.resize_virtual (newNmax);
  type.resize_preserve(newNmax+2);   // is type 2*fl (for creation) or 2*fl+1 (for destruction) operator, where fl is the index of the operator type (for one band model either spin-up or spin-down)
  type.resize_virtual (newNmax+2);
  tr_index.resize_preserve(newNmax); // index for trial moves (not yet accepted moves)
  tr_index.resize_virtual (newNmax);
  p_ifl.resize_preserve(newNmax+2);     // link between kink in this class and kink in intervals
  p_ifl.resize_virtual (newNmax+2);
  for (int i=Nmax; i<newNmax; i++)
    Empty.push_back(i);              // empty space in the containers

  Exps.resize_preserve(newNmax+2, cluster.Enes.size());   // contains factors : -dt*E_m
  Exps.resize_virtual(newNmax+2, cluster.Enes.size());
  for (int m=0; m<cluster.Enes.size(); m++){// last entry should be moved to the end of the array by convention.
    Exps(newNmax,m)   = Exps(Nmax,m);    // the last time interval up to beta
    Exps(newNmax+1,m) = Exps(Nmax+1,m);// the same for the trial step
  }
  /*
  cout<<"sizes= "<<Exps.size_N()<<" "<<Exps.size_Nd()<<endl;
  for (int i=0; i<newNmax+2; i++){
    cout<<"isize="<<i<<" "<<Exps[i].size()<<endl;
    for (int j=0; j<cluster.Enes.size(); j++){
      cout<<i<<" "<<Exps[i][j]<<endl;
    }
  }
  */
}

void NOperators::ReleaseSomeTempMemory()
{
  time.ReleaseMemory();//~function1D();
  index.ReleaseMemory();//~function1D();
  type.ReleaseMemory();//~function1D();
  tr_index.ReleaseMemory();//~function1D();
  p_ifl.ReleaseMemory();//~function1D();
  Exps.ReleaseMemory();//~function2D();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Some auxiliary routines needed by the main routines below. They are located here          ///
/// because the compiler can inline them only if they are defined before they are used.       ///
/////////////////////////////////////////////////////////////////////////////////////////////////
inline int NOperators::FindSuccessive(int opera, int which, int istart=0) const
{
  int nsize = size();
  int i=0;
  int it=0;
  while (i<nsize){
    if (type[index[i]]==opera){
      it++;
      if (it==which+1) goto label;
    }
    i++;
  }
  if (i==nsize){cout<<"Could not find in FindSuccessive!! "<<nsize<<" "<<opera<<" "<<which<<" "<<istart<<endl;}
 label:
  return i;
}

inline int NOperators::FindIndex(int opera, double tc, int lo=0, int istep=-1) const
{
  const double small=common::beta*2e-16;
  // bisection
  int hi = size();
  while (lo < hi){
    int mid = (lo+hi)/2;
    double midval = time[index[mid]];
    if (midval < tc) lo = mid+1;
    else if (midval > tc) hi = mid;
    else{
      if (type[index[mid]]==opera) return mid;
      if (mid>0 && fabs(time[index[mid-1]]-tc)<small && type[index[mid-1]]==opera) return mid-1;
      if (mid+1<size() && fabs(time[index[mid+1]]-tc)<small && type[index[mid+1]]==opera) return mid+1;
      {
	cout<<"ERROR : Could not find time "<<tc<<" in FindIndex: t="<<tc<<" op="<<opera<<" at step istep="<<istep<<endl;
	cout<<"Currently have times : "<<endl;
	for (int i=0; i<size(); i++){
	  cout<<i<<" "<<t(i)<<" "<<typ(i)<<endl;
	}
      }
      return FindIndexSlow(opera,tc);
    }
  }
  if (lo<size() && fabs(time[index[lo]]-tc)<small) return lo;
  if (lo>0 && fabs(time[index[lo-1]]-tc)<small) return lo-1;
  if (lo+1<size() && fabs(time[index[lo+1]]-tc)<small) return lo+1;
  {
    cout<<"ERROR : Could not find time "<<tc<<" in FindIndex: t="<<tc<<" op="<<opera<<" at step istep="<<istep<<endl;
    cout<<"Currently have times : "<<endl;
    for (int i=0; i<size(); i++){
      cout<<i<<" "<<t(i)<<" "<<typ(i)<<endl;
    }
  }
  return FindIndexSlow(opera,tc);
}

inline int NOperators::FindIndexSlow(int opera, double tc) const
{
  const double small=common::beta*4e-16;
  ofstream eout("error.ctqmc");
  eout<<"myrank="<<common::my_rank<<endl;
  for (int i=0; i<size(); i++){
    eout<<i<<" typ="<<type[index[i]]<<" dt="<<time[index[i]]-tc<<endl;
    if (type[index[i]]==opera && fabs(time[index[i]]-tc)<small) return i;
  }
  eout<<"ERROR : Could not find time "<<tc<<" in FindIndexSlow: t="<<tc<<" op="<<opera<<endl;
  cout<<"ERROR : Could not find time "<<tc<<" in FindIndexSlow: t="<<tc<<" op="<<opera<<endl;
  return 0;
}

inline int NOperators::FindExactTime(double tc, int lo=0) const
{
  // bisection
  int hi = size();
  while (lo < hi){
    int mid = (lo+hi)/2;
    double midval = time[index[mid]];
    if (midval < tc) lo = mid+1;
    else if (midval > tc) hi = mid;
    else return mid;
  }
  return -1;
}

inline void NOperators::ChangeEntry(int isort, double time_, int type_, double dt, int new_place, function1D<int>& x_index, const IntervalIndex& p_ifl_)
{
  tr_index[isort] = new_place;
  time[new_place] = time_;
  type[new_place] = type_;
  p_ifl[new_place] = p_ifl_;
  /*
  for (int i=1; i<Expn[0].size_N(); i++){
    for (int m=0; m<Expn[0][i].size(); m++){
      Expn[new_place](i,m) = -cluster.Ene(i,m)*dt;
    }
  }
  */
  for (int m=0; m<cluster.Enes.size(); m++) Exps(new_place,m) = -cluster.Enes[m]*dt;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Routines to add kinks /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
inline pair<int,int> NOperators::Add_Cd_C(int fls, double t_start, int fle, double t_end, const IntervalIndex& ifls, const IntervalIndex& ifle)
{ // The step to add two kinks succeded -> actually insert the two kinks
  int is1, is2;       // place to insert the new kink
  int count1, count2; // size of the rest of the interval
  Add(2*fls, t_start, is1, count1, ifls);
  Add(2*fle+1, t_end, is2, count2, ifle);
  if (t_end<t_start) is1++;
  time_order_sign = 1-2*((count1+count2)%2);
  return make_pair(is1,is2);
}

inline bool NOperators::Try_Add_Cd_C_(int fls, double& t_start, int fle, double& t_end, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right)
{ // Only trial step // migh succeed or might not. This step should not modify the status of the system
  // This is very optimized routine which is called at every trial move when adding a kink.
  // With simple lookups in the table tries to determin if the matrix element is zero or not.
  double t0, t1; int op0, op1;
  if (t_start<t_end){
    t0 = t_start;
    op0 = 2*fls;
    t1 = t_end;
    op1 = 2*fle+1;
  }else{
    t1 = t_start;
    op1 = 2*fls;
    t0 = t_end;
    op0 = 2*fle+1;
  }
  int nsize = size();
  int i0=0;
  while (i0<nsize && time[index[i0]]<=t0) i0++;
  int i1=i0;
  while (i1<nsize && time[index[i1]]<=t1) i1++;

  if (i0>0 && time[index[i0-1]]==t0){
    // Happens very rearly that the new time is equal an old time. But it needs to be avoided
    if (t_start<t_end) t_start += common::smallest_dt;
    else t_end += common::smallest_dt;
    //    clog<<"Found a case of equal times!"<<endl;
  }
  if (i1>0 && time[index[i1-1]]==t1){
    // Happens very rearly that the new time is equal an old time. But it needs to be avoided
    if (t_start<t_end) t_end += common::smallest_dt;
    else t_start += common::smallest_dt;
    //    clog<<"Found a case of equal times!"<<endl;
  }
  bool succn;
  for (int i=0; i<cluster.size(); i++){
    int st0 = cluster.praState(i);
    succn = true;
    int stn=st0;
    if (i0>0){
      if (state_evolution_left[i].size()<i0) {succn=false; goto loops_out;}
      stn = state_evolution_left(i,i0-1).istate;
      if (stn==0) {succn=false; goto loops_out;}
    }
    if (i1<nsize && state_evolution_right[i].size()<=nsize-1-i1) {succn=false; goto loops_out;}
    stn = cluster.Fi(op0,stn);
    
    if (stn==0) {succn=false; goto loops_out;}
    for (int l=i0; l<i1; l++){
      int op = type[index[l]];
      stn = cluster.Fi(op,stn);
      if (stn==0) {succn=false; goto loops_out;}
    }
    stn = cluster.Fi(op1,stn);
    if (stn==0) {succn=false; goto loops_out;}
    
    if (i1<nsize && stn!=state_evolution_right[i][nsize-1-i1].istate) {succn=false; goto loops_out;}
    if (succn) break;
  loops_out:;
  }
  return succn;
}


inline pair<int,int> NOperators::Try_Add_Cd_C(int fls, double t_start, int fle, double t_end)
{ // Only trial step 
  // After the lookups were successful, the matrix element is probably nonzero and the move has a finite probability.
  // This is still a trial step and therefore should not modify the status of the system. However, it has to
  // give the matrix element of the trial "new" configuration.
  int is1, ie1;
  TryAdd_basic(0, 2*fls, t_start, is1);
  TryAdd_basic(1, 2*fle+1, t_end, ie1);
  if (t_end<t_start) is1++;
  return make_pair(is1,ie1);
}



////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Routines to remove kinks //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
inline void NOperators::Remove_C_Cd(int ipe, int ips)
{// The step to remove two kinks succeded -> actually remove the two kinks
  int count1 = Remove(ipe);
  int count2 = Remove(ips);
  time_order_sign = 1-2*((count1+count2)%2);
}

inline bool NOperators::Try_Remove_C_Cd_(int fle, int ie, int fls, int is, int& ipe, int& ips, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right, int istep)
{// Only trial step // migh succeed or might not. This step should not modify the status of the system
  // This is very optimized routine which is called at every trial move when removing a kink.
  // With simple lookups in the table tries to determin if the matrix element is zero or not.
   
  ipe = FindSuccessive(2*fle+1, ie);
  ips = FindSuccessive(2*fls,   is);
  int i0, i1;
  if (ips<ipe){
    i0 = ips;
    i1 = ipe;
  }else{
    i1 = ips;
    i0 = ipe;
  }
  int nsize = size();
  bool succn = true;
  int st0;

  for (int i=0; i<cluster.size(); i++){
    st0 = cluster.praState(i);
    succn=true;
    int stn=st0;
    if (i0>0){
      if (state_evolution_left[i].size()<i0) {succn=false; goto loops_out;}
      stn = state_evolution_left(i,i0-1).istate;
      if (stn==0) {succn=false; goto loops_out;}
    }
    if (i1<nsize-1 && state_evolution_right[i].size()<=nsize-2-i1) {succn=false; goto loops_out;}
    
    for (int l=i0+1; l<i1; l++){
      int op = type[index[l]];
      stn = cluster.Fi(op,stn);
      if (stn==0) {succn=false; goto loops_out;}
    }
    if (i1<nsize-1){
      if (stn!=state_evolution_right[i][nsize-2-i1].istate) {succn=false; goto loops_out;}
    }else{
      if (stn!=st0) {succn=false; goto loops_out;}
    }
    if (succn) break;
  loops_out:;
  }
  maxi_ = max(ips,ipe);
  if (ips>ipe) ips--;
  mini_ = min(ips,ipe);

  return succn;
}

inline void NOperators::Try_Remove_C_Cd(int ipe, int ips)
{ // Only trial step 
  // After the lookup was successful, the matrix element is probably nonzero and the move has a finite probability.
  // This is still a trial step and therefore should not modify the status of the system. However, it has to
  // give the matrix element of the trial "new" configuration.
  TryRemove_basic(0, ipe);
  TryRemove_basic(1, ips);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Routines to move a kinks //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

inline void NOperators::Move(int ip_old, int ip_new)
{ // Routine which actually moves the kink
  // after the step is accepted
  int nsize = size();
  for (int i=0; i<Empty_used; i++) Empty.pop_front();
  for (int i=0; i<nsize; i++){
    int ii = index[i]; bool used=false;
    for (int j=0; j<nsize; j++) // is this index still used?
      if (ii==tr_index[j]){ used=true; break;}
    if (!used) Empty.push_back(ii);
  }
  if (size()!=nsize) cerr<<"Did not correctly fill up Empty array!"<<endl;
  for (int i=0; i<nsize; i++) index[i] = tr_index[i];
  ChangeEntry(nsize+1,  common::beta,  -1,  common::beta-time[index[nsize-1]], index.size(), index, -1);
  time_order_sign = 1-2*((ip_new+ip_old)%2);
}

inline bool NOperators::Try_Move_(int opera, double t_old, int i_old, double& t_new, const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right)
{// Only trial step // migh succeed or might not. This step should not modify the status of the system
  // This is very optimized routine which is called at every trial to move a kink.
  // With simple lookups in the table tries to determin if the matrix element is zero or not.
  int ip_old = FindSuccessive(opera, i_old);
  if (time[index[ip_old]]!=t_old) cerr<<"TryMove_basic failed, did not find the right time "<<t_old<<" "<<time[index[ip_old]]<<endl;
  
  int nsize = size();
  int in=0;
  while (in<nsize && time[index[in]]<=t_new) in++;

  if (in>0 && time[index[in-1]]==t_new){// Change March 2013
    // Happens very rearly that the new time is equal an old time. But it needs to be avoided
    //cout<<"3:Found a case of equal times at in="<<in<<": Changing ";
    t_new += common::smallest_dt;
  }

  bool succn;
  if (t_new>t_old){
    for (int i=0; i<cluster.size(); i++){
      int st0 = cluster.praState(i);
      succn = true;
      int stn=st0;
      if (ip_old>0){
	if (state_evolution_left[i].size()<ip_old) {succn=false; goto loops_out1;}
	stn = state_evolution_left(i,ip_old-1).istate;
	if (stn==0) {succn=false; goto loops_out1;}
      }
      if (in<nsize && state_evolution_right[i].size()<=nsize-1-in) {succn=false; goto loops_out1;}
      for (int l=ip_old+1; l<in; l++){
	int op = type[index[l]];
	stn = cluster.Fi(op,stn);
	if (stn==0) {succn=false; goto loops_out1;}
      }
      stn = cluster.Fi(opera,stn);
      if (stn==0) {succn=false; goto loops_out1;}

      if (in<nsize && stn!=state_evolution_right[i][nsize-1-in].istate) {succn=false; goto loops_out1;}
      if (succn) break;
    loops_out1:;
    }
  }else{
    for (int i=0; i<cluster.size(); i++){
      int st0 = cluster.praState(i);
      succn = true;
      int stn=st0;
      if (in>0){
	if (state_evolution_left[i].size()<in) {succn=false; goto loops_out2;}
	stn = state_evolution_left(i,in-1).istate;
	if (stn==0) {succn=false; goto loops_out2;}
      }
      if (state_evolution_right[i].size()<nsize-1-ip_old) {succn=false; goto loops_out2;}
	
      stn = cluster.Fi(opera,stn);
      if (stn==0) {succn=false; goto loops_out2;}
      for (int l=in; l<ip_old; l++){
	int op = type[index[l]];
	stn = cluster.Fi(op,stn);
	if (stn==0) {succn=false; goto loops_out2;}
      }
      if (ip_old<nsize-1 && stn!=state_evolution_right[i][nsize-2-ip_old].istate) {succn=false; goto loops_out2;}
      if (succn) break;
    loops_out2:;
    }
  }
  return succn;
}

inline void NOperators::TryMove(int opera, double t_old, int i_old, int& ip_old, int& ipo, double t_new, int& ip_new)
{ // Only trial step 
  // After the lookup was successful, the matrix element is probably nonzero and the move has a finite probability.
  // This is still a trial step and therefore should not modify the status of the system. However, it has to
  // give the matrix element of the trial "new" configuration.
  ip_old = FindSuccessive(opera, i_old);
  if (time[index[ip_old]]!=t_old){
    cerr<<"TryMove_basic failed, did not find the right time "<<t_old<<" "<<time[index[ip_old]]<<endl;
  }
    
  Empty_used=0;
  int nsize = size();
  // prepare index array where the choose operator (ip_old) is missing
  for (int i=0; i<ip_old; i++) tr_index[i] = index[i];
  for (int i=ip_old+1; i<nsize; i++) tr_index[i-1] = index[i];
  // finds the new place for this operator according to new time t_new
  int is=0;
  while (is<nsize-1 && time[tr_index[is]]<t_new) is++;
  // makes place for the new operator to insert
  for (int i=nsize-1; i>is; i--) tr_index[i] = tr_index[i-1];
  // changes the entry for this particular moved operator
  double dt = is>0 ? t_new - time[tr_index[is-1]] : t_new;
  ChangeEntry(is, t_new, opera, dt, Empty[Empty_used++], tr_index, p_ifl[index[ip_old]]);
  // its neighbor needs to be changes as well. Makes a copy of time and type, and changes exponents
  if (is<nsize-1)
    ChangeEntry(is+1, time[tr_index[is+1]], type[tr_index[is+1]], time[tr_index[is+1]]-t_new, Empty[Empty_used++], tr_index, p_ifl[tr_index[is+1]]);
  ipo = ip_old;
  if (ip_old<nsize-1){      // the operator closest to the old place needs to be changed as well      
    ipo = ip_old>1 ? ip_old-2 : 0;
    while (!(time[tr_index[ipo]]==time[index[ip_old+1]] && type[tr_index[ipo]]==type[index[ip_old+1]])) ipo++;
    if (ipo!=is && ipo!=is+1){
      double dt = ipo>0 ? time[tr_index[ipo]]-time[tr_index[ipo-1]] : time[tr_index[ipo]];
      ChangeEntry(ipo,  time[tr_index[ipo]],  type[tr_index[ipo]], dt, Empty[Empty_used++], tr_index, p_ifl[tr_index[ipo]]);
    }
  }
  ChangeEntry(nsize+1,  common::beta,  -1,  common::beta-time[tr_index[nsize-1]], index.size()+1, tr_index, -1);

  ip_new = is;
}

inline void NOperators::ExchangeTwo(int ip_a, int ip_b)
{ // Routine which actually exchanges two kinks
  // We keep exactly the same time. We just change the type pf kink
  int op_a = type[index[ip_a]];
  int op_b = type[index[ip_b]];
  type[index[ip_a]] = op_b;
  type[index[ip_b]] = op_a;
  IntervalIndex ifla = p_ifl[index[ip_a]];
  IntervalIndex iflb = p_ifl[index[ip_b]];
  p_ifl[index[ip_a]]= iflb;
  p_ifl[index[ip_b]]= ifla;
}

////////////////////////////////////////////////////////////////////////////////////////
/// Lower level routines which peform some operation and are called by above routines //
////////////////////////////////////////////////////////////////////////////////////////

inline void NOperators::Add(int op, double t, int& is, int& count, const IntervalIndex& iflx)
{ 
  int nsize = size();
  int to_insert = Empty.front(); // first empty space in the array
  time[to_insert] = t;           // just insert into container in the first empy space
  type[to_insert] = op;          // will reorder index rather than data
  p_ifl[to_insert] = iflx;       // save the index to the same kink in the intervals
  is=0;                          // because index rearangment is cheap, data rearangmenet is expensive
  while (is<nsize && time[index[is]]<=t) is++; // these times are before the current time and are already ordered
  count = (nsize-is);            // found the place where the order will need to be changed
  for (int i=nsize; i>is; i--) index[i] = index[i-1]; // first make a room for new fellow
  index[is] = to_insert;         // and finally insert it in.
  Empty.pop_front();             // now mark that this place in the container is occupied
  
  bool last = is == nsize;       // we need now the kink before this kink and the kink after this kink
  double dt1 = !last ? time[index[is+1]]-time[index[is]] : common::beta-time[index[is]]; // if the kink is the last, next time is beta
  double dt0 = (is>0) ? time[index[is]]-time[index[is-1]] : time[index[is]]; // if kink is the first, previous time is zero
  // We need to modify two time differences: [tau_new,t_previous] and [tau_next,tau_new]
  // In most of the cases, we modify the content of Exp[t_new]==Exp[index[is]] and Exp[t_next]==Exp[index[is+1]]
  // If this new kink is the last kink, we need to modify the last time interval which is stored in
  // index.size()
  int istore = !last ? index[is+1] : index.size(); 

  /*
  for (int i=1; i<Expn[0].size_N(); i++){
    for (int m=0; m<Expn[0][i].size(); m++){
      Expn[index[is]](i,m) = -cluster.Ene(i,m)*dt0;
      Expn[istore](i,m) = -cluster.Ene(i,m)*dt1;
    }
  }
  */
  for (int m=0; m<cluster.Enes.size(); m++){
    Exps(index[is],m) = -cluster.Enes[m]*dt0;
    Exps(istore,m) = -cluster.Enes[m]*dt1;
  }
}

inline int NOperators::Remove(int to_remove)
{
  int nsize = size();
  Empty.push_back(index[to_remove]);
  for (int i=to_remove; i<nsize-1; i++) index[i] = index[i+1];
  int count = nsize-to_remove-1;

  bool last = to_remove==nsize-1;
  double t1 = !last ? time[index[to_remove]] : common::beta;
  double t0 = to_remove>0 ? time[index[to_remove-1]] : 0;
  int istore = !last ? index[to_remove] : index.size();

  /*
  for (int i=1; i<Expn[0].size_N(); i++){
    for (int m=0; m<Expn[0][i].size(); m++)
      Expn[istore](i,m) = -cluster.Ene(i,m)*(t1-t0);
  }
  */
  for (int m=0; m<cluster.Enes.size(); m++) Exps[istore][m] = -cluster.Enes[m]*(t1-t0);
  return count;
}

inline void NOperators::TryAdd_basic(int which, int op, double t, int& is)
{
  int nsize = size();
  nsize += which;
  int to_insert = Empty[2*which];

  time[to_insert] = t;
  type[to_insert] = op;
  is=0;
  if (which==0){
    while (is<nsize && time[index[is]]<=t) { tr_index[is]=index[is]; is++;}
    for (int i=nsize; i>is; i--) tr_index[i] = index[i-1];
  }else{
    while (is<nsize && time[tr_index[is]]<=t) is++;
    for (int i=nsize; i>is; i--) tr_index[i] = tr_index[i-1];
  }
  tr_index[is] = to_insert;
  int to_neigh = Empty[2*which+1];
  bool last = is == nsize;
  double dt1 = !last ? time[tr_index[is+1]]-time[tr_index[is]] : common::beta-time[tr_index[is]];
  double dt0 = (is>0) ? time[tr_index[is]]-time[tr_index[is-1]] : time[tr_index[is]];
  int istore = to_neigh;

  if (!last){
    time[to_neigh] = time[tr_index[is+1]];
    type[to_neigh] = type[tr_index[is+1]];
    tr_index[is+1] = to_neigh;
  }else{
    istore = index.size()+1;
  }
  /*
  for (int i=1; i<Expn[0].size_N(); i++){
    for (int m=0; m<Expn[0][i].size(); m++){
      Expn[tr_index[is]](i,m) = -cluster.Ene(i,m)*dt0;
      Expn[istore](i,m) = -cluster.Ene(i,m)*dt1;
    }
  }
  if (!last && which==0){
    for (int i=1; i<Expn[0].size_N(); i++){
      for (int m=0; m<Expn[0][i].size(); m++){
	Expn[index.size()+1](i,m) = Expn[index.size()](i,m);
      }
    }
  }
  */
  for (int m=0; m<cluster.Enes.size(); m++){
    Exps(tr_index[is],m) = -cluster.Enes[m]*dt0;
    Exps(istore,m)       = -cluster.Enes[m]*dt1;
  }
  if (!last && which==0){
    for (int m=0; m<cluster.Enes.size(); m++){
      Exps(index.size()+1,m) = Exps(index.size(),m); // the last dt is now in Nmax+1 for trial step
    }
  }
}

inline void NOperators::TryRemove_basic(int which, int to_remove)
{
  int nsize = size();
  nsize -= which;
  if (which==0){
    for (int i=0; i<to_remove; i++) tr_index[i] = index[i];
    for (int i=to_remove; i<nsize-1; i++) tr_index[i] = index[i+1];
  }else{
    for (int i=to_remove; i<nsize-1; i++) tr_index[i] = tr_index[i+1];
  }

  int to_insert = Empty[which];
    
  bool last = to_remove==nsize-1; // Should it be to_remove>nsize-1?
  double t1 = !last ? time[tr_index[to_remove]] : common::beta;
  double t0 = to_remove>0 ? time[tr_index[to_remove-1]] : 0;
  int istore = !last ? tr_index[to_remove] : tr_index.size()+1;
  if (!last){
    istore = to_insert;
    time[istore] = time[tr_index[to_remove]];
    type[istore] = type[tr_index[to_remove]];
    tr_index[to_remove] = to_insert;
  }
  /*
  for (int i=1; i<Expn[0].size_N(); i++){
    for (int m=0; m<Expn[0][i].size(); m++){
      Expn[istore](i,m) = -cluster.Ene(i,m)*(t1-t0);
    }
  }
  if (!last && which==0){
    for (int i=1; i<Expn[0].size_N(); i++){
      for (int m=0; m<Expn[0][i].size(); m++){
	Expn[index.size()+1](i,m) = Expn[index.size()](i,m);
      }
    }
  }
  */
  
  for (int m=0; m<cluster.Enes.size(); m++) Exps(istore,m) = -cluster.Enes[m]*(t1-t0);
  if (!last && which==0)
    for (int m=0; m<cluster.Enes.size(); m++) Exps(index.size()+1,m) = Exps(index.size(),m);
}



inline pair<int,int> NOperators::Try_Add_2_Cd_2_C(int fls1, double t_start1, int fle1, double t_end1, int fls2, double t_start2, int fle2, double t_end2)
{
  vector<int> isx(4);
  TryAdd_basic(0, 2*fls1, t_start1, isx[0]);
  TryAdd_basic(1, 2*fle1+1, t_end1, isx[1]);
  if (t_end1<t_start1) isx[0]++;
  TryAdd_basic(2, 2*fls2, t_start2, isx[2]);
  if (t_start2<t_start1) isx[0]++;
  if (t_start2<t_end1) isx[1]++;
  TryAdd_basic(3, 2*fle2+1, t_end2, isx[3]);
  if (t_end2<t_start1) isx[0]++;
  if (t_end2<t_end1) isx[1]++;
  if (t_end2<t_start2) isx[2]++;
  sort(isx.begin(), isx.end());
  return make_pair(isx[0],isx[3]);
}

/////////////////////////////////////////////////////////////////////////////////////
class Cmp{
  const function1D<double>& tm;
public:
  Cmp(const function1D<double>& tm_) : tm(tm_){};
  bool operator()(int a, int b)
  { return tm[a] < tm[b];}
};

void SortTimes(double t_start1, int fls1, double t_end1, int fle1, double t_start2, int fls2, double t_end2, int fle2,
	       function1D<double>& tm, function1D<int>& op)
{
  static function1D<double> tm_o(4);
  static function1D<int> op_o(4);
  static function1D<int> index(4);
  tm.resize(4); op.resize(4);
  
  for (int i=0; i<index.size(); i++) index[i]=i;
  tm_o[0] = t_start1; op_o[0] = fls1;
  tm_o[1] = t_end1;   op_o[1] = fle1;
  tm_o[2] = t_start2; op_o[2] = fls2;
  tm_o[3] = t_end2;   op_o[3] = fle2;
  
  Cmp cmp(tm_o);
  sort(index.begin(),index.end(), cmp);
  for (int i=0; i<index.size(); i++){
    tm[i] = tm_o[index[i]];
    op[i] = op_o[index[i]];
  }
  
}

bool NOperators::Try_Add_2_Cd_2_C_(int fls1, double t_start1, int fle1, double t_end1, int fls2, double t_start2, int fle2, double t_end2,
				   const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right)
{
  static function1D<double> tm(4);
  static function1D<int> op(4);
  SortTimes(t_start1, 2*fls1, t_end1, 2*fle1+1, t_start2, 2*fls2, t_end2, 2*fle2+1, tm, op);
    
  int nsize = size();
  static function1D<int> isx(4);
  for (int l=0; l<4; l++){
    int ii = l>0 ? isx[l-1] : 0;
    while (ii<nsize && time[index[ii]]<=tm[l]) ii++;
    isx[l] = ii;
  }
    
  bool succn;
  for (int i=0; i<cluster.size(); i++){
    int st0 = cluster.praState(i);
    succn = true;
    int stn=st0;
    if (isx[0]>0){
      if (state_evolution_left[i].size()<isx[0]) {succn=false; goto loops_out;}
      stn = state_evolution_left(i,isx[0]-1).istate;
      if (stn==0) {succn=false; goto loops_out;}
    }
    if (isx[3]<nsize && state_evolution_right[i].size()<=nsize-1-isx[3]) {succn=false; goto loops_out;}
	
    stn = cluster.Fi(op[0],stn);
    if (stn==0) {succn=false; goto loops_out;}
    for (int l=isx[0]; l<isx[1]; l++){
      int op = type[index[l]];
      stn = cluster.Fi(op,stn);
      if (stn==0) {succn=false; goto loops_out;}
    }
    stn = cluster.Fi(op[1],stn);
    if (stn==0) {succn=false; goto loops_out;}
    for (int l=isx[1]; l<isx[2]; l++){
      int op = type[index[l]];
      stn = cluster.Fi(op,stn);
      if (stn==0) {succn=false; goto loops_out;}
    }
    stn = cluster.Fi(op[2],stn);
    if (stn==0) {succn=false; goto loops_out;}
    for (int l=isx[2]; l<isx[3]; l++){
      int op = type[index[l]];
      stn = cluster.Fi(op,stn);
      if (stn==0) {succn=false; goto loops_out;}
    }
    stn = cluster.Fi(op[3],stn);
    if (stn==0) {succn=false; goto loops_out;}

    if (isx[3]<nsize && stn!=state_evolution_right[i][nsize-1-isx[3]].istate) {succn=false; goto loops_out;}
    if (succn) break;
  loops_out:;
  }
  return succn;
}

bool NOperators::Try_Remove_2_C_2_Cd_(int fle1, int ie1, int fls1, int is1, int fle2, int ie2, int fls2, int is2,
				      int& ipe1, int& ips1, int& ipe2, int& ips2,
				      const function2D<NState>& state_evolution_left, const function2D<NState>& state_evolution_right)
{
  ipe1 = FindSuccessive(2*fle1+1, ie1);
  ips1 = FindSuccessive(2*fls1,   is1);
  ipe2 = FindSuccessive(2*fle2+1, ie2);
  ips2 = FindSuccessive(2*fls2,   is2);
  static function1D<int> it(4);
  it[0]=ipe1; it[1]=ips1; it[2]=ipe2; it[3]=ips2;
  sort(it.begin(),it.end());
  
  int nsize = size();
  bool succn = true;
  int st0;
  for (int i=0; i<cluster.size(); i++){
    st0 = cluster.praState(i);
    succn=true;
    int stn=st0;
    if (it[0]>0){
      if (state_evolution_left[i].size()<it[0]) {succn=false; goto loops_out;}
      stn = state_evolution_left(i,it[0]-1).istate;
      if (stn==0) {succn=false; goto loops_out;}
    }
    if (it[3]<nsize-1 && state_evolution_right[i].size()<=nsize-2-it[3]) {succn=false; goto loops_out;}
	
    for (int l=it[0]+1; l<it[1]; l++){
      int op = type[index[l]];
      stn = cluster.Fi(op,stn);
      if (stn==0) {succn=false; goto loops_out;}
    }
    for (int l=it[1]+1; l<it[2]; l++){
      int op = type[index[l]];
      stn = cluster.Fi(op,stn);
      if (stn==0) {succn=false; goto loops_out;}
    }
    for (int l=it[2]+1; l<it[3]; l++){
      int op = type[index[l]];
      stn = cluster.Fi(op,stn);
      if (stn==0) {succn=false; goto loops_out;}
    }
    if (it[3]<nsize-1){
      if (stn!=state_evolution_right[i][nsize-2-it[3]].istate) {succn=false; goto loops_out;}
    }else{
      if (stn!=st0) {succn=false; goto loops_out;}
    }
    if (succn) break;
  loops_out:;
  }
  maxi_ = max(max(ips1,ipe1),max(ips2,ipe2));
  if (ips1>ipe1) ips1--;
  if (ipe2>ipe1) ipe2--;
  if (ips2>ipe1) ips2--;
  if (ipe2>ips1) ipe2--;
  if (ips2>ips1) ips2--;
  if (ips2>ipe2) ips2--;
  mini_ = min(min(ips1,ipe1),min(ips2,ipe2));
  return succn;
}

// bool NOperators::Slow_Try_Add_2_Cd_2_C_(int fls1, double t_start1, int fle1, double t_end1, int fls2, double t_start2, int fle2, double t_end2)
// {
//   static function1D<double> tm(4);
//   static function1D<int> op(4);
//   SortTimes(t_start1, 2*fls1, t_end1, 2*fle1+1, t_start2, 2*fls2, t_end2, 2*fle2+1, tm, op);
    
//   int nsize = size();
//   static function1D<int> isx(4);
//   for (int l=0; l<4; l++){
//     int ii = l>0 ? isx[l-1] : 0;
//     while (ii<nsize && time[index[ii]]<=tm[l]) ii++;
//     isx[l] = ii;
//   }
    
//   bool succn;
//   for (int i=0; i<cluster.size(); i++){
//     int st0 = cluster.praState(i);
//     succn = true;
//     int stn=st0;
//     for (int l=0; l<isx[0]; l++){
//       int op = type[index[l]];
//       stn = cluster.Fi(op,stn);
//       if (stn==0) {succn=false; goto loops_out;}
//     }
//     stn = cluster.Fi(op[0],stn);
//     if (stn==0) {succn=false; goto loops_out;}
//     for (int l=isx[0]; l<isx[1]; l++){
//       int op = type[index[l]];
//       stn = cluster.Fi(op,stn);
//       if (stn==0) {succn=false; goto loops_out;}
//     }
//     stn = cluster.Fi(op[1],stn);
//     if (stn==0) {succn=false; goto loops_out;}
//     for (int l=isx[1]; l<isx[2]; l++){
//       int op = type[index[l]];
//       stn = cluster.Fi(op,stn);
//       if (stn==0) {succn=false; goto loops_out;}
//     }
//     stn = cluster.Fi(op[2],stn);
//     if (stn==0) {succn=false; goto loops_out;}
//     for (int l=isx[2]; l<isx[3]; l++){
//       int op = type[index[l]];
//       stn = cluster.Fi(op,stn);
//       if (stn==0) {succn=false; goto loops_out;}
//     }
//     stn = cluster.Fi(op[3],stn);
//     if (stn==0) {succn=false; goto loops_out;}
//     for (int l=isx[3]; l<nsize; l++){
//       int op = type[index[l]];
//       stn = cluster.Fi(op,stn);
//       if (stn==0) {succn=false; goto loops_out;}
//     }
//     if (stn!=st0) {succn=false; goto loops_out;}
//     if (succn) break;
//   loops_out:;
//   }
//   return succn;
// }

// bool NOperators::Slow_Try_Remove_2_C_2_Cd_(int fle1, int ie1, int fls1, int is1, int fle2, int ie2, int fls2, int is2)
// {
//   static function1D<int> it(4);
//   it[0] = FindSuccessive(2*fle1+1, ie1);
//   it[1] = FindSuccessive(2*fls1,   is1);
//   it[2] = FindSuccessive(2*fle2+1, ie2);
//   it[3] = FindSuccessive(2*fls2,   is2);
//   sort(it.begin(),it.end());
  
//   int nsize = size();
//   bool succn = true;
//   for (int i=0; i<cluster.size(); i++){
//     int st0 = cluster.praState(i);
//     succn=true;
//     int stn=st0;
//     for (int l=0; l<it[0]; l++){
//       int op = type[index[l]];
//       stn = cluster.Fi(op,stn);
//       if (stn==0) {succn=false; goto loops_out;}
//     }
//     for (int l=it[0]+1; l<it[1]; l++){
//       int op = type[index[l]];
//       stn = cluster.Fi(op,stn);
//       if (stn==0) {succn=false; goto loops_out;}
//     }
//     for (int l=it[1]+1; l<it[2]; l++){
//       int op = type[index[l]];
//       stn = cluster.Fi(op,stn);
//       if (stn==0) {succn=false; goto loops_out;}
//     }
//     for (int l=it[2]+1; l<it[3]; l++){
//       int op = type[index[l]];
//       stn = cluster.Fi(op,stn);
//       if (stn==0) {succn=false; goto loops_out;}
//     }
//     for (int l=it[3]+1; l<nsize; l++){
//       int op = type[index[l]];
//       stn = cluster.Fi(op,stn);
//       if (stn==0) {succn=false; goto loops_out;}
//     }
//     if (stn!=st0) {succn=false; goto loops_out;}
//     if (succn) break;
//   loops_out:;
//   }
//   return succn;
// }

void NOperators::Try_Remove_C_Cd(int ipe1, int ips1, int ipe2, int ips2)
{
  TryRemove_basic(0, ipe1);
  TryRemove_basic(1, ips1);
  TryRemove_basic(2, ipe2);
  TryRemove_basic(3, ips2);
}

void NOperators::GlobalFlipChange(const functionb<int>& gflip_fl, int ifa, int ifb)
{
  for (int i=0; i<size(); i++){
    int iop = type[index[i]];
    int nfl = gflip_fl[iop/2];
    int nop = 2*nfl + iop%2;
    type[index[i]] = nop;

    if(p_ifl[index[i]].ifl == ifa){
      p_ifl[index[i]].ifl = ifb;
    } else if (p_ifl[index[i]].ifl == ifb){
      p_ifl[index[i]].ifl = ifa;
    }
  }
}

void NOperators::GlobalFlipChange(const functionb<int>& gflip_fl, const function1D<pair<int,int> >& gflip_ifl)
{
  for (int i=0; i<size(); i++){
    int iop = type[index[i]];
    int nfl = gflip_fl[iop/2];
    int nop = 2*nfl + iop%2;
    type[index[i]] = nop;

    for (int iu=0; iu<gflip_ifl.size(); iu++){
      int ifa = gflip_ifl[iu].first;
      int ifb = gflip_ifl[iu].second;
      if(p_ifl[index[i]].ifl == ifa){
	p_ifl[index[i]].ifl = ifb;
      } else if (p_ifl[index[i]].ifl == ifb){
	p_ifl[index[i]].ifl = ifa;
      }
    }
  }
}
