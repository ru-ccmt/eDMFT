#ifndef _TOPOLOGY
#define _TOPOLOGY


int PPconfigs2(const Array<int,2>& diag, const Array<int,1>& dirspin, const Cix& cix){
  // Pauli exclusion principle checking given a diagram and a dirspin.
  // Using lookup method with table Fi
  //
  int n = diag.extent(0);  // Perturbation order
  Array<int,1> op(2*n);
  for (int ih = 0; ih < n; ++ih) {
    int i = diag(ih,0);
    int j = diag(ih,1);
    int ds = dirspin(ih);
    op(j) = ds;
    op(i) = 1-(ds & 1) + 2*(ds >> 1);  // operator opposite to ds, just like in dirspin2ops. Works for diagonal hybridization.
  }
  //cout<<"op="<< op<<endl;
  int PPs=0;
  for (int ib=0; ib<cix.Nblock; ++ib){
    int pp = ib;
    for (int i=0; i<2*n; ++i){
      pp = cix.Fbc(op(i),pp);
      if (pp==-1) break;
    }
    PPs += (pp!=-1); // Equivalent to incrementing if pp != 0
  }
  return PPs;
}

void dirspin2ops(Array<int,1>& op,
		 const Array<int,2>& diag, const Array<int,1>& dirspin){
  // Converts dirspin into operators Psi_{s,orbital} or Psi_^\dagger_{s,orbital}
  // Note that we are currently restricting to diagonal hybridization only, hence dirspin(n), which gives unique operators.
  int n = diag.extent(0);  // Perturbation order
  op.resize(2*n);  // Operators
  // Compute occupation states
  for (int ih = 0; ih<n; ++ih){ // over all hybridizations
    int i = diag(ih,0);
    int j = diag(ih,1);
    int ds = dirspin(ih);
    op(j) = ds;
    op(i) = 1-(ds & 1) + 2*(ds>>1);  // operator opposite to ds: from number (2n+1) gives 2n, and from 2n gives (2n+1)
  }
}
void PPstates1(Array<int,2>& sts,
	       int PPs, const Array<int,2>& diag, const Array<int,1>& ops, const Cix& cix) {
  // After we run PPconfigs2 to determine how many pseudo states survive Pauli principle, we compute them all
  //
  //
  int n = diag.extent(0);  // Perturbation order
  sts.resize(PPs,2*n);  // Stores pseudo-particle sequences
  //int Npseudo = Fi.extent(1)-1;
  // Iterate over possible pseudo-particle states
  int ir = 0;  // Which row to write to
  for (int ib=0; ib<cix.Nblock; ++ib){
    int pp = ib;  // Start pseudo-particle index at ib
    for (int i=0; i<2*n; ++i){
      pp = cix.Fbc(ops(i),pp);  // Transition using lookup table
      sts(ir, i) = pp;  // Store the state (adjusted for zero-based indexing)
      //cout<<"ib="<<ib<<",i="<<i<<" ir="<<ir<<" pp="<<pp<<" op["<<i<<"]="<<ops(i)<<" PPs="<<PPs<<endl;
      if (pp == -1) break;
    }
    if (pp != -1) {  // If a valid state sequence is found
      ir += 1;
      if (ir==PPs) break;
    }
  }
}

void PPstates2(Array<int,2>& sts, int PPs, const Array<int,2>& diag, const Array<int,1>& dirspin, const Array<int,2>& Fi) {
  // After we run PPconfigs2 to determine how many pseudo states survive Pauli principle, we compute them all
  //
  //
  int n = diag.extent(0);  // Perturbation order
  sts.resize(PPs,2*n);  // Stores pseudo-particle sequences
  //cout<<"n="<<n<<" PPs="<<PPs<<endl;
  Array<int,1> op(2*n);  // Operators
  // Compute occupation states
  for (int ih = 0; ih<n; ++ih){ // over all hybridizations
    int i = diag(ih,0);
    int j = diag(ih,1);
    int ds = dirspin(ih);
    op(j) = ds;
    op(i) = 1-(ds & 1) + 2*(ds>>1);  // operator opposite to ds
  }
  //cout<<"op="<<op<<endl;
  int Npseudo = Fi.extent(1)-1;
  //cout<<"Npseudo="<<Npseudo<<endl;
  // Iterate over possible pseudo-particle states
  int ir = 0;  // Which row to write to
  for (int ip=0; ip<Npseudo; ++ip) {
    int pp = ip + 1;  // Start pseudo-particle index at ip + 1
    for (int i=0; i<2*n; ++i){
      pp = Fi(op(i),pp);  // Transition using lookup table
      //cout<<"ip="<<ip+1<<" i="<<i<<" pp="<<pp<<endl;
      sts(ir, i) = pp-1;  // Store the state (adjusted for zero-based indexing)
      if (pp == 0) {
	break;
      }
    }
    if (pp != 0) {  // If a valid state sequence is found
      ir += 1;
      if (ir==PPs) break;
    }
  }
}

/*
Array<double,1> SignMatrixElements(const Array<int,2>& diag, const Array<int,1>& dirspin, const Array<int,2>& sts, const Array<double,2>& Fm)
{
  int n = diag.extent(0);  // Perturbation order
  Array<int,1> op(2*n);
  for (int ih = 0; ih < n; ++ih) {
    int i = diag(ih,0);
    int j = diag(ih,1);
    int ds = dirspin(ih);
    op(j) = ds;
    op(i) = 1-(ds & 1) + 2*(ds >> 1);  // operator opposite to ds
  }
  Array<double,1> mem(sts.extent(0));
  for (int ir=0; ir<sts.extent(0); ++ir){
    int pp = sts(ir,2*n-1)+1;
    double mm = Fm(op(0),pp);
    for (int i=0; i<2*n-1; ++i){
      pp = sts(ir,i)+1;
      mm *= Fm(op(i+1),pp);
    }
    mem(ir)=mm;
  }
  return mem;
}
*/

template<int Nt>
int weighted_choice(const TinyVector<double,Nt>& prob_dist, double r){
  double cum = 0.0;
  for (int i = 0; i<Nt; ++i){
    cum += prob_dist[i];
    if (r < cum) {
      return i;
    }
  }
  return Nt-1;  // Fallback in case of rounding errors
}
/*
void Add_Hyb(int& ind1, int& ind2, Array<int,2>& diag, Array<double,1>& times, double beta,
	     const TinyVector<double,2>& random)
{
  int n = diag.extent(0);
  int size_2n = 2*n;
  // randomly select two new times between 0-beta, and ensure that nt1<nt2
  double nt1 = random[0]*beta;
  double nt2 = random[1]*beta;
  if (nt1>nt2) std::swap(nt1,nt2);
  // Find insertion indices using std::lower_bound (equivalent to np.searchsorted)
  auto ind1_it = std::lower_bound(times.begin(), times.end(), nt1);
  auto ind2_it = std::lower_bound(ind1_it, times.end(), nt2);
  ind1 = std::distance(times.begin(), ind1_it);
  ind2 = ind1 + std::distance(ind1_it, ind2_it);

  //cout<<"nt1="<<nt1<<" nt2="<<nt2<<" ind1="<<ind1<<" ind2="<<ind2<<endl;
  
  // Resize diag for one more hybridization adding at the end.
  diag.resizeAndPreserve(TinyVector<int,2>(n+1,2));
  diag(n,0) = ind1;
  diag(n,1) = ind2+1;
  // Adjust indices in diag
  for (int i=0; i<n; ++i){
    if (diag(i,0)>=ind1) diag(i,0)+=1;
    if (diag(i,1)>=ind1) diag(i,1)+=1;
    if (diag(i,0)>ind2) diag(i,0)+=1;
    if (diag(i,1)>ind2) diag(i,1)+=1;
  }
  // Create new times with extra time points inserted
  times.resizeAndPreserve(size_2n+3);
  for (int i=size_2n+2; i>=ind2+2; i--) times(i) = times(i-2);
  times(ind2+1)=nt2;
  for (int i=ind2; i>=ind1+1; i--) times(i) = times(i-1);
  times(ind1) = nt1;
}

int Remove_Hyb(Array<int,2>& diag, Array<double,1>& times, double beta, double random)
{
    int n = diag.extent(0);
    int hyb_to_remove = static_cast<int>(random*n); // Randomly select an element to remove
    int ind1 = diag(hyb_to_remove,0);
    int ind2 = diag(hyb_to_remove,1);
    // Remove the selected elements from diag
    for (int i=hyb_to_remove; i<n-1; ++i){
      diag(i,0) = diag(i+1,0);
      diag(i,1) = diag(i+1,1);
    }
    diag.resizeAndPreserve(TinyVector<int,2>(n-1,2));
    // Adjust indices in diag
    for (int i=0; i<n-1; ++i){
      if (diag(i,0)>ind1) diag(i,0) -= 1;
      if (diag(i,1)>ind1) diag(i,1) -= 1;
    }
    for (int i=0; i<n-1; ++i){
      if (diag(i,0)>=ind2) diag(i,0) -= 1;
      if (diag(i,1)>=ind2) diag(i,1) -= 1;
    }
    // push values of the times array for one entry lower up to ind1, so that ind1 is eliminated
    for (int i=ind1; i<2*n; ++i)   times(i) = times(i+1);
    // and now push all entries for one lower up to ind2, so that ind2 is eliminetaed too
    for (int i=ind2-1; i<2*n-1; ++i) times(i) = times(i+1);
    times.resizeAndPreserve(2*n-1);
    // Adjust times if necessary
    if (times(0) != 0){
        double shift = times(0);
	times(Range::all()) -= shift;
	times(2*n-2) = times(0)+beta;
    }
    return hyb_to_remove;
}
*/

class exp_distribution{
public:
  double lambda;
  exp_distribution(double x0) : lambda(1/x0) {}
  double w(double x){ return exp(-fabs(x)/lambda)/(2*lambda); }
  double random(double W){ return W<0.5 ? lambda*log(2*W) : -lambda*log(2*(1-W));}
};


template<bool ImagAxis>
class ChangeHybridization{
  // Class that adds and removes hybridization on both real and imaginary axis.
  // While the topology of the diagram and dirspin array is treated in the same way in both axis
  // the frequency and time array needs to be treated very differently.
  // The time array needs to be sorted and uniformly distributed on the ring.
  // The frequency array should be uniform (but not sorted) or distributed by a known distribution,
  // that can be used for trial step probablity.
public:
  double beta, cutoffw;
  deque<double> added_freq, removed_freq; // for real axis we need to remember all added/removed frequencies
  exp_distribution distribution; // for real axis we propose frequency with this distribution
  
  ChangeHybridization(double _beta_, double _cutoffw_, double distr_width) :
    beta(_beta_), cutoffw(_cutoffw_), distribution(distr_width) {};
    
  void reset_add_rm(){
    if (not ImagAxis) {
      added_freq.clear();
      removed_freq.clear();
    }
  }

  double distribution_x_jump()
  {
    if constexpr (ImagAxis){
      return 1.0;
    }else{
      double P=1.0;
      for (int i=0; i<added_freq.size(); ++i) P *= 2*distribution.w(added_freq[i]);
      for (int i=0; i<removed_freq.size(); ++i) P /= 2*distribution.w(removed_freq[i]);
      return P;
    }
  }
  std::tuple<int,int,double,double,double> 
  try_insert_hyb(const Array<int,2>& diag, const Array<double,1>& times, const TinyVector<double,4>& random){
    int n = diag.extent(0);
    int size_2n = 2*n;
    // Randomly pick two insertion indices (1 to size_2n)
    int ind1 = 1+static_cast<int>(random[0]*size_2n);
    int ind2 = 1+static_cast<int>(random[1]*size_2n);
    // Ensure ind1 <= ind2
    if (ind1 > ind2) std::swap(ind1, ind2);
    if (ind2>size_2n) ind2=size_2n;
    if (ind1>size_2n) ind1=size_2n;
    
    if constexpr (ImagAxis){
      // Retrieve time intervals
      double t0 = times(ind1-1);
      double t1 = times(ind1);
      // Generate new time points within the selected intervals
      double t1_new = t0 + (t1 - t0)*random[2];
      double average_spacing = beta/size_2n;
      double P1 = (t1-t0)/average_spacing;  // Ratio of actual length to average length
      double t2_new, P2, t2, t3;
      if (ind2 != ind1) {
	t2 = times(ind2-1);
	t3 = times(ind2);
	t2_new = t2 + (t3 - t2)*random[3];
	P2 = (t3-t2)/average_spacing;
      } else {
	t2_new = t0 + (t1 - t0)*random[3];
	if (t1_new > t2_new) std::swap(t1_new, t2_new);
	P2 = P1;
      }
      //cout<<"....ind1="<<ind1<<",ind2="<<ind2<<",t0="<<t0<<",t1="<<t1<<",t1_new="<<t1_new<<",P1="<<P1<<",t2="<<t2<<",t3="<<t3<<",t2_new="<<t2_new<<",P2="<<P2<<endl;
      return {ind1, ind2, t1_new, t2_new, P1*P2};
    }else{
      double new_freq = cutoffw*(2*random[2]-1);
      //double P = 0.125*n*n;  // we use only the 1/8*n^2 from the total (2*n)^2 * (N_bath*2*L)/(n+1)
      double P = 1.0;
      if (ind1==ind2) P *= 2.0;
      return {ind1, ind2, new_freq, 0, P};
    }
  }
  std::tuple<int,double> 
  try_remove_hyb(const Array<int,2>& diag, const Array<double,1>& times, double random)
  {
    int n = diag.extent(0);
    int size_2n = n*2;
    // Randomly pick a hybridization to remove
    int hyb_to_remove = random*n;
    // Retrieve corresponding indices
    int i1 = diag(hyb_to_remove,0);
    int i2 = diag(hyb_to_remove,1);
    // Ensure i1 <= i2
    if (i1 > i2) std::swap(i1, i2);

    if constexpr (ImagAxis){
      //cout<<"n="<<n<<" size_2n="<<size_2n<<" htr="<<hyb_to_remove<<" r="<<random<<" i1="<<i1<<" i2="<<i2<<endl;
      double average_size = beta/(size_2n - 2);
      double P1, P2;
      if (i2-i1!=1 && i2-i1!=size_2n-1) {
	double dl1 = times((i1+1)%size_2n) - times((i1-1+size_2n)%size_2n);
	if (i1==0 || i1==size_2n-1) dl1 += beta;
	double dl2 = times((i2+1)%size_2n) - times((i2-1+size_2n)%size_2n);
	if (i2==0 || i2==size_2n-1) dl2 += beta;
	//cout<<"dl1="<<dl1<<" dl2="<<dl2<<" as="<<average_size<<endl;
	P1 = dl1/average_size;
	P2 = dl2/average_size;
	if (P1<1e-50 or P2<1e-50){
	  cout<<"WARNING P1="<<P1<<" P2="<<P2<<" dl1="<<dl1<<" dl2="<<dl2<<" i1="<<i1<<" i2="<<i2<<" average_size="<<average_size<<endl;
	}
      } else {
	double t0, t1;
	if (i2 - i1 == 1) {
	  t1 = times(i2+1);
	  t0 = (i1 > 0) ? times(i1-1) : times(size_2n-1) - beta;
	} else { // i1 = 0, i2 = 2*n - 1
	  t1 = times(1);
	  t0 = times(size_2n-2) - beta;
	}
	P1 = (t1-t0)/average_size;
	P2 = P1;
	if (P1<1e-50){
	  cout<<"WARNING P1<1e-50 t0="<<t0<<" t1="<<t1<<" i1="<<i1<<" i2="<<i2<<" average_size="<<average_size<<endl;
	}
      }
      if (P1<1e-50) P1=1e-50;
      if (P2<1e-50) P2=1e-50;
      return {hyb_to_remove, 1.0/(P1*P2)};
    }else{
      //double P= 1.0/(0.125*(n-1.)*(n-1.)); // we use only the (2n)^2 from the total (2*n)^2 * (N_bath*2*L)/(n+1)
      double P = 1.0;
      if (i2-i1==1 or i2-i1==size_2n-1) P *= 0.5;
      return {hyb_to_remove, P};
    }
  }
  void execute_insert_hyb(Array<int,2>& diag, Array<double,1>& times, int ind1, int ind2, double nt1, double nt2)
  {
    int n = diag.extent(0);
    int size_2n = 2*n;
    // Resize diag for one more hybridization adding at the end.
    diag.resizeAndPreserve(TinyVector<int,2>(n+1,2));
    diag(n,0) = ind1;
    diag(n,1) = ind2+1;
    // Adjust indices in diag
    for (int i=0; i<n; ++i){
      if (diag(i,0)>=ind1) diag(i,0)+=1;
      if (diag(i,1)>=ind1) diag(i,1)+=1;
      if (diag(i,0)>ind2) diag(i,0)+=1;
      if (diag(i,1)>ind2) diag(i,1)+=1;
    }
    if constexpr (ImagAxis){
      // Create new times with extra time points inserted
      times.resizeAndPreserve(size_2n+3);
      for (int i=size_2n+2; i>=ind2+2; i--) times(i) = times(i-2);
      times(ind2+1)=nt2;
      for (int i=ind2; i>=ind1+1; i--) times(i) = times(i-1);
      times(ind1) = nt1;
    }else{
      times.resizeAndPreserve(n+2);
      times(n+1) = nt1;
      added_freq.push_back(nt1);
    }
  }
  void execute_remove_hyb(Array<int,2>& diag, Array<double,1>& times, int hyb_to_remove)
  {
    int n = diag.extent(0);
    int size_2n = 2*n;
    // Get indices to remove
    int ind1 = diag(hyb_to_remove,0);
    int ind2 = diag(hyb_to_remove,1);
    // Remove the selected elements from diag
    for (int i=hyb_to_remove; i<n-1; ++i){
      diag(i,0) = diag(i+1,0);
      diag(i,1) = diag(i+1,1);
    }
    diag.resizeAndPreserve(TinyVector<int,2>(n-1,2));
    // Adjust indices in diag
    for (int i=0; i<n-1; ++i){
      if (diag(i,0)>ind1) diag(i,0) -= 1;
      if (diag(i,1)>ind1) diag(i,1) -= 1;
    }
    for (int i=0; i<n-1; ++i){
      if (diag(i,0)>=ind2) diag(i,0) -= 1;
      if (diag(i,1)>=ind2) diag(i,1) -= 1;
    }
    if constexpr (ImagAxis){
      // push values of the times array for one entry lower up to ind1, so that ind1 is eliminated
      for (int i=ind1; i<2*n; ++i)   times(i) = times(i+1);
      // and now push all entries for one lower up to ind2, so that ind2 is eliminetaed too
      for (int i=ind2-1; i<2*n-1; ++i) times(i) = times(i+1);
      times.resizeAndPreserve(2*n-1);
      // Adjust times if necessary
      if (times(0) != 0){
	double shift = times(0);
	times(Range::all()) -= shift;
	times(2*n-2) = times(0)+beta;
      }
    }else{
      removed_freq.push_back(times(hyb_to_remove+1));
      for (int i=hyb_to_remove+1; i<n; ++i) times(i) = times(i+1);
      times.resizeAndPreserve(n);
    }
  }
};

template<class RNG, bool ImagAxis>
int propose_new_config4(Array<int,2>& diag, Array<int,1>& dirspin, Array<double,1>& times, int& ind1, int& ind2, ChangeHybridization<ImagAxis>& ch,
			const TinyVector<double,3>& prob_move, double some_factor, double prob_exchange,
			int Nbaths, int minorder, int maxorder, /*double beta,*/ RNG& drand)
{
  int n = diag.extent(0);
  double xi0=drand(), xi1=drand();
  double c_old=1.0;
  int jump = 0;
  //////// ASK Yueyi ????? ADDED Jun 14 2026
  if (maxorder==minorder){
    jump=0;
  }else if (n == maxorder) {
    c_old = prob_move[0]+prob_move[1];
    jump = weighted_choice(TinyVector<double,2>(prob_move[0]/c_old,prob_move[1]/c_old), xi0)-1;
  } else if (n == minorder) {
    c_old = prob_move[1]+prob_move[2];
    jump = weighted_choice(TinyVector<double,2>(prob_move[1]/c_old,prob_move[2]/c_old), xi0);
  } else {
    jump = weighted_choice(prob_move, xi0)-1;
  }
  int n_new = n + jump;
  double c_new=1.0;
  if (n_new==maxorder)
    c_new = prob_move[0]+prob_move[1];
  else if (n_new==minorder)
    c_new = prob_move[1]+prob_move[2];
  //if (n_new==1)
  //  c_new=prob_move[2];
  
  double P_to_change = c_old/c_new;
  int hyb_to_remove=-1;
  double t1_new=0, t2_new=0;
  if (jump == 1) {
    P_to_change *= some_factor;
    //if (ioption!=1){
    double Pchange;
    std::tie(ind1,ind2,t1_new,t2_new,Pchange) = ch.try_insert_hyb(diag,times,TinyVector<double,4>(drand(),drand(),drand(),drand()));
    P_to_change *= Pchange;
    //}
  } else if (jump == -1) {
    P_to_change *= 1/some_factor;
    //if (ioption!=1){
    double Pchange;
    std::tie(hyb_to_remove,Pchange) = ch.try_remove_hyb(diag,times,drand());
    P_to_change *= Pchange;
    //}
  } else { // jump==0
    P_to_change = 1.0;
  }
  if (P_to_change < xi1){
    if (DEBUG) cout<<"... P_to_change="<<P_to_change<<" jump="<<jump<<" xi1="<<xi1<<" xi0="<<xi0<<endl;//" icounter="<<drand.icounter<<endl;
    return 0;
  }
  // proposing a new config
  if (jump == 1) {
    //if (ioption==1)
    //  Add_Hyb(ind1,ind2,diag,times,beta,TinyVector<double,2>(drand(),drand()));
    //else
    ch.execute_insert_hyb(diag, times, ind1, ind2, t1_new, t2_new);
    // modify dirspin, and randomly assign the inserted hybridization
    dirspin.resizeAndPreserve(n+1);
    dirspin(n) = static_cast<int>(drand()*Nbaths*2);
  } else if (jump == -1) {
    //if (ioption==1)
    //hyb_to_remove = Remove_Hyb(diag,times,beta,drand());
    //else
    ch.execute_remove_hyb(diag, times, hyb_to_remove);
    // correct dirspin
    for (int i=hyb_to_remove; i<n-1; ++i) dirspin(i) = dirspin(i+1);
    dirspin.resizeAndPreserve(n-1);
    ind1=hyb_to_remove; ind2=0;
  } else {
    if (n>2 and drand()<prob_exchange){
      // swap two hybridizations, but keep times. Change dirspin to a random value for the two hybridizations exchanged.
      int ver1 = drand()*n;
      int ver2 = drand()*(n-1);
      if (ver2>=ver1) ver2+=1;
      int ifdagger1 = drand()*2;
      int ifdagger2 = drand()*2;
      std::swap(diag(ver1,ifdagger1), diag(ver2,ifdagger2));
      if (diag(ver1,0)>diag(ver1,1)) std::swap(diag(ver1,0),diag(ver1,1));
      if (diag(ver2,0)>diag(ver2,1)) std::swap(diag(ver2,0),diag(ver2,1));
      // choose completely random dirspins
      //for (int i=0; i<n; ++i) dirspin(i) = static_cast<int>(drand()*Nbaths*2);
      dirspin(ver1) = static_cast<int>(drand()*Nbaths*2);
      dirspin(ver2) = static_cast<int>(drand()*Nbaths*2);
      ind1 = 2*ver1+ifdagger1; ind2 = 2*ver2+ifdagger2;
    }else{
      for (int i=0; i<n; ++i) dirspin(i) = static_cast<int>(drand()*Nbaths*2);
      ind1=ind2=-1;
    }
  }
  if (DEBUG) cout<<"...jump="<<jump<<" diag="<<diag<<" times="<<times<<" dirspin="<<dirspin<<" xi0="<<xi0<<" xi1="<<xi1<<endl;//" icounter="<<drand.icounter<<endl;
  return jump+2;
}

//Array<fsec,1> t_propose_new(3);
//fsec t_Pauli;
//fsec t_skelet;
template<class RNG, bool ImagAxis>
int Next_valid_config_bold4(int itt, Array<int,2>& diag, Array<int,1>& dirspin, Array<double,1>& times, int& PPs, ChangeHybridization<ImagAxis>& ch,
			    const Cix& cix, //int Nbaths, const Array<int,2>& Fi,
			    const TinyVector<double,3>& prob_move, double some_factor, double prob_exchange, int minorder, int maxorder, RNG& drand)
{
  bool is_skeleton = true; // Original diag is valid, hence PPs should be nonzero when entering
  ch.reset_add_rm(); // all previously attempted frequencies that were added or removed are deleted,
                   // so that we only save those from this iteration.
  Skeletons skl;
  int itd=0;
  while (true) {
    int ind1,ind2;
    //auto t1 = high_resolution_clock::now();  // Start time tracking
    // Propose a new configuration
    int new_found = propose_new_config4(diag,dirspin,times,ind1,ind2,ch,
					prob_move,some_factor,prob_exchange,cix.Nbath,minorder,maxorder,drand);
    //auto t2 = high_resolution_clock::now();  // End time tracking
    int jump = new_found - 2;
    //t_propose_new(jump+1) += (t2 - t1);  // Time in seconds
    if (new_found) {  // A new valid diagram is found
      is_skeleton = false;
      if (jump==1 and ind1==ind2) {
	continue;  // Skip skeleton check for NCA insertion
      }
      //auto t3 = high_resolution_clock::now();
      PPs = PPconfigs2(diag, dirspin, cix);
      //is_Pauli_allowed = PPs;
      //auto t4 = high_resolution_clock::now();
      //t_Pauli += (t4 - t3);
      if (PPs) {
	//auto t5 = high_resolution_clock::now();
	// we could check ind1==ind2 when adding, or abs(diag(i,0)-diag(i,1))==1 when removing.
	is_skeleton = (diag.extent(0)==1) or (not skl.T_Is_Reducible(diag));
	//auto t6 = high_resolution_clock::now();
	//t_skelet += (t6 - t5);
	if (is_skeleton) {
	  break;
	}
      }
    } else {  // Old configuration remains valid
      if (is_skeleton and PPs) {
	//if (DEBUG1) cout<<"itt="<<itt<<" itd="<<itd<<" going back to old dig and apparently is_skeleton="<<is_skeleton<<" Pauli_allowed="<<PPs<<endl;
	return 0;
      }
    }
    //if (DEBUG1) cout<<"itt="<<itt<<"itd="<<itd<<" PPs="<<PPs<<" is_skeleton="<<is_skeleton<<endl;
    itd+=1;
  }
  if (DEBUG) cout<<"itt="<<itt<<" trying newdiag="<<diag<<" newdirspin="<<dirspin<<" ntimes="<<times<<" PPs="<<PPs<<endl;
  // Check if new configuration is same as old
  //bool new_equal_old = all(cdiag == diag) && all(cdirspin == dirspin) && all(abs(ctimes - times) < 1e-12);
  
  return 1;
}
#endif //_TOPOLOGY
