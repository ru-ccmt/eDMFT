#include <blitz/array.h>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>  // For std::lower_bound (binary search)
#define _GSL
#include "Random2.h"
#include <deque>
#include "basics.h"
using namespace blitz;
using namespace std;
const bool DEBUG=false;
//#include "frequency_loops_real_axis.h"
#include "skeletons.h"
#include "topologymove.h"

Array<int,1> Skeletons::Freq(16);
double cutoffw=10;

template<class RNG, bool ImagAxis>
int Simple_propose_new_diag(Array<int,2>& diag, Array<int,1>& dirspin, Array<double,1>& freq, int& ind1, int& ind2, ChangeHybridization<ImagAxis>& ch,
			    const TinyVector<double,3>& prob_move, double prob_exchange,
			    int Nbaths, int minorder, int maxorder, RNG& drand)
{
  int n = diag.extent(0);
  double xi0=drand(), xi1=drand();
  double c_old=1.0;
  int jump = 0;
  if (n == maxorder) {
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
  if (n_new==1)
    c_new=prob_move[2];

  double P_to_change = c_old/c_new;
  int hyb_to_remove=-1;
  double t1_new=0, t2_new=0;
  if (jump == 1) {
    double Pchange;
    std::tie(ind1,ind2,t1_new,t2_new,Pchange) = ch.try_insert_hyb(diag,freq,TinyVector<double,4>(drand(),drand(),drand(),drand()));
    P_to_change *= Pchange;
    //P_to_change *= 8*n*n/(n+1);
    P_to_change *= 64./(n+1.); // because 0.125*n^2 is inside Pchange
  } else if (jump == -1) {
    double Pchange;
    std::tie(hyb_to_remove,Pchange) = ch.try_remove_hyb(diag,freq,drand());
    P_to_change *= Pchange;
    //P_to_change *= n/(8.*(n-1.)*(n-1.));
    P_to_change *= n/64.;  // becatse 1/(0.125*(n-1)^2)
  } else { // jump==0
    P_to_change = 1.0;
  }
  if (P_to_change < xi1){
    if (DEBUG) cout<<"... P_to_change="<<P_to_change<<" jump="<<jump<<" xi1="<<xi1<<" xi0="<<xi0<<endl;//" icounter="<<drand.icounter<<endl;
    return 0;
  }
  // proposing a new config
  if (jump == 1) {
    ch.execute_insert_hyb(diag, freq, ind1, ind2, t1_new, t2_new);
    // modify dirspin, and randomly assign the inserted hybridization
    dirspin.resizeAndPreserve(n+1);
    dirspin(n) = static_cast<int>(drand()*Nbaths*2);
  } else if (jump == -1) {
    ch.execute_remove_hyb(diag, freq, hyb_to_remove);
    // correct dirspin
    for (int i=hyb_to_remove; i<n-1; ++i) dirspin(i) = dirspin(i+1);
    dirspin.resizeAndPreserve(n-1);
    ind1=hyb_to_remove; ind2=0;
  } else {
    if (n>1 and drand()<prob_exchange){
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
  if (DEBUG) cout<<"...jump="<<jump<<" diag="<<diag<<" freq="<<freq<<" dirspin="<<dirspin<<" xi0="<<xi0<<" xi1="<<xi1<<endl;//" icounter="<<drand.icounter<<endl;
  return jump+2;
}

template<class RNG, bool ImagAxis>
int Simple_Next_valid_diag(int itt, Array<int,2>& diag, Array<int,1>& dirspin, Array<double,1>& times, int& PPs, ChangeHybridization<ImagAxis>& ch,
			   int Nbaths, const Array<int,2>& Fi, const TinyVector<double,3>& prob_move, double prob_exchange,
			   int minorder, int maxorder, RNG& drand)
{
  bool is_skeleton = true; // Original diag is valid, hence PPs should be nonzero when entering
  ch.reset_add_rm(); // all previously attempted frequencies that were added or removed are deleted,
                   // so that we only save those from this iteration.
  Skeletons skl;
  int itd=0;
  while (true) {
    int ind1,ind2;
    //auto t1 = high_resolution_clock::now();  // Start time tracking
    int new_found = Simple_propose_new_diag(diag,dirspin,times,ind1,ind2,ch,
					    prob_move,prob_exchange,Nbaths,minorder,maxorder,drand);
    //auto t2 = high_resolution_clock::now();  // End time tracking
    int jump = new_found - 2;
    //t_propose_new(jump+1) += (t2 - t1);  // Time in seconds
    if (new_found) {  // A new valid diagram is found
      is_skeleton = false;
      //if (jump==1 and ind1==ind2) {
      //continue;  // Skip skeleton check for NCA insertion
      //}
      //auto t3 = high_resolution_clock::now();
      PPs = PPconfigs2(diag, dirspin,Fi);
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

void diag_resort(Array<int,2>& tdiag, Array<int,1>& tdirspin,
		 const Array<int,2>& diag, const Array<int,1>& dirspin)
{// Given a diagram  diag=[(4,1),(6,3),(0,2),(7,5)] it resorts it such that
  //  we always have: diag(i,1)>diag(i,0) and that diag(i+1,0)>diag(i,0).
  // In above example, we would get: tdiag=[(0,2),(1,4),(3,6),(5,7)]
  // It also resorts dirspin so that it is compatible with resorting, i.e.,
  // tdiag+tdirspin is equivalent to diag+dirspin.
  int n = diag.extent(0);
  std::vector<int> indices(n);
  for (int i=0; i<n; ++i) indices[i] = i;
  std::sort(indices.begin(), indices.end(), [&](int i, int j){return diag(i,0) < diag(j,0);});
  tdiag.resize(diag.shape());
  tdirspin.resize(dirspin.size());
  for (int i=0; i<n; ++i){
    tdiag(i,0) = diag(indices[i],0);
    tdiag(i,1) = diag(indices[i],1);
    tdirspin(i) = dirspin(indices[i]);
  }
}

void quickrotate(Array<int,2>& rdiag, Array<int,1>& rdirspin,
		 const Array<int,2>& diag, const Array<int,1>& dirspin)
{// It rotates diagram around the ring in clockwise direction for one kink, and than resorts it
  // to get equivalent diagra, i.e., rdiag+rdirspin is equivalent to diag+dirspin
  rdiag.resize(diag.shape());
  rdirspin.resize(dirspin.size());
  int n = diag.extent(0);
  int size_2n = 2*n;
  for (int i=0; i<n; ++i){
    rdiag(i,0) = (diag(i,0)-1+size_2n)%size_2n; // rotate for one unit 
    rdiag(i,1) = (diag(i,1)-1+size_2n)%size_2n;
    rdirspin(i) = dirspin(i);     // dirspin is kept constant at this step
    if (rdiag(i,1)<rdiag(i,0)){   // but can change direction if second index is smaller than first index of hybridization
      swap(rdiag(i,1),rdiag(i,0));// and therefore hybridization is flipped and dirspin with it.
      rdirspin(i) = rdirspin(i)^1;// flipping direction in dirspin is equivalent to ^1 operation
    }
  }
  // Here we resort array so that vertex with smaller first index comes first
  Array<int,2> tdiag = rdiag.copy();
  Array<int,1> tdirspin = rdirspin.copy();
  diag_resort(rdiag,rdirspin,tdiag,tdirspin);
}
int lexicographic_smaller(const Array<int,2>& diag1, const Array<int,2>& diag2){
  // returns -1,0,1 if diag1<diag2, diag1==diag2 or diag1>diag2 in lexicographic order.
  if (diag1.extent(0) != diag2.extent(0)) cout<<"ERROR: Comparing diags not of equal size"<<endl;
  for (int i=0; i<diag1.extent(0); ++i){
    for (int j=0; j<diag1.extent(1); ++j){
      if (diag1(i,j) < diag2(i,j)) return -1;
      else if (diag1(i,j) > diag2(i,j)) return 1;
    }
  }
  return 0;
}
int lexicographic_smaller(const Array<int,1>& arr1, const Array<int,1>& arr2) {
  // returns -1,0,1 if arr1<arr2, arr1==arr2 or arr1>arr2 in lexicographic order.
  if (arr1.extent(0) != arr2.extent(0)) cout<<"ERROR: Comparing arrays not of equal size"<<endl;
  int n = arr1.extent(0);
  for (int i = 0; i < n; ++i) {
    if (arr1(i) < arr2(i)) return -1;
    else if (arr1(i) > arr2(i)) return 1;
  }
  return 0;
}
double best_rotation_config(Array<int,2>& bestdiag, Array<int,1>& bestdirspin,
			  const Array<int,2>& diag, const Array<int,1>& dirspin)
{// It rotates diagram diag+dirspin for 2n-times to generate all equivalent
  // representations and then sorts each diagram to its ordered way (diag(i,0)<diag(i,1) and diag(i,0)<diag(i+1,0)).
  // Finally it compares all such equivalent representations of a diagram+dirspin,
  // and returns the one which is lexicogrphically the smallest.
  int n = diag.extent(0);
  bestdiag.resize(n,2); 
  bestdirspin.resize(n);
  
  diag_resort(bestdiag,bestdirspin, diag,dirspin);
  
  Array<int,2> tdiag(diag.shape());
  Array<int,1> tdirspin(dirspin.size());
  Array<int,2> rdiag = diag.copy();
  Array<int,1> rdirspin = dirspin.copy();
  int deg_best=1;
  //cout<<"start with "<<diag<<","<<dirspin<<" best="<<bestdiag<<","<<bestdirspin<<endl;
  for (int i=0; i<2*n-1; ++i){
    quickrotate(tdiag, tdirspin, rdiag, rdirspin);
    int is_current_smaller = lexicographic_smaller(tdiag,bestdiag);
    if (is_current_smaller==-1){
      bestdiag = tdiag;
      bestdirspin = tdirspin;
      deg_best = 1;
    }else if(is_current_smaller==0){ // diagrams equal
      int is_dir_smaller = lexicographic_smaller(tdirspin,bestdirspin);
      //cout<<"  dir_smaller="<<is_dir_smaller<<endl;
      if (is_dir_smaller==-1){
	bestdiag = tdiag;
	bestdirspin = tdirspin;
	deg_best = 1;
      }else if (is_dir_smaller==0){ // both diag and dirspin equal, hence higher symmetry
	deg_best += 1;
      }
    }
    rdiag = tdiag;
    rdirspin = tdirspin;
    //cout<<tdiag<<" "<<tdirspin<<" is="<<is_current_smaller<<" best="<<bestdiag<<","<<bestdirspin<<" deg_best="<<deg_best<<endl;
  }
  return deg_best;
}
unsigned long long diag_spin_to_key(const Array<int,2>& diag, const Array<int,1>& dirspin, int maxorder)
{// given diag+dirspin it produces a has, i.e., long integer that encodes the diagram and dirspin.
  int base = 2 * maxorder;
  int base2 = base * base;
  int full_base = base2 * 4;
  unsigned long long key = 0;
  for (int i = 0; i < diag.extent(0); ++i) {
    int unit = ((diag(i,0) * base + diag(i,1)) * 4 + dirspin(i));
    key = key * full_base + unit;
  }
  return key;
}
void key_to_diag_spin(Array<int,2>& diag, Array<int,1>& dirspin,
		      unsigned long long key, int maxorder)
{// given a unique key obtained by diag_spin_to_key() returns diag+dirspin from the key, i.e.,
  // it is the inverse of diag_spin_to_key()
  int base = 2 * maxorder;
  int base2 = base * base;
  int full_base = base2 * 4;
  // First determine the number of hybridizations
  int n = 0;
  unsigned long long temp = key;
  while (temp > 0) {
    temp /= full_base;
    n++;
  }
  diag.resize(n, 2);
  dirspin.resize(n);
  for (int i = n-1; i >= 0; --i) {
    int unit = key % full_base;
    key /= full_base;
    dirspin(i) = unit % 4;
    int combined = unit / 4;
    diag(i,1) = combined % base;
    diag(i,0) = combined / base;
  }
}
int main() {
  int minorder=1, maxorder=5;
  TinyVector<double,3> prob_move(0.5,0,0.5);
  
  const int max_order=20;
  Array<double,1> factorial(max_order+1);
  factorial(0)=1;
  for (int i=1; i<factorial.size(); ++i) factorial(i) = factorial(i-1)*i;
  
  int Nbath=2, Npseudo=4;
  Array<int,2> Fi(2*Nbath,Npseudo+1);
  //  1=|0> 2=|up> 3=|dn> 4=|2>
  //  0=c^+_up 1=c_up 2=c^+_dn 3=c_dn
  //  c^+_up|0>=|up> =>  F[0,1]=2
  //  c^+_up|dn>=|2> =>  F[0,3]=4
  //  c_up|up>=|0>   =>  F[1,2]=1
  //  c_up|2> = |dn> =>  F[1,4]=3
  //  c^+_dn|0>=|dn> =>  F[2,1]=3
  //  c^+_dn|up>=-|2> => F[2,2]=4
  //  c_dn|dn>=|0>    => F[3,3]=1
  //  c_dn|2> =-|up>  => F[3,4]=2
  Fi = 0, 2, 0, 4, 0,  // c^+_up
       0, 0, 1, 0, 3,  // c_up
       0, 3, 4, 0, 0,  // c^+_dn
       0, 0, 0, 1, 2;  // c_dn
  {
    int two_2_i = 1;
    for (int i=0; i<Skeletons::Freq.extent(0); ++i){
      Skeletons::Freq(i) = two_2_i;
      two_2_i *= 2;
    }
  }

  int real_iseed=0;//time(0);
  RanGSL drand(real_iseed);
  
  int n=1;
  Array<int,2> diag(n,2);
  Array<int,1> dirspin(n);
  diag = 0,1;
  dirspin = 0;
  Array<double,1> freq(n+1);
  freq=0;

  double prob_exchange=0.5;
  double distr_width=1.0;
  ChangeHybridization<false> ch(25., cutoffw, distr_width);
  int ind1,ind2;
  //int Nbin=10;
  //Array<double,2> hist(maxorder-minorder+1,Nbin); hist=0;
  
  Array<double,1> Vphys(maxorder-minorder+1); Vphys=0;
  Array<int,1> ops;
  Array<int,2> ppstates;
  
  int PPs = PPconfigs2(diag,dirspin,Fi);
  dirspin2ops(ops,diag, dirspin);
  PPstates1(ppstates,PPs,diag,ops,Fi);

  Array<int,2> bestdiag;
  Array<int,1> bestdirspin;
  std::unordered_map<unsigned long long, double> diagMap;
  double sym = best_rotation_config(bestdiag,bestdirspin,diag,dirspin);
  unsigned long long key = diag_spin_to_key(diag,dirspin,maxorder);
  
  int Nitt=1000000;
  for (int itt=0; itt<Nitt; ++itt){
    int n = diag.extent(0);
    
    //int new_found=Simpler_propose_new_diag(diag,dirspin,freq,ind1,ind2,
    //    				     prob_move,prob_exchange,Nbath,minorder,maxorder,drand);


    //int new_found = Simple_propose_new_diag(diag,dirspin,freq,ind1,ind2,ch,
    //					      prob_move,1.0,prob_exchange,Nbath,minorder,maxorder,drand);

    //int new_found = propose_new_config4(diag, dirspin, freq, ind1, ind2, ch,
    //					  prob_move,1.0,prob_exchange,Nbath,minorder,maxorder,drand);

    //int new_found = Next_valid_config_bold4(itt, tdiag, tdirspin, tfreq, tPPs, ch,
    //					    Nbath, Fi, prob_move, 1.0, prob_exchange, minorder, maxorder, drand);


    int new_found = Simple_Next_valid_diag(itt,diag,dirspin,freq,PPs,ch,
					   Nbath,Fi,prob_move,prob_exchange,minorder,maxorder,drand);
    
    if (new_found){
      // finds lexicographically minimal representation for each diagram by rotating diagram and resorting it.
      sym = best_rotation_config(bestdiag,bestdirspin,diag,dirspin);
      // The best diagram is then encoded into key to save it in map (Python dictionary).
      key = diag_spin_to_key(bestdiag,bestdirspin,maxorder);
    }
    diagMap[key] += sym; // here hashtable (map) is updated. We add symmetry factor (1,2,3,..) to ensure that
    // diagrams with 1/3 symmetry factors get sym=3 and hence appear equally probably in the hash table
    int nn=diag.extent(0);
    Vphys(nn-minorder) += PPs*sym;
  }
  Vphys *= 4.0/Vphys(0);
  cout<<"Vphys="<<Vphys<<endl;
  
  vector<unsigned long long> keys;
  for (const auto& pair : diagMap) keys.push_back(pair.first);
  
  std::sort(keys.begin(), keys.end());
  
  Array<int,2> xdiag;
  Array<int,1> xdirspin;
  double mean;
  for (const auto& [key, value] : diagMap) mean += value;
  mean /= diagMap.size();
  int no=1, indx=0;
  for (int i=0; i<keys.size(); ++i){
    unsigned long long key = keys[i];
    double value = diagMap[key];
    key_to_diag_spin(xdiag,xdirspin,key,maxorder);
    if (xdiag.extent(0)!=no){
      no = xdiag.extent(0);
      indx = 1;
    }else indx+=1;
    cout << left<<setw(3)<<indx <<" "<<left<<setw(8)<<value/mean<<" : "<<xdiag << ", "<< xdirspin <<std::endl;
  }
  return 0;
}
