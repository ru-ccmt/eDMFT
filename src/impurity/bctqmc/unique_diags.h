/*
#include <map>
#include <numeric>
#include "basics.h"
#include <blitz/array.h>
#include "wmpi.h"
#define FMT_HEADER_ONLY
#include "fmt/core.h"
#include "fmt/format.h"
#include "wdiags.h"
MPIcls *mpi;
using namespace std;
using namespace blitz;
*/

deque<string> ReadDiags(const string& filename)
{
  ifstream inp(filename);
  deque<string> accumulated_new_diags;
  string str_diag;
  while(inp){
    inp>>str_diag;
    accumulated_new_diags.push_back(str_diag);
    if (not inp or inp.eof()) break;
  }
  return accumulated_new_diags;
}
void Convert_Diags_2_Array(const deque<string>& accumulated_new_diags, int N_min_order, int N_max_order, deque<Array<unsigned char,2>>& diagsv)
{
  vector<deque<vector<unsigned char>>> diags_by_size(N_max_order-N_min_order+1);
  for (int i=0; i<accumulated_new_diags.size(); ++i){
    vector<unsigned char> diag_pseudo = String2diagarray(accumulated_new_diags[i]);
    int Norder = diag_pseudo.size()/4;
    if (Norder>=N_min_order and Norder<=N_max_order)
      diags_by_size[Norder-N_min_order].push_back(diag_pseudo);
  }
  int jit=0;
  diagsv.resize(N_max_order-N_min_order+1);
  for (int Norder=N_min_order; Norder<=N_max_order; ++Norder){
    int iorder=Norder-N_min_order;
    deque<vector<unsigned char>>& dd = diags_by_size[iorder];
    diagsv[iorder].resize(dd.size(),4*Norder);
    for (int i=0; i<dd.size(); ++i)
      for (int j=0; j<dd[i].size(); ++j) diagsv[iorder](i,j)=dd[i][j];
  }
}
void FindWhichPseudoAndJumps(const Array<unsigned char,2>& diags1, vector<string>& p_repeat)
{
  int Norder = diags1.extent(1)/4;
  p_repeat.resize(diags1.extent(0));
  Array<unsigned char,1> ps(2*Norder);
  Array<signed char,1> dp(2*Norder);
  Array<signed char,1> ddp(2*Norder);
  for (int i=0; i<diags1.extent(0); ++i){
    ps = diags1(i,Range(2*Norder,4*Norder-1));
    for (int j=0; j<2*Norder-1; ++j) dp(j) = ps(j+1)-ps(j);
    dp(2*Norder-1) = ps(0)-ps(2*Norder-1);
    //
    for (int j=0; j<2*Norder-2; ++j) ddp(j) = ps(j+2)-ps(j);
    ddp(2*Norder-2) = ps(0)-ps(2*Norder-2);
    ddp(2*Norder-1) = ps(1)-ps(2*Norder-1);
    //
    sort(ps.data(), ps.data()+ps.size());
    sort(dp.data(), dp.data()+dp.size());
    sort(ddp.data(), ddp.data()+ddp.size());
    //
    //p_repeat[i] = fmt::format("{0};{1}", fmt::join(ps.data(), ps.data()+ps.size(), ","), fmt::join(dp.data(), dp.data()+dp.size(), ","));
    p_repeat[i] = fmt::format("{0};{1};{2}", fmt::join(ps.data(), ps.data()+ps.size(), ","), fmt::join(dp.data(), dp.data()+dp.size(), ","), fmt::join(ddp.data(), ddp.data()+ddp.size(), ","));
  }
}
void FindWhichPseudoAndJumps(const deque<Array<unsigned char,1>>& diags1, vector<string>& p_repeat)
{
  int Norder = diags1[0].extent(0)/4;
  p_repeat.resize(diags1.size());
  Array<unsigned char,1> ps(2*Norder);
  Array<signed char,1> dp(2*Norder);
  Array<signed char,1> ddp(2*Norder);
  for (int i=0; i<diags1.size(); ++i){
    ps = diags1[i](Range(2*Norder,4*Norder-1));
    for (int j=0; j<2*Norder-1; ++j) dp(j) = ps(j+1)-ps(j);
    dp(2*Norder-1) = ps(0)-ps(2*Norder-1);
    //
    for (int j=0; j<2*Norder-2; ++j) ddp(j) = ps(j+2)-ps(j);
    ddp(2*Norder-2) = ps(0)-ps(2*Norder-2);
    ddp(2*Norder-1) = ps(1)-ps(2*Norder-1);
    //
    sort(ps.data(), ps.data()+ps.size());
    sort(dp.data(), dp.data()+dp.size());
    sort(ddp.data(), ddp.data()+ddp.size());
    //
    p_repeat[i] = fmt::format("{0};{1};{2}", fmt::join(ps.data(), ps.data()+ps.size(), ","), fmt::join(dp.data(), dp.data()+dp.size(), ","), fmt::join(ddp.data(), ddp.data()+ddp.size(), ","));
  }
}

void FindDifferent(const Array<unsigned char,2>& diags1, const Array<unsigned char,2>& diags2, Array<int,1>& equal, deque<int>& new_diags)
{
  if (diags1.extent(1)!=diags2.extent(1)) cout<<"ERROR The diags size should be the same to compare"<<endl;
  equal.resize(diags2.extent(0));
  new_diags.clear();
  equal=-1;
  int Norder = diags1.extent(1)/4;
  int two_n = Norder*2;
  vector<string> p_repeat1, p_repeat2;
  FindWhichPseudoAndJumps(diags1, p_repeat1);
  FindWhichPseudoAndJumps(diags2, p_repeat2);
#ifdef _DEBUG  
  cout<<"p_repeat1="<<p_repeat1<<" p_repeat2="<<p_repeat2<<endl;
#endif  
  double aver_size=0;
  vector<int> idx(Norder);
  Array<unsigned char,1> d2n(4*Norder);
  Array<unsigned char,1> d2m(2*Norder);
  for (int i2=0; i2<diags2.extent(0); ++i2){
    // first check if the same number of pseudoparticles appear in time evolution.
    // Only those diagrams with the same number of repetitions of pseudoparticles can be equal.
    deque<int> possible;
    for (int i1=0; i1<diags1.extent(0); ++i1)
      if (p_repeat2[i2]==p_repeat1[i1]) possible.push_back(i1); // these could possible be equal
    
    Array<unsigned char,1> d2(diags2(i2,Range::all()));
    bool found=false;
    for (int i1 : possible){
      Array<unsigned char,1> d1(diags1(i1,Range::all()));
      if (d2==d1){
	equal(i2)=i1;
	found=true;
#ifdef _DEBUG	
	if (found) cout<<"i2="<<i2<<" size="<<possible.size()<<" ir=0 i1="<<i1<<" d2="<<d2<<" d1="<<d1<<endl;
#endif	
	break;
      }
    }
    int ir=1;
    while (not found and ir<two_n){
      // this rotates pseudoparticles on the ring
      for(int j=0; j<two_n; ++j) d2n(two_n+j) = d2(two_n + (j-ir+two_n)%two_n);
      // this changes hybridization propagators on the ring
      for(int j=0; j<Norder; ++j){
	d2m(2*j  ) = (d2(2*j  )+ir)%two_n;
	d2m(2*j+1) = (d2(2*j+1)+ir)%two_n;
	if (d2m(2*j)>d2m(2*j+1)) swap(d2m(2*j),d2m(2*j+1));
      }
      // initialize original index locations
      iota(idx.begin(), idx.end(), 0);
      stable_sort(idx.begin(), idx.end(),
		  [&d2m](int i1, int i2) {return d2m(2*i1) < d2m(2*i2);});
      for (int j=0; j<Norder; ++j){
	d2n(2*j)   = d2m(2*idx[j]);
	d2n(2*j+1) = d2m(2*idx[j]+1);
      }
      for (int i1 : possible){
	Array<unsigned char,1> d1(diags1(i1,Range::all()));
	if (d2n==d1){
	  equal(i2)=i1;
	  found=true;
#ifdef _DEBUG	  
	  if (found) cout<<"i2="<<i2<<" size="<<possible.size()<<" ir="<<ir<<" i1="<<i1<<" d2="<<d2<<" d1="<<d1<<" d2n="<<d2n<<endl;
#endif	  
	  break;
	}
      }
      ir+=1;
    }
    aver_size += possible.size();
    if (not found) new_diags.push_back(i2);
#ifdef _DEBUG    
    if (not found) cout<<"WARN: i2="<<i2<<" d2="<<d2<<" did not match any diag in d1"<<endl;
#endif    
  }
  aver_size *= 1./diags2.extent(0);
#ifdef _DEBUG  
  cout<<"aver_size="<<aver_size<<endl;
  cout<<"equal="<<equal<<endl;
  cout<<"new_diags.size="<<new_diags.size()<<endl;
#endif
}
void FindDifferent(const deque<Array<unsigned char,1>>& diags1, const Array<unsigned char,2>& diags2, Array<int,1>& equal, deque<int>& new_diags)
{
  //#define _DEBUG  
  if (diags1[0].extent(0)!=diags2.extent(1)){
    (*mpi->out)<<"ERROR The diags size should be the same to compare"<<endl;
    (*mpi->out)<<"diags1[0].extent(0)="<<diags1[0].extent(0)<<" diags1[0]="<<diags1[0]<<" diags2.extent(1)="<<diags2.extent(1)<<endl;
  }
  equal.resize(diags2.extent(0));
  new_diags.clear();
  equal=-1;
  int Norder = diags1[0].extent(0)/4;
  int two_n = Norder*2;
  vector<string> p_repeat1, p_repeat2;
  FindWhichPseudoAndJumps(diags1, p_repeat1);
  FindWhichPseudoAndJumps(diags2, p_repeat2);
#ifdef _DEBUG  
  (*mpi->out)<<"p_repeat1="<<p_repeat1<<endl;
  (*mpi->out)<<"p_repeat2="<<p_repeat2<<endl;
#endif  
  double aver_size=0;
  vector<int> idx(Norder);
  Array<unsigned char,1> d2n(4*Norder);
  Array<unsigned char,1> d2m(2*Norder);
  for (int i2=0; i2<diags2.extent(0); ++i2){
    // first check if the same number of pseudoparticles appear in time evolution.
    // Only those diagrams with the same number of repetitions of pseudoparticles can be equal.
    deque<int> possible;
    for (int i1=0; i1<diags1.size(); ++i1)
      if (p_repeat2[i2]==p_repeat1[i1]) possible.push_back(i1); // these could possible be equal
    
    Array<unsigned char,1> d2(diags2(i2,Range::all()));
    bool found=false;
    for (int i1 : possible){
      const Array<unsigned char,1>& d1=diags1[i1];
      if (d2==d1){
	equal(i2)=i1;
	found=true;
#ifdef _DEBUG	
	if (found) (*mpi->out)<<"i2="<<i2<<" size="<<possible.size()<<" ir=0 i1="<<i1<<" d2="<<d2<<" d1="<<d1<<endl;
#endif	
	break;
      }
    }
    int ir=1;
    while (not found and ir<two_n){
      // this rotates pseudoparticles on the ring
      for(int j=0; j<two_n; ++j) d2n(two_n+j) = d2(two_n + (j-ir+two_n)%two_n);
      // this changes hybridization propagators on the ring
      for(int j=0; j<Norder; ++j){
	d2m(2*j  ) = (d2(2*j  )+ir)%two_n;
	d2m(2*j+1) = (d2(2*j+1)+ir)%two_n;
	if (d2m(2*j)>d2m(2*j+1)) swap(d2m(2*j),d2m(2*j+1));
      }
      // initialize original index locations
      iota(idx.begin(), idx.end(), 0);
      stable_sort(idx.begin(), idx.end(),
		  [&d2m](int i1, int i2) {return d2m(2*i1) < d2m(2*i2);});
      for (int j=0; j<Norder; ++j){
	d2n(2*j)   = d2m(2*idx[j]);
	d2n(2*j+1) = d2m(2*idx[j]+1);
      }
      for (int i1 : possible){
	const Array<unsigned char,1>& d1=diags1[i1];
	if (d2n==d1){
	  equal(i2)=i1;
	  found=true;
#ifdef _DEBUG	  
	  if (found) (*mpi->out)<<"i2="<<i2<<" size="<<possible.size()<<" ir="<<ir<<" i1="<<i1<<" d2="<<d2<<" d1="<<d1<<" d2n="<<d2n<<endl;
#endif	  
	  break;
	}
      }
      ir+=1;
    }
    aver_size += possible.size();
    if (not found) new_diags.push_back(i2);
#ifdef _DEBUG    
    if (not found) (*mpi->out)<<"WARN: i2="<<i2<<" d2="<<d2<<" did not match any diag in d1"<<endl;
#endif    
  }
  aver_size *= 1./diags2.extent(0);
#ifdef _DEBUG  
  (*mpi->out)<<"aver_size="<<aver_size<<endl;
  (*mpi->out)<<"equal="<<equal<<endl;
#endif
}

void Unique_Diagrams_From_all_nodes(const deque<string>& accumulated_new_diags2, int N_min_order, int N_max_order,
				    vector<deque<Array<unsigned char,1>>>& combined2)
{//  Combines diagrams across all MPI nodes to set of unique diagrams.
  //    combines2[iorder][i](4*Norder)
  //
  deque<Array<unsigned char,2>> diagsv2;
  Convert_Diags_2_Array(accumulated_new_diags2, N_min_order, N_max_order, diagsv2);
  Array<int,1> equal;
  deque<int> new_diags;
  //combined2.resize(diagsv2.size());
  if (combined2.size()!=N_max_order-N_min_order+1) combined2.resize(N_max_order-N_min_order+1);
  for (int i=0; i<diagsv2.size(); ++i){// over diagrams of the same order
    if (mpi->rank==mpi->master){
      // Now we prepare diagrams from master as the first part of the output
      for (int j=0; j<diagsv2[i].extent(0); ++j) combined2[i].push_back( diagsv2[i](j,Range::all()) );
    }
#ifdef _MPI    
    int four_Norder=diagsv2[i].extent(1); // We have diag and pseudo combined, hence 4*Norder numbers
    int Norder = four_Norder/4;
    // We first sent diagrams accumulated on master to all other nodes.
    // All other notes will checked which diagrams they have which are not present on master.
    int sizes[] = {diagsv2[i].extent(0),diagsv2[i].extent(1)}; // sizes of diagrams to be broadcasted
    MPI_Bcast(sizes, 2, MPI_INT, mpi->master, MPI_COMM_WORLD);
    Array<unsigned char,2>* p_diagsv1; // pointer to the diagrams on master
    if (mpi->rank==mpi->master)
      p_diagsv1 = &(diagsv2[i]); // here it is pointer to data to be sent
    else // here is new array, which is a copy of diagrams in master
      p_diagsv1 = new Array<unsigned char,2>(sizes[0],sizes[1]); 
    // Now sending all diagrams from master to all other notes.
    MPI_Bcast(p_diagsv1->data(), sizes[0]*sizes[1], MPI_CHAR, mpi->master, MPI_COMM_WORLD);
    if (mpi->rank!=mpi->master){
      // Now we check how many diagrams were found on other nodes, which are different from master node.
      FindDifferent(*p_diagsv1,diagsv2[i],equal,new_diags); // new diags are those that are not found on master
      // new_diags is only index array, but now we create actual copy of these diagrams, and sent them back to master
      Array<unsigned char,2> new_d(new_diags.size(),four_Norder);
      for (int ii=0; ii<new_diags.size(); ++ii)
	new_d(ii,Range::all()) = diagsv2[i](new_diags[ii],Range::all());
      // Now send sizes for these new diagrams found on this node
      int size=new_diags.size()*four_Norder;
      Array<int,1> sizes(mpi->size); // We will use Allgather, so that we know when there is no new diagram found anywhere.
      MPI_Allgather(&size,1,MPI_INT,sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
      int total_size = sum(sizes)/four_Norder;// On each node we know how many new diagrams are found on all nodes combined.
      if (total_size>0){
#ifdef _DEBUG	
	(*mpi->out)<<Norder<<" size="<<size<<" realsize="<<new_d.size()<<endl<<" new_d="<<endl;
	(*mpi->out)<<new_d<<endl;
#endif	
	// Now sending the different diagrams found on this node.
	MPI_Gatherv(new_d.data(),size,MPI_CHAR,NULL,NULL,NULL,MPI_CHAR,mpi->master,MPI_COMM_WORLD);
      }
    }else{
      int size=0; // In the first part, nothing needs to be done on master, as we use diagrams for master to compare
      Array<int,1> sizes(mpi->size); // hence by definition there is no new diagrams on master.
      MPI_Allgather(&size,1,MPI_INT,sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
      // How many different/new diagrams are found on other nodes?
      int total_size = sum(sizes)/four_Norder;
      if (total_size>0){ // if none, we are done
	Array<int,1> rank_sizes(sizes/four_Norder); // similar than sizes, except counts diagrams, not amount of data transfered
	// for allgatherv we need to know displacements for data gathered.
	Array<int,1> displs(mpi->size);
	displs(0)=0;
	for (int i=1; i<mpi->size; ++i) displs(i)=displs(i-1)+sizes(i-1);
	// this will contain all gathered diagrams which were found to be new on any other node.
	Array<unsigned char,2> new_d(total_size,four_Norder);
#ifdef _DEBUG	
	(*mpi->out)<<Norder<<" sizes="<<rank_sizes<<" total_size="<<total_size<<" mpi_sizes="<<sizes<<" displs="<<displs<<" size(new_d)="<<new_d.size()<<endl;
#endif	
	// Now gathering all new diagrams on master
	MPI_Gatherv(new_d.data(),0,MPI_CHAR,new_d.data(),sizes.data(),displs.data(),MPI_CHAR,mpi->master, MPI_COMM_WORLD);
	//
#ifdef _DEBUG	
	(*mpi->out)<<new_d<<endl;
#endif	
	// next we find which node has maximum number of new diagrams. We than compare other sets with this maximal set.
	// Note that maximal set has to be nonempty, and total_size>0.
	int argmax=0; // node with maximal number of new diagrams
	for (int i=1; i<sizes.extent(0); ++i)
	  if (sizes(i)>sizes(argmax)) argmax=i;
	// displacements in terms of number of diagrams in array new_d, which contains all new diagrams from all other nodes.
	Array<int,1> disp(mpi->size+1);
	disp(0)=0;
	for (int i=1; i<mpi->size+1; ++i) disp(i)=disp(i-1)+rank_sizes(i-1);
	// In terms of displacements, we know that diagrams in the interval [is1,ie1] should be used for comparison.
	int is1 = disp(argmax), ie1=disp(argmax+1);
	//
#ifdef _DEBUG	
	(*mpi->out)<<"argmax="<<argmax<<" is1="<<is1<<" ie1="<<ie1<<endl;
#endif	
	// currents will contain all unique new diagrams. Namely, accumuated diagrams in new_d are different from those on master (*p_diagsv1)
	// but can be equivalent. Hence now we need to find unique diagrams obtained by Gatherv.
	deque<Array<unsigned char,1>> currents;
	// At the start, the maximal node will be added to currents. These diags are definitely non-equivalent to those in *p_diagsv1
	for (int i=is1; i<ie1; ++i) currents.push_back( new_d(i,Range::all()) );
#ifdef _DEBUG	
	(*mpi->out)<<"currents are "<<endl;
	for (int i=0; i<currents.size(); ++i) (*mpi->out)<<currents[i]<<endl;
	(*mpi->out)<<endl;
#endif	
	// Now we check diagrams obtained from node [1...mpi.size] and compare them to those on argmax node.
	for (int ik=1; ik<mpi->size; ++ik){
	  if (ik==argmax) continue; // of course argmax node should be skept, because we already accepted these diagrams.
	  int is2 = disp(ik), ie2 = disp(ik+1); // where in array new_d those diagrams appear?
	  Array<unsigned char,2> nd(new_d(Range(is2,ie2-1),Range::all())); // These are diagrams from node k
#ifdef _DEBUG	  
	  (*mpi->out)<<ik<<") considering nds is2="<<is2<<" ie2="<<ie2<<endl;
	  for (int i=0; i<nd.extent(0); ++i) (*mpi->out)<<nd(i,Range::all())<<endl;
#endif	  
	  // Do we have any different diagram from node k than those already accmulated?
	  FindDifferent(currents,nd,equal,new_diags);
#ifdef _DEBUG	  
	  (*mpi->out)<<ik<<") Found equal="<<equal<<" and new_diags=";
	  if (new_diags.size()>0) (*mpi->out)<<new_diags;
	  (*mpi->out)<<endl;
#endif	  
	  // Now we add the diagrams which were found to be unique to the currents. These are to be added to *p_diagsv1 from master
	  for (int i=0; i<new_diags.size(); ++i) currents.push_back( new_d(is2+new_diags[i],Range::all()) );
#ifdef _DEBUG
	  (*mpi->out)<<"currents are "<<endl;
	  for (int i=0; i<currents.size(); ++i) (*mpi->out)<<currents[i]<<endl;
	  (*mpi->out)<<endl;
#endif	  
	}
	// Now we add the unique new diagrams
	combined2[i].insert(combined2[i].end(), currents.begin(), currents.end());
      }
    }
    if (mpi->rank!=mpi->master) delete p_diagsv1;
#endif    
  }
}

auto combined2diag(Array<unsigned char,1>& cb)
{
  struct retVals {
    Array<int,2> diag;
    Array<int,1> pseud;
  };
  int norder = cb.extent(0)/4;
  retVals r;
  r.diag.resize(norder,2);
  r.pseud.resize(2*norder);
  for (int i=0; i<norder; ++i){
    r.diag(i,0) = cb(2*i);
    r.diag(i,1) = cb(2*i+1);
  }
  for (int i=0; i<2*norder; ++i) r.pseud(i) = cb(2*norder+i);
  return r;
}

/*
int main(int argc, char *argv[])
{
  setvbuf (stdout, NULL, _IONBF, BUFSIZ);
  mpi = new MPIcls(argc, argv);

  int N_min_order=2, N_max_order=6;
  
  deque<string> accumulated_new_diags2 = ReadDiags("debug_diags."s+to_string(mpi->rank));

  vector<deque<Array<unsigned char,1>>> combined2;
  Unique_Diagrams_From_all_nodes(accumulated_new_diags2, N_min_order, N_max_order, combined2);
  
  if (mpi->rank==mpi->master){
    for (int i=0; i<combined2.size(); ++i){
      (*mpi->out)<<"Norder="<< i+N_min_order<<endl;
      for (int j=0; j<combined2[i].size(); ++j)
	(*mpi->out)<<j<<" "<<combined2[i][j]<<endl;
    }
  }
  delete mpi;
  return 0;
}
*/
