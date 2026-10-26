#ifndef _CIX_
#define _CIX_
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
#include <map>
#include <unordered_map>
#include <string>
#include <blitz/array.h>
#include "basics.h"

using namespace blitz;
using namespace std;

void inline_trim(std::string& s);
struct NjmResult{
  int njm_c;
  double njm_r;
  const Array<double,2>* njm;   // pointer to array
  bool transposed;        // whether it need to be transposed
  NjmResult() : njm_c(0), njm_r(0), njm(NULL), transposed(false){}
};
struct BathBlocks{
  std::vector<std::vector<int>> blocks;   // blocks[b] = list of orbitals in block b
  std::vector<bool> is_block_sc;          // is_sc_block[b] = true if block b has SC (Nambu) terms
  std::vector<TinyVector<int,3>> orb_ind; // orb_ind[Nbath] = {ib,j,k}
};
inline std::ostream& operator<<(std::ostream& os, const BathBlocks& bb)
{
  os << "BathBlocks\n";
  os << "==========\n";
  // Blocks
  os << "blocks (" << bb.blocks.size() << "):\n";
  for (std::size_t b = 0; b < bb.blocks.size(); ++b) {
    os << "  block " << b << " : [ ";
    for (std::size_t i = 0; i < bb.blocks[b].size(); ++i) {
      os << bb.blocks[b][i];
      if (i + 1 < bb.blocks[b].size()) os << ", ";
    }
    os << " ]\n";
  }
  // SC flags
  os << "is_block_sc (" << bb.is_block_sc.size() << "):\n";
  for (std::size_t b = 0; b < bb.is_block_sc.size(); ++b) {
    os << "  block " << b << " : "
       << (bb.is_block_sc[b] ? "SC" : "normal") << "\n";
  }
  // Orbital indices
  os << "orb_ind (" << bb.orb_ind.size() << "):\n";
  for (std::size_t i = 0; i < bb.orb_ind.size(); ++i) {
    os << "  orb " << i << " : { "
       << bb.orb_ind[i][0] << ", "
       << bb.orb_ind[i][1] << ", "
       << bb.orb_ind[i][2] << " }\n";
  }
  return os;
}

class Cix{
public:
  bool kugler_estimators=false; // Optional analytic physical endpoint correlators.
  string cix_filename; // Also identifies an optional <cix>.nlsj sidecar.
  int Nsite, Natst, Npsi, Nbath, Nblock, Nmaxsize; // cluster size, number of cluster states, number of creation operators, number of hybridization entries, number of Hamiltonian blocks, the largest dimension of Hamiltonian
  Array<int,1> bath_sigind0, bath_sigind;//, bath_gflip;
  Array<double,1> bath_spin;
  Array<bool,1> bath_SC, SC_bath;
  Array<bool,1> bath_conj;
  //Array<int,2> bath_psip_psi;
  Array<string,1> bath_leg;
  Array<double,1> Ebath;
  Array<int,1> iBlck;
  Array<int,1> iN;
  Array<double,1> iSz;
  Array<int,2> Fpsi, Fi;
  Array<int,2> Fbpsi, Fbc;
  Array<int,1> Fi_indx;
  vector<Array<double,2>> Ham;
  vector<Array<int,1>> Block, nmBlock;
  Array<int,1> inBlock;
  vector<vector<Array<double,2>>> Njm;
  Array<int,2> Njm_c;
  unordered_map<int,double> Njm_r;
  vector<Array<int,2>> Sgm_symmetry;
  vector<deque<TinyVector<int,3>>> Sgm_indx;
  int Npseudo_unique, Nbath_unique, Nbath_unique_unique;
  double Emin, Emax, F_atom;
  double F0_input = 0.0, mu_input = 0.0; // Retain the CIX energy correction for Probability2.dat.
  Array<int,1> deg_pp;
  Array<int,1> deg_bth;
  Array<int,1> deg_bath;
  vector<vector<int>> unique_blocks;
  int first_anomalous;
  //
  //  pp = cix.Fbc(op(i),pp);
  //  pp = cix.Fi(op(i),pp);
  //
  //   ih=1...n
  //     int i = diag(ih,0);
  //     int j = diag(ih,1);
  //     int ds = dirspin(ih);
  //     op(j) = ds;
  //     op(i) = 1-(ds & 1) + 2*(ds>>1);  // operator opposite to ds: from number (2n+1) gives 2n, and from 2n gives (2n+1)
  // 
  //
  bool parse_bath_entry(std::string& line,
			int& ii, int& psi_dagger, int& psi, int& op_int, bool& op_conj, double& spin_value, std::string& moment_legend)
  {
    if (auto pos = line.find('#'); pos != std::string::npos) line.erase(pos);   // drop everything after '#'
    std::istringstream iss(line);
    std::string ops[3];  // token like "-0*"
    // 1) Read five whitespace-separated tokens
    if (!(iss >> ii >> ops[0] >> ops[1] >> ops[2] >> spin_value)){
      return false; // parsing failed
    }
    std::getline(iss, moment_legend);
    // Remove leading whitespace from moment_legend
    moment_legend.erase(0, moment_legend.find_first_not_of(" \t"));
    
    int ops_ints[3]={0,0,0};
    int ops_conj[3]={0,0,0};
    for (int i=0; i<3; ++i){
      // 2) Handle the 1-3th token: integer with optional '-' and optional '*'
      ops_conj[i] = 0;
      if (!ops[i].empty() && ops[i].back() == '*') {
        ops_conj[i] = 1;
        ops[i].pop_back();      // remove the trailing '*'
      }
      // Now op_str should be something like "0", "-1", "3", etc.
      ops_ints[i] = std::stoi(ops[i]); // std::stoi handles leading '-'
    }
    psi_dagger = 2*ops_ints[0]+ops_conj[0];
    psi        = 2*ops_ints[1]+ops_conj[1];
    op_int = ops_ints[2];
    op_conj= ops_conj[2];
    return true;
  }

  int encode_block_bath(int ib, int ia) const { return ia+Nbath*ib; } // second index in bath index

  NjmResult Njm_ind(int ib, int two_ia) const
  {
    bool creation= (Fi_indx(two_ia) % 2==0); // Is Fi(two_ia) is creation operator
    int ia_op = Fi_indx(two_ia)/2;           // which operator is Fi(two_ia) in the list psi^+_{ia_op}
    int ibs = creation ? ib : Fbpsi(2*ia_op+1,ib); //Fbc(2*ia_op+1,ib); // if creation, Njm is already stored as Njm[ib][ia_op]. If anhilation, we should find the block for which creation operator brings us to ib
    int ibn = creation ? Fbpsi(2*ia_op,ib) : ib;
    NjmResult nx;
    nx.njm_c = Njm_c(ibs,ia_op);
    //if (not (nx.njm_c==0 or nx.njm_c==1 or nx.njm_c==2) ) cout<<"ERROR nx.njm_c set to "<<nx.njm_c<<endl;
    nx.transposed = not creation;
    if (nx.njm_c==1){
      //nx.njm_r = Njm_r[encode_block_bath(ibs,ia_op)];
      nx.njm_r = Njm_r.at(encode_block_bath(ibs,ia_op));
    }else if (nx.njm_c==2){
      nx.njm = &Njm[ibs][ia_op];
    }
    return nx;
  }

  void Read_cix(const string& _fcix_, double F0, double mu, double beta)
  {
    cix_filename = _fcix_;
    F0_input = F0;
    mu_input = mu;
    bool DEBUG=false;
    ifstream fcix(_fcix_);
    if (not fcix.good()){
      cout<<"ERROR cix file '"<<_fcix_<<"' does not exists or is corrupt."<<endl;
      exit(1);
    }
    Npseudo_unique=0; // This is only determined through Diag2Matrix.Determine_Symmetry
    string line;
    getline(fcix, line);
    getline(fcix, line);
    fcix >>Nsite>>Natst>>Npsi>>Nbath>>Nblock>>Nmaxsize;
    if (not fcix or fcix.eof() or Nbath<=0 or Natst<=0){
      cout<<"ERROR reading cix file '"<<_fcix_<<"' first line"<<endl;
      exit(1);
    }
    //cout<<"Nsite="<<Nsite<<" Nmaxsize="<<Nmaxsize<<endl;
    getline(fcix, line); // finish line
    getline(fcix, line); // skip next line
    //bath_gflip.resize(Nbath);
    bath_sigind0.resize(Nbath);
    bath_sigind.resize(Nbath); bath_sigind=-1;
    bath_conj.resize(Nbath); bath_conj=false;
    bath_spin.resize(Nbath);
    bath_SC.resize(Nbath);
    SC_bath.resize(Nbath);
    bath_leg.resize(Nbath);
    //bath_psip_psi.resize(Nbath,2);
    Fi_indx.resize(2*Nbath);
    int i_bath_last=0;

    for (int i = 0; i < Nbath; ++i) {
      int ii, psi_dagger, psi, op_int;
      bool op_conj;
      double spin_value;
      string line, leg;
      getline(fcix, line);
      bool pars_ok = parse_bath_entry(line,
				      ii, psi_dagger, psi, op_int, op_conj, spin_value,leg);
      //fcix >>ii>>idim>>gflp>>isig>>spin_value;
      //getline(fcix,leg);
      //inline_trim(leg);
      if (ii!=i or not pars_ok) cout<<"WARNING reading cix="<<_fcix_<<" i="<<i<<"!="<<ii<<" or parsing of baths failed"<<endl;
      //cout<<"ii="<<ii<<" dim="<<idim<<" gflp="<<gflp<<" isig="<<isig<<" spin="<<isp<<" leg='"<<leg<<"'"<<endl;
      //bath_gflip(i) = gflp;
      //bath_psip_psi(i,0) = psi_dagger;
      //bath_psip_psi(i,1) = psi;
      Fi_indx(2*i  ) = psi_dagger;
      Fi_indx(2*i+1) = psi;
      bath_sigind0(i) = op_int;
      bath_conj(i) = op_conj;
      bath_spin(i) = spin_value;
      bath_leg(i) = leg;
      bath_SC(i) = (psi_dagger+psi)%2==0;  // if we have psi^+ psi it is fermionic, but psi^+ psi^+ and psi psi is SC and bosonic
      // determine sign from bath_sigind0(i).  Assume bath_sigind0(i) != 0.
      int sgn = (bath_sigind0(i) >= 0) ? +1 : -1;
      bool found = false;
      for (int j = 0; j < i; ++j) {
        // same abs index and same conjugation => same family
        if (bath_conj(i) == bath_conj(j) && bath_sigind0(i)==bath_sigind0(j)){
	  //int base = abs(bath_sigind(j));  // family label (positive)
	  bath_sigind(i) = bath_sigind(j);
	  found = true;
	  break;
        }
      }
      if (!found) {
        // new family for this (abs(index), conj) combination
        bath_sigind(i) = i_bath_last;
        ++i_bath_last;
      }
      /*
      for (int j=0; j<i; ++j){
	if (bath_sigind0(i)==bath_sigind0(j) and bath_conj(i)==bath_conj(j)){
	  bath_sigind(i)=bath_sigind(j);
	  break;
	}
      }
      if (bath_sigind(i)<0){
	bath_sigind(i) = i_bath_last;
	i_bath_last+=1;
      }
      */
    }
    Nbath_unique_unique = max(abs(bath_sigind0))+1;
    Nbath_unique = max(bath_sigind)+1;
    Ebath.resize(Nbath_unique_unique);
    deg_bath.resize(Nbath_unique_unique); deg_bath=0;
    getline(fcix, line); // skip
    //cout<<"before skipping line "<<line<<endl;
    for (int ie=0; ie<Nbath_unique_unique; ++ie) fcix >> Ebath(ie);
    for (int iorb=0; iorb<Nbath; ++iorb) deg_bath(abs(bath_sigind0(iorb))) += 1;

    SC_bath=false;
    for (int iorb = 0; iorb < Nbath; ++iorb){
      if (bath_SC(iorb)){
	for (int jorb=0; jorb<Nbath; ++jorb){
	  if (Fi_indx(2*jorb)==Fi_indx(2*iorb) or Fi_indx(2*jorb)==Fi_indx(2*iorb+1)) SC_bath(jorb)=true;
	  if (Fi_indx(2*jorb+1)==Fi_indx(2*iorb+1) or Fi_indx(2*jorb+1)==Fi_indx(2*iorb)) SC_bath(jorb)=true;
	}
      }
    }
    //Array<bool,1> _bath_SC_(Nbath_unique); _bath_SC_=false;
    //for (int iorb=0; iorb<Nbath; ++iorb) _bath_SC_(bath_sigind(iorb)) = SC_bath(iorb);
    //cout << "_bath_SC_="<< _bath_SC_ << endl;
    if (DEBUG){
      cout<< "Fi_indx="<< Fi_indx << endl;
      cout<< "sigind0="<< bath_sigind0 << endl;
      cout << "sigind=" << bath_sigind << endl;
      cout << "spin="<< bath_spin << endl;
      cout << "bath_conj="<< bath_conj << endl;
      cout << "Nbath_unique="<<Nbath_unique<<endl;
      cout << "Nbath_unique_unique="<< Nbath_unique_unique << endl;
      cout << "leg="<< bath_leg << endl;
      cout << "Ebath="<<Ebath<<endl;
      cout << "bath_SC=" << bath_SC << endl;
      cout << "SC_bath="<< SC_bath << endl;
    }

    //cout<<"bath_fermion="<<bath_fermion<<endl;
    std::getline(fcix, line); // finish line
    //cout<<"skipping line "<<line<<endl;
    std::getline(fcix, line); // skip
    //cout<<"skipping line "<<line<<endl;
    //cout<<"line="<<line<<endl;
    iBlck.resize(Natst+1);iBlck=0;
    iN.resize(Natst);
    iSz.resize(Natst);
    Array<int,1> iK(Natst);
    Fpsi.resize(2*Npsi,Natst+1);
    Fpsi=0; // These are psi^+_{ib} stored in Fpsi[2*ib,:], which are specified in cix file
    for (int i=0; i<Natst; ++i) {
      int ii;
      fcix >> ii>>iBlck(i+1)>>iN(i)>>iK(i)>>iSz(i);
      for (int ib=0; ib<Npsi; ++ib) fcix >> Fpsi(2*ib,i+1);
      if (i+1 != ii){
	cout<<"WARNING reading cix file "<<_fcix_<<" i="<<i<<" ii="<<ii<<" reading atomic states last line="<<line<<endl;
      }
      std::getline(fcix, line); // finish line
      //cout<<" Fi["<<i+1<<"]="<<Fi(Range::all(),i+1)<<endl;
    }
    // From Fpsi[2*ib,:]=psi^+_{ib} we compute psi_{ib}, which is stored in Fpsi[2*ib+1,:]
    for (int i=1; i<=Natst; ++i){
      for (int ib=0; ib<Npsi; ++ib){
	int ifi = Fpsi(2*ib,i);
	Fpsi(2*ib+1,abs(ifi)) = i*sign(ifi);
      }
    }
    Fi.resize(2*Nbath, Natst+1);
    for (int ib=0; ib<2*Nbath; ++ib){
      Fi(ib,Range::all()) = Fpsi(Fi_indx(ib),Range::all());
      if (DEBUG) cout<<"Setting Fi["<<ib<<",:]=Fpsi["<<Fi_indx(ib)<<",:]"<<endl;
    }
    std::getline(fcix, line); // skip
    //cout<<"skiping line"<<line<<endl;
    Ham.resize(Nblock);
    vector<Array<double,2>> Spn;
    Spn.resize(Nblock);
    Block.resize(Nblock);
    inBlock.resize(Natst+1);
    nmBlock.resize(Nblock);
    Emin=1000, Emax=-1000;
    double U = (Nsite==1) ? F0 : 0;
    double Zatom = 0;
    for (int i = 0; i < Nblock; ++i) {
      int ii, idim;
      fcix >> ii >> idim;
      //cout<<"ii="<<ii<<" and idim="<<idim<<" was read"<<endl;
      if (ii!=i+1) {
	cout<<"WARNING ii="<<ii<<" i="<<i+1<<" in reading Hamiltonian blocks"<<endl;
	if (i>0){
	  cout<<"last Ham="<<Ham[i-1]<<endl;
	  cout<<"last S2="<<Spn[i-1]<<endl;
	}
	DEBUG=true;
      }
      Block[i].resize(idim);
      for (int j=0; j<idim; ++j) fcix>>Block[i](j);
      for (int j=0; j<idim; ++j) inBlock(Block[i](j))=j;

      nmBlock[i].resize(idim);
      nmBlock[i](0)=1;
      for (int j=1; j<idim; ++j){
	if (Block[i](j)==Block[i](j-1)+1)
	  nmBlock[i](j)=nmBlock[i](j-1);
	else
	  nmBlock[i](j)=nmBlock[i](j-1)+1;
      }

      if (DEBUG){
	cout<<"Block["<<i<<"]="<<Block[i]<<" idim="<<idim<<endl;
	cout<<"nmBlock[i]="<<nmBlock[i]<<endl;
      }
      std::getline(fcix, line);
      // remove comments
      if (auto p=line.find('#'); p != string::npos) line.erase(p);
      if (DEBUG) cout<<"1.line="<<line<<endl;
      // 2) turn '[', ']', ',' into plain spaces
      replace_if(line.begin(), line.end(),
		 [](char c) { return c == '[' || c == ']' || c == ','; }, ' ');
      if (DEBUG) cout<<"2.line="<<line<<endl;
      std::istringstream iss(line);
      Array<double,1> energies(idim), spin(idim);
      spin=0;
      for (int j=0; j<idim; ++j) iss >> energies(j);
      bool spin_present=false;
      if (iss >> spin(0)){
	for (int j=1; j<idim; ++j) iss >> spin(j);
	spin_present=true;
      }
      if (DEBUG) cout<<"energies="<<energies<<" spin_present="<<spin_present<<endl;
      Ham[i].resize(idim,idim);
      Spn[i].resize(idim,idim); Spn[i]=0;
      for (int j = 0; j<idim; ++j) {
	for (int k = 0; k < idim; ++k) fcix>>Ham[i](j,k);
	if (spin_present){
	  for (int k = 0; k < idim; ++k) fcix>>Spn[i](j,k);
	}
      }
      if (DEBUG) cout<<"Ham["<<i<<"]="<<Ham[i]<<endl;
      // correct energies as should be corrected in the main part of the code.
      for (int j = 0; j<idim; ++j) {
	double N = iN(Block[i](j)-1);
	double dE = U*N*(N-1)/2. - mu*N;
	//cout<<"Correct E0["<<i<<","<<j<<"]="<<energies(j)<<" N="<<N<<" dE="<<dE<<" E="<<energies(j)+dE<<endl;
	energies(j) += dE;
      }
      double Emin_ = min(energies);
      double Emax_ = max(energies);
      if (Emin_ < Emin) Emin = Emin_;
      if (Emax_ > Emax) Emax = Emax_;
      for (int k=0; k<energies.size(); k++) Zatom += exp(-beta*energies(k));
    }
    
    F_atom = -log(Zatom)/beta;
    //if (mpi.rank==mpi.master) (*mpi.out) << "log(Zatom)=" << log(Zatom) << " F_atom="<<F_atom<<endl;
    
    std::getline(fcix, line); // to the end of previous line
    std::getline(fcix, line); // skip
    std::unordered_map<int,double> Uc;
    while(fcix.good() and not fcix.eof()){
      TinyVector<int,4> ms; double Ur;
      fcix >> ms[0] >> ms[1] >> ms[2] >> ms[3] >> Ur;
      std::getline(fcix, line); // to the end if comments
      if (not fcix.good()) break;
      Uc[encode(ms)]=Ur;
    }
    if (Nsite==1){
      for (int a=0; a<Npsi; a++)
	for (int b=0; b<Npsi; b++)
	  if (a!=b) Uc[encode(TinyVector<int,4>(a,b,b,a))]+=F0;
    }
    if (DEBUG){
      /*
      for (const auto& [key,value] : Uc){
	TinyVector<int,4> ms = decode(key);
	cout<<"Uc"<<ms<<"="<<value<<endl;
      }
      */
      // it is more useful to see sorted Uc
      std::vector<std::pair<int,double>> items(Uc.begin(), Uc.end());
      std::sort(items.begin(), items.end(),
		[](const auto& a, const auto& b) {
		  return a.first < b.first;
		});
      for (const auto& [key, value] : items) {
	TinyVector<int,4> ms = decode(key);
	cout << "Uc" << ms << "=" << value << endl;
      }
    }
    int maxBlockSize=0;
    for (int i=0; i<Nblock; ++i)
      if (Block[i].size()>maxBlockSize) maxBlockSize=Block[i].size();
    if (Nmaxsize!=maxBlockSize){
      cout<<"WARNING cix.Nmaxsize="<<Nmaxsize<<" while actual Nmaxsize="<<maxBlockSize<<" correcting it..."<<endl;
      Nmaxsize=maxBlockSize;
    }
    if (DEBUG){
      // printing out the content of cix file
      for (int i=0; i<Natst; ++i) {
	cout << setw(3) <<i+1<<" "<<setw(2)<<iBlck(i+1)<<" "<<setw(2)<<iN(i)<<" "<<setw(2)<<iK(i)<<" "<<setw(4)<<iSz(i);
	for (int j=0; j<Nbath; ++j) cout << setw(3) << Fi(2*j,i+1) <<" ";
	cout<<endl;
      }
      cout<<"# Hamiltonians"<<endl;
      for (int i=0; i<Nblock; ++i){
	cout<<" Block["<<i<<"]="<<Block[i]<<endl;
	cout<<"H["<<setw(3)<<i<<"]=";
	for (int j=0; j<Ham[i].extent(0); ++j){
	  if (j>0) cout<<"       ";
	  for (int k=0; k<Ham[i].extent(0); ++k) cout<<setw(8)<<Ham[i](j,k)<<" ";
	  cout<<endl;
	}
      }
    }
    // This is psi^+_{ib} and psi_{ib} for blocks but for Fpsi specified in cix file
    Fbpsi.resize(2*Npsi,Nblock);
    Fbpsi = -1;
    for (int i=0; i<Nblock; ++i){
      for (int ib=0; ib<Npsi; ++ib){
	for (int j=0; j<Block[i].size(); ++j){
	  int ifi = Fpsi(2*ib,Block[i](j));
	  if (ifi!=0){
	    int jj = iBlck(abs(ifi))-1;
	    if (Fbpsi(2*ib,i)==-1){
	      Fbpsi(2*ib,i) = jj;
	      Fbpsi(2*ib+1,jj) = i;
	    }else{
	      if (Fbpsi(2*ib,i) != jj){
		cout<<"WARNING psi^+_{"<<ib<<"}|iblock="<<i+1<<">=|"<<"iblock="<<jj+1<<"> but also|iblock="<<Fbpsi(2*ib,i)+1<<">"<<endl;
		cout<<"Your cix file has likely an error"<<endl;
	      }
	    }
	  }
	}
      }
    }
    // This is psi^+_{ib} and psi_{ib} for blocks but for psis in all hybriduzation baths
    Fbc.resize(2*Nbath,Nblock);
    Fbc = -1;
    for (int i=0; i<Nblock; ++i){
      for (int ib=0; ib<Nbath; ++ib){
	for (int j=0; j<Block[i].size(); ++j){
	  int ifi = Fi(2*ib,Block[i](j));
	  if (ifi!=0){
	    int jj = iBlck(abs(ifi))-1;
	    if (Fbc(2*ib,i)==-1){
	      Fbc(2*ib,i) = jj;
	      Fbc(2*ib+1,jj) = i;
	    }else{
	      if (Fbc(2*ib,i) != jj){
		cout<<"WARNING psi^+_{"<<ib<<"}|iblock="<<i+1<<">=|"<<"iblock="<<jj<<"> but also|iblock="<<Fbc(2*ib,i)<<">"<<endl;
		cout<<"Your cix file has likely an error"<<endl;
	      }
	    }
	    //break;
	  }
	}
      }
    }
    if (DEBUG){
      for (int i=0; i<Nblock; ++i){
	for (int ib=0; ib<Nbath; ++ib){
	  if (Fbc(2*ib,i)>=0)
	    cout << "Psi^+_{"<<ib<<"} |"<<i<<"> -> |"<<Fbc(2*ib,i)<<">"<<endl;
	}
      }
      for (int i=0; i<Nblock; ++i){
	for (int ib=0; ib<Nbath; ++ib){
	  if (Fbc(2*ib+1,i)>=0)
	    cout << "Psi_{"<<ib<<"} |"<<i<<"> -> |"<<Fbc(2*ib+1,i)<<">"<<endl;
	}
      }
    }
    
    for (int i=1; i<=Natst; ++i){
      int ib = iBlck(i)-1;         // which block
      int in_block = inBlock(i);   // where in the block
      int ii = Block[ib](in_block);// global index for state in the block
      if (ii!=i){
	cout<<"ERROR expecting the entries in blocks to be consistent.";
	cout<<" state i="<<i<<" is in block ib="<<ib<<" and its place inside the block is "<<in_block<<" which holds global index "<<ii<<" which is different from "<<i<<endl;
      }
      //cout<<i<<" "<<ii<<" Blck="<<ib<<" inBlock="<<in_block<<endl;
    }
    // Remember that Njm is transposed as compared to old ctqmc code.
    bool Print=DEBUG;
    Njm_c.resize(Nblock,Npsi); Njm_c=0; //Njm_c.resize(Nblock,Nbath); Njm_c=0;
    Njm.resize(Nblock);
    for (int ib=0; ib<Nblock; ++ib){ // over all blocks of states
      Njm[ib].resize(Npsi); //Njm[ib].resize(Nbath);
      if (Print) cout<<" Block="<< Block[ib] << endl;
      for (int ia=0; ia<Npsi; ++ia){// over all baths, and expect Njm to be proportional to psi^_{ia}|state0>
	//
	bool creation= true; //(Fi_indx(2*ia) % 2==0); // means Fi(2*ia,:) is creation operator
	int ia_op = ia;      // Fi_indx(2*ia)/2;
	//int st_bl =  Fbc(2*ia,ib);    // |st_bl> = psi^+_{ia}|block0>
	int st_bl =  Fbpsi(2*ia,ib);    // |st_bl> = psi^+_{ia}|block0>
	//                              For operators we have Fi[2*ia,:]=Fpsi[2*ia_op,:], which is at ia_op-position in the list of specified operators,
	//                              and Fi[2*ia,:] is Fpsi[Fi_indx[2*ia],:], hence 2*ia_op=Fi_indx(2*ia)
	if (Print) cout<<"ib="<<ib<<" ia="<<ia<<" Fi_indx["<<ia<<"]="<<ia_op<<" creation="<<creation<<" st_bl="<<st_bl<<endl;
	if (st_bl==-1) continue;      // if psi^+_{ia}|block0>=0 we know that Njm should also vanish
	if (Print) cout<<"ib="<<ib<<" ia="<<ia<<" Block["<<ib<<"].size="<<Block[ib].size()<<" Block["<<st_bl<<"].size="<<Block[st_bl].size()<<endl;
	Array<double,2> njm(Block[ib].size(),Block[st_bl].size()); njm=0;// This is expected size of njm in going from block ib to st_bl
	for (int i_st=0; i_st<Block[ib].size(); ++i_st){// now we go over all atomic states in block ib
	  int i=Block[ib](i_st); // i is state |state0>==|block=ib,state=i>
	  if (Print) cout<<"   ia="<<ia<<" i="<<i<<" i_st="<<i_st<<endl;
	  for (int ik=0; ik<Npsi; ++ik){
	    int i1 = Fpsi(2*ik+1,i);        // |i1> = psi_{ik}|state0>. Here we rely on the first Npsi operators to be psi^dagger operators
	    if (i1==0) continue;
	    for (int ij=0; ij<Npsi; ++ij){
	      int i2 = creation ? Fpsi(2*ij,abs(i1)) : Fpsi(2*ij+1,abs(i1));  // |i2> = psi^+_{ij}|i1> or |i2>=psi_{ij}|i1> if anhilation
	      if (i2==0) continue;
	      for (int ii=0; ii<Npsi; ++ii){
		int i3 = Fpsi(2*ii,abs(i2)); // |i3> = psi^+_{ii}|i2>
		if (i3==0) continue;
		int key1 = creation ? encode(TinyVector<int,4>(ii,ij,ik,ia_op)) : encode(TinyVector<int,4>(ik,ij,ii,ia_op));
		int key2 = creation ? encode(TinyVector<int,4>(ij,ii,ik,ia_op)) : encode(TinyVector<int,4>(ij,ik,ii,ia_op));
		double UU=0;
		if (Uc.contains(key1)) UU+=0.5*Uc[key1];
		if (Uc.contains(key2)) UU-=0.5*Uc[key2];
		if (UU!=0){
		  int ib3 = iBlck(abs(i3));
		  int tots = sign(i1)*sign(i2)*sign(i3);
		  if (st_bl+1==ib3){
		    int n2 = inBlock(abs(i3));
		    njm(i_st,n2) += UU*tots;
		    if (false and ia>=8){
		      int ib3 = iBlck(abs(i3));
		      int tots = sign(i1)*sign(i2)*sign(i3);
		      cout<<"   ia="<<ia<<"  i="<<i_st<<" j="<<n2<<"="<<abs(i3)<<" njm="<<UU*tots<<" psi^+_{"<<ii<<"}psi^+_{"<<ij<<"}psi_{"<<ik<<"}|"<<i<<">="<<tots<<"|"<<abs(i3)<<"> ";
		      cout<<"Uc["<<ii<<","<<ij<<","<<ik<<","<<ia_op<<"]="<<UU<<" st_bl="<<st_bl+1<<" ib3="<<ib3<<" Fi[i]="<<Fi(2*ia,i)<<": i3="<<i3<<",i2="<<i2<<",i1="<<i1<<endl;
		    }
		  }else{
		    cout<<"WARNING mode S is likely not working with this cix-file. Check that G mode gives the same answer!"<<endl;
		    cout<<"WARN : Probably UCoulomb in cix file is wrong : ib="<<ib<<" i="<<i<<" b1="<<st_bl<<" i3="<<i3<<" ib3="<<ib3<<endl;
		    cout<<"WARN : ia="<<ia<<"  i="<<i_st<<" psi^+_{"<<ii<<"}psi^+_{"<<ij<<"}psi_{"<<ik<<"}|"<<i<<">="<<tots<<"|"<<abs(i3)<<"> ";
		    cout<<"WARN : Uc["<<key1<<"="<<ii<<","<<ij<<","<<ik<<","<<ia_op<<"]="<<UU<<" st_bl="<<st_bl+1<<" ib3="<<ib3<<" Fi[i]="<<Fi(2*ia,i)<<": i3="<<i3<<",i2="<<i2<<",i1="<<i1<<endl;
		    cout<<"WARN : psi_{"<<ik<<"}|"<<i<<">=|"<<i1<<">"<<endl;
		    cout<<"WARN : psi^+_{"<<ij<<"}|"<<i1<<">=|"<<i2<<">"<<endl;
		    cout<<"WARN : psi^+_{"<<ii<<"}|"<<i2<<">=|"<<i3<<">"<<endl;
		  }
		}
	      }
	    }
	  }
	}
	if (sum(abs(njm))<1e-13){
	  Njm_c(ib,ia)=0;
	  if (Print) cout<<"Njm["<<setw(3)<<ib+1<<","<<ia<<"]=0"<<endl;
	}else{
	  bool Q_proportional=true;
	  double r=1;
	  for (int i_st=0; i_st<Block[ib].size(); ++i_st){ // normally this loop should use only i_st=0 and exit, unless Fi gives 0, and we cannot compute ratio r.
	    int i4 = Fi(2*ia,Block[ib](i_st)); // psi_{1,a}, normally psi^+_a this is definitely creation operator now.
	    if (i4==0){ // we can not compute ratio, but we still need to check that njm[i_st,:]==0.
	      double rest = sum(abs(njm(i_st,Range::all())));
	      if (rest>1e-12){
		Q_proportional=false;
		break;
	      }
	      //if (Print){
	      //cout<<"   i4="<<i4<<" rest="<<njm(i_st,Range::all())<<" st_bl="<<st_bl<<" Block["<<ib+1<<"]="<<Block[ib]<<" ia="<<ia<<endl;
	      //cout<<"   njm="<<njm<<endl;
	      //}
	    }else{
	      int ii4 = inBlock(abs(i4));
	      r = njm(i_st,ii4)/sign(i4);
	      double rest=0;
	      for (int k=0; k<njm.extent(1); ++k)
		if (k!=ii4) rest += abs(njm(i_st,k));
	      if (rest>1e-12) Q_proportional=false;
	      //if (Print){
	      //cout<<"   i4="<<i4<<" ii4="<<ii4<<" r="<<r<<" rest="<<njm(i_st,Range::all())<<" st_bl="<<st_bl<<" Block["<<ib+1<<"]="<<Block[ib]<<" ia="<<ia<<" Q_prop="<<Q_proportional<<endl;
	      //cout<<"   njm="<<njm<<endl;
	      //}
	      break;
	    }
	  }
	  if (Q_proportional){
	    for (int i_st=1; i_st<Block[ib].size(); ++i_st){// now we go over all atomic states in block ib
	      int i=Block[ib](i_st); // i is state |state0>==|block=ib,state=i>
	      int i4 = Fi(2*ia,i);     // psi_{ia}^+|state0==i>=|i4>
	      if (i4!=0){
		int ii4 = inBlock(abs(i4));
		double dr = fabs(r - njm(i_st,ii4)/sign(i4)); // is njm(i_st,:) equal to F^+_ia[i_str,:]
		//double rest = sum(abs(njm(i_st,SRange(0,i_st-1," n2")))) + sum(abs(njm(i_st,SRange(i_st+1,toEnd," n3"))));
		double rest=0;
		for (int k=0; k<njm.extent(1); ++k)
		  if (k!=ii4) rest += abs(njm(i_st,k));
		if (dr>1e-12 or rest>1e-12){
		  Q_proportional=false;
		  break;
		}
		//if (Print){
		// cout<<"i4="<<i4<<" ii4="<<ii4<<" r="<<r<<" dr="<<dr<<" rest="<<njm(i_st,Range::all())<<endl;
		//  cout<<"  njm="<<njm<<endl;
		//}
	      }else{
		if (sum(abs(njm(i_st,Range::all())))>1e-10){
		  Q_proportional=false;
		  break;
		}
	      }
	    }
	  }
	  //cout<<"Q_proportional="<<Q_proportional<<endl;
	  if (Q_proportional){
	    Njm_c(ib,ia)=1;
	    Njm_r[encode_block_bath(ib,ia)]=r;
	    if (Print) cout<<"Njm["<<setw(3)<<ib+1<<","<<ia<<"]="<<r<<"*F^+_{"<<ia<<"}"<<"   proportional"<<endl;
	  }else{
	    Njm_c(ib,ia)=2;
	    Njm[ib][ia].resize(njm.shape());
	    Njm[ib][ia]=njm;
	    if (Print){
	      cout<<"Njm["<<setw(3)<<ib+1<<","<<ia<<"]=";
	      for (int n1=0; n1<Njm[ib][ia].extent(0); ++n1){
		if (n1==0) cout<<"";
		else cout<<"           ";
		for (int n2=0; n2<Njm[ib][ia].extent(1); ++n2)
		  cout<<setw(3)<<Njm[ib][ia](n1,n2)<<" ";
		cout<<endl;
	      }

	      cout<<"Psi["<<setw(3)<<ib+1<<","<<ia<<"]=";
	      for (int n1=0; n1<Njm[ib][ia].extent(0); ++n1){
		if (n1==0) cout<<"";
		else cout<<"           ";
		int i4 = Fi(2*ia,Block[ib](n1));
		int ii4 = i4!=0 ? inBlock(abs(i4)) : -1;
		for (int n2=0; n2<Njm[ib][ia].extent(1); ++n2){
		  double val = (n2==ii4) ? sign(i4) : 0;
		  cout<<setw(3)<<val<<" ";
		}
		cout<<endl;
	      }
	    }
	  }
	}
      }
    }
    if (Print){
      for (int ib=0; ib<Nblock; ++ib){ // over all blocks of states
	for (int ia=0; ia<Nbath; ++ia){// over all baths, and expect Njm to be proportional to psi^_{ia}|state0>
	  NjmResult nx = Njm_ind(ib,2*ia);
	  cout<<"ib="<<ib<<", ia="<<ia<<" njm_c="<<nx.njm_c<<endl;
	  if (nx.njm_c==1){
	    cout<<"   r="<< nx.njm_r << "  "<<endl;
	  }else if (nx.njm_c==2){
	    const Array<double,2>& njm1 = *nx.njm;
	    int dim0 = nx.transposed ? njm1.extent(1) : njm1.extent(0);
	    int dim1 = nx.transposed ? njm1.extent(0) : njm1.extent(1);
	    cout<<"Njm["<<setw(3)<<ib+1<<","<<ia<<"]=";
	    for (int n1=0; n1<dim0; ++n1){
	      if (n1==0) cout<<"";
	      else cout<<"           ";
	      for (int n2=0; n2<dim1; ++n2){
		double nn = nx.transposed ? njm1(n2,n1) : njm1(n1,n2);
		cout<<setw(3)<< nn <<" ";
	      }
	      cout<<endl;
	    }
	  }
	}
      }
    }
  }

  
  TinyVector<int,2> decode_block_bath(int key) const { return TinyVector<int,2>(key/Nbath,key%Nbath); } // second index is bath index
  int encode(const TinyVector<int,4>& ms) const { return ms[0] + Nbath * (ms[1] + Nbath * (ms[2] + Nbath * ms[3]));} // all four are bath indices
  TinyVector<int,4> decode(int key) const { // all four are bath indices
    int m0 = key % Nbath; key = key/Nbath;
    int m1 = key % Nbath; key = key/Nbath;
    int m2 = key % Nbath; key = key/Nbath;
    int m3 = key % Nbath;
    return TinyVector<int,4>(m0, m1, m2, m3);
  }

  void PrintSymmetry(ostream& cout, bool pindx=true)
  {
    cout<<"# Sindex version 2"<<endl;
    cout<<"# Nblock "<<Nblock<<endl;
    cout<<"# Npseudo_unique "<<Npseudo_unique<<endl;
    cout<<"# first_anomalous "<<first_anomalous<<endl;
    cout<<"# block_sizes";
    for (int ib=0; ib<Nblock; ++ib) cout<<" "<<Block[ib].size();
    cout<<endl;
    for (int ib=0; ib<Nblock; ++ib){
      cout<<"S["<<ib+1<<"]="<<endl;
      for (int i=0; i<Block[ib].size(); ++i){
	cout<<"       ";
	for (int j=0; j<Block[ib].size(); ++j)
	  cout<<setw(4)<<Sgm_symmetry[ib](i,j)<<", ";
	cout<<endl;
      }
      cout<<endl;
    }
    if (pindx){
      for (int ii=0; ii<Npseudo_unique; ++ii){
	cout<<"Sgm_indx["<<ii+1<<"]: ";
	for (int j=0; j<Sgm_indx[ii].size(); ++j)
	  cout << "["<<Sgm_indx[ii][j][0]+1<<","<<Sgm_indx[ii][j][1]<<","<<Sgm_indx[ii][j][2]<<"], ";
	cout<<endl;
      }
    }
  }
  BathBlocks FindBlocks()
  {
    Array<int,1> InBlock(Npsi);
    for (int iorb=0; iorb<Npsi; ++iorb) InBlock(iorb)=iorb;
    Array<bool,1> _SC_(Npsi); _SC_=false; // is block Nambu block describing SC?
    for (int iorb=0; iorb<Nbath; ++iorb){
      int psi_dagger = Fi_indx(2*iorb);
      int psi = Fi_indx(2*iorb+1);
      int alpa = psi/2;
      int beta = psi_dagger/2;
      if (InBlock(alpa) < InBlock(beta)){
	InBlock(beta)=InBlock(alpa);
      }else if (InBlock(alpa) > InBlock(beta)){
	InBlock(alpa)=InBlock(beta);
      }
      bool sc = (psi_dagger+psi)%2==0;// we have two psi^+ or two psi, hence Nambu
      _SC_(alpa) = _SC_(alpa) or sc; // remember this indices
      _SC_(beta) = _SC_(beta) or sc; // are in Nambu block
    }
    BathBlocks bb;
    //vector<bool> SC;
    //vector<vector<int>> blocks;
    std::map<int, int> label_to_block_idx;
    for (int i = 0; i < Npsi; ++i) {
      int lbl = InBlock(i);  // label of this orbital
      auto it = label_to_block_idx.find(lbl);
      if (it == label_to_block_idx.end()) {
	// First time we see this label → create a new block
	int new_idx = static_cast<int>(bb.blocks.size());
	label_to_block_idx[lbl] = new_idx;
	bb.blocks.push_back(vector<int>());      // add empty block
	bb.blocks.back().push_back(i);           // put i into new block
	bb.is_block_sc.push_back(_SC_(i));
      } else {
	// Label already known → append to existing block
	int b = it->second;
	bb.blocks[b].push_back(i);
	bb.is_block_sc[b] = bb.is_block_sc[b] or _SC_(i);
      }
    }
    bb.orb_ind.resize(Nbath);
    for (int ib=0; ib<bb.blocks.size(); ++ib){
      int N = bb.blocks[ib].size();
      for (int j=0; j<N; ++j){
	for (int k=0; k<N; ++k){
	  int alpa = bb.blocks[ib][j];
	  int beta = bb.blocks[ib][k];
	  bool found=false;
	  int iorb=0;
	  for (iorb=0; iorb<Nbath; ++iorb)
	    if (Fi_indx(2*iorb+1)/2 == alpa and Fi_indx(2*iorb)/2 == beta){
	      found=true;
	      break;
	    }
	  if (found) bb.orb_ind[iorb]=TinyVector<int,3>(ib,j,k);
	}
      }
    }
    return bb;
  }
};

// Trim from both ends
void inline_trim(std::string& s) {
    // Left trim
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            [](unsigned char ch) { return !std::isspace(ch); }));
    // Right trim
    s.erase(std::find_if(s.rbegin(), s.rend(),
            [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}


void diag2op(Array<int,1>& op,
	     const Array<int,2>& diag, const Array<int,1>& dirspin)
{
  op.resize(2*diag.extent(0));
  for (int ibth=0; ibth<diag.extent(0); ++ibth){
    int i = diag(ibth,0);
    int j = diag(ibth,1);
    int ds = dirspin(ibth);
    op(j) = ds;
    op(i) = 1-(ds & 1) + 2*(ds >> 1);  // operator opposite to ds
  }
}

string op2str(int iop){
  if (iop%2==0)
    return "psi^+_{"+to_string(iop/2)+"}";
  else
    return "psi_{"+to_string(iop/2)+"}";
}

class Diag2Matrix{
public:
  map<int,vector<Array<int,2>>> allowed;
  map<int,Array<int,2>> diagonal;
  Array<int,1> current, starts;
  Array<int,1> op;
  // 
  Diag2Matrix(int N_max_order, const Cix& cix)
  {
    current.resize(cix.Nmaxsize);
    starts.resize(cix.Nmaxsize);
  }
  Diag2Matrix& operator=(const Diag2Matrix& rhs) {
    if (this == &rhs) return *this;
    // Helper: deep-copy a Blitz 2D array with correct shape
    auto deep_copy = [](const blitz::Array<int,2>& src) {
      blitz::Array<int,2> dst(src.shape());
      dst = src;  // elementwise copy after resize
      return dst; // NRVO/move elision
    };
    // Build temporaries first for strong exception safety
    std::map<int, blitz::Array<int,2>>              new_diagonal;
    std::map<int, std::vector<blitz::Array<int,2>>> new_allowed;
    //It turns out "diagonal = rhs.diagonal;" does not work, because Array<int,2> needs to be resized before copying.    
    // Deep-copy diagonal
    for (const auto& [key, arr] : rhs.diagonal) {
      new_diagonal.emplace(key, deep_copy(arr));
    }
    // allowed = rhs.allowed; // Copy std::map<int, vector<Array<int,2>>>
    // Deep-copy allowed (vector of Blitz arrays per key)
    for (const auto& [key, vec] : rhs.allowed) {
      std::vector<blitz::Array<int,2>> v;
      v.reserve(vec.size());
      for (const auto& arr : vec) {
        v.push_back(deep_copy(arr));
      }
      new_allowed.emplace(key, std::move(v));
    }
    // Commit all at once (no partial state if an allocation throws)
    diagonal.swap(new_diagonal);
    allowed.swap(new_allowed);
    
    //diagonal.clear();
    //for (const auto& [key, arr] : rhs.diagonal) {
    //  blitz::Array<int,2> tmp(arr.shape());
    //  tmp = arr; // Safe: tmp is resized first
    //  diagonal.emplace(key, std::move(tmp));
    //}
    
    // Copy Blitz arrays
    current.resize(rhs.current.shape());
    current = rhs.current;
    op.resize(rhs.op.shape());
    op = rhs.op;

    return *this;
  }
  auto Cmp_both(const Array<int,2>& diag, const Array<int,1>& dirspin, const Cix& cix)
    -> std::tuple< map<int,vector<Array<int,2>>>&, /*unordered_*/map<int,Array<int,2>>& >
  {
    diag2op(op,diag,dirspin);
    int Nd = 2*diag.extent(0);
    allowed.clear();
    diagonal.clear();
    for (int ib=0; ib<cix.Nblock; ++ib){
      int ib0 = ib;
      int ii=0;
      for (; ii<Nd; ++ii){ // going around the ring and see if this block survives.
        int ib1 = cix.Fbc(op(ii),ib0);
        if (ib1<0) break; // if we break, it did not survive
        ib0 = ib1;
      }
      if (ib0>=0 and ii==Nd){// we survived, hence this block contributes
        ib0 = ib;
	//cout<<"cix.Block["<<ib0<<"]"<<endl;
        int current_size = cix.Block[ib0].size();
        //current(SRange(0,current_size-1," n4")) = cix.Block[ib0];
        current(Range(0,current_size-1)) = cix.Block[ib0];
	vector<Array<int,2>> allowed_(Nd);
	int max_size=0;
        for (int i=0; i<Nd; ++i){
	  int ib1 = cix.Fbc(op(i),ib0);
	  cout << op2str(op(i)) <<" |"<<ib0<<"> -> |"<<ib1<<">"<<" Blck["<<ib0<<"]= "<<cix.Block[ib0]<<" Blck["<<ib1<<"]= "<<cix.Block[ib1]<<endl;
	  int il=0;
	  deque<TinyVector<int,2>> rr;
	  for (int j=0; j<cix.Block[ib0].size(); ++j){
	    int iq = cix.Block[ib0](j);
	    int ir = cix.Fi(op(i),iq);
	    if (ir!=0){
	      rr.push_back(TinyVector<int,2>(iq,ir));
	      cout<<"   "<<iq<<" "<<ir<<endl;
	      il+=1;
	    }
	  }
	  allowed_[i].resize(rr.size(),2);
	  if (rr.size()>max_size) max_size = rr.size();
	  for (int il=0; il<rr.size(); ++il){
	    allowed_[i](il,0) = rr[il][0];
	    allowed_[i](il,1) = rr[il][1];
	  }
	  //#############
	  il=0;
	  for (int j=0; j<current_size; ++j){
	    int iq=current(j);
	    int ir = cix.Fi(op(i),iq);
	    if (ir!=0)
	      current(il++)=abs(ir);
	  }
	  current_size=il;
	  ib0 = ib1;
        }
	if (max_size>0)
	  allowed.emplace(ib, std::move(allowed_));
        cout<<"....................."<<endl;
        if (current_size>0){
	  diagonal[ib].resize(Nd,current_size);
	  for (int i=0; i<Nd; ++i){
	    for (int ii=0; ii<current_size; ++ii) current(ii) = abs(cix.Fi(op(i),current(ii)));
	    //diagonal[ib](i,Range::all())=current(SRange(0,current_size-1, " n5"));
	    diagonal[ib](i,Range::all())=current(Range(0,current_size-1));
	  }
	  cout << "diagonal="<< diagonal[ib] << endl;
        }else{
	  cout<< "diagonal= []"<<endl;
        }
        cout << "---------------------" << endl;
	if (allowed.contains(ib)){
	  vector<Array<int,2>>& alw = allowed[ib];
	  for (int i=0; i<Nd; ++i){
	    cout<<"[";
	    int allowed_size=alw[i].extent(0);
	    int j=0;
	    for (; j<allowed_size-1; ++j) cout << alw[i](j,1) <<", ";
	    j=allowed_size-1;
	    cout << alw[i](j,1)<<"] x [";
	    j=0;
	    allowed_size = alw[(i+1)%Nd].extent(0);
	    for (; j<allowed_size-1; ++j) cout << alw[(i+1)%Nd](j,0) << ", ";
	    cout << alw[(i+1)%Nd](allowed_size-1,0) << "]" << endl;
	  }
	}
        cout << "---------------------" << endl;
      }
    }
    return std::tie(allowed, diagonal);
  }
  auto Cmp_matrix(const Array<int,2>& diag, const Array<int,1>& dirspin, const Cix& cix)
  //-> std::tuple< map<int,vector<Array<int,2>>>& >
    -> map<int,vector<Array<int,2>>>&
  { // allowed[iblock][2*Norder](npossible,2)
    // (iq,ir) PP-propagator just before and after operator op(i)
    // allowed[iblock][2*Norder](i,0) = iq // just before op(i)
    // allowed[iblock][2*Norder](i,1) = ir // just after op(i)
    diag2op(op,diag,dirspin);
    int Nd = 2*diag.extent(0);
    
    allowed.clear();
    for (int ib=0; ib<cix.Nblock; ++ib){
      int ib0 = ib;
      int ii=0;
      for (; ii<Nd; ++ii){// going around the ring and see if this block survives.
        int ib1 = cix.Fbc(op(ii),ib0);
        if (ib1<0) break;// if we break, it did not survice
        ib0 = ib1;
      }
      if (ib0>=0 and ii==Nd){// we survived, hence this block contributes
        ib0 = ib;
	vector<Array<int,2>> allowed_(Nd);
	int max_size=0;
        for (int i=0; i<Nd; ++i){
	  int ib1 = cix.Fbc(op(i),ib0);
	  int il=0;
	  deque<TinyVector<int,2>> rr;
	  for (int j=0; j<cix.Block[ib0].size(); ++j){
	    int iq = cix.Block[ib0](j); // PP-propagator just before operator op(i)
	    int ir = cix.Fi(op(i),iq);  // PP-propagator just after operator op(i)
	    if (ir!=0){
	      rr.push_back(TinyVector<int,2>(iq,ir));
	      //cout<<"   "<<iq<<" "<<ir<<endl;
	      il+=1;
	    }
	  }
	  // storing deque rr into allowed_, which is a vector
	  allowed_[i].resize(rr.size(),2);
	  if (rr.size()>max_size) max_size = rr.size();
	  for (int il=0; il<rr.size(); ++il){
	    allowed_[i](il,0) = rr[il][0];
	    allowed_[i](il,1) = rr[il][1];
	  }
	  ib0 = ib1;
        }
	if (max_size>0)
	  allowed.emplace(ib, std::move(allowed_));
      }
    }
    //return std::tie(allowed);
    return allowed;
  }
  /*
  auto Cmp_diagonal(const Array<int,2>& diag, const Array<int,1>& dirspin, const Cix& cix)
  //-> std::tuple<map<int,Array<int,2>>&>
    -> map<int,Array<int,2>>&
  {
    diag2op(op,diag,dirspin);
    diagonal.clear();
    int Nd = 2*diag.extent(0);
    for (int ib=0; ib<cix.Nblock; ++ib){
      int ib0 = ib;
      int ii=0;
      for (; ii<Nd; ++ii){ // going around the ring and see if this block survives.
        int ib1 = cix.Fbc(op(ii),ib0);
        if (ib1<0) break;  // if we break, it did not survice
        ib0 = ib1;
      } 
      if (ib0>=0 and ii==Nd){    // we survived, hence this block contributes
        int current_size = cix.Block[ib].size();
        current(SRange(0,current_size-1, " n6")) = cix.Block[ib]; // all states in this block to start with
	//cout<<"start with current="<<current<<endl;
        for (int i=0; i<Nd; ++i){
	  int il=0;
	  for (int j=0; j<current_size; ++j){
	    int iq=current(j);
	    int ir = cix.Fi(op(i),iq);
	    if (ir!=0)
	      current(il++)=abs(ir); // only ir!=0 contribute, hence il is smaller than current_size at ech iteration
	  }
	  //cout<<"  current["<<i<<"]="<<current(Range(0,il-1))<<endl;
	  current_size=il; // reducing it for only those states that survive.
        }
        if (current_size>0){ // should survive
	  diagonal[ib].resize(Nd,current_size); // going through once more and now remember it
	  for (int i=0; i<Nd; ++i){
	    for (int ii=0; ii<current_size; ++ii) current(ii) = cix.Fi(op(i),abs(current(ii)));
	    diagonal[ib](i,Range::all())=current(SRange(0,current_size-1, " n7"));
	  }
	  //cout<<"current_final="<<diagonal[ib]<<endl;
        }
      }
    }
    //return std::tie(diagonal);
    return diagonal;
  }
  */
  //auto Cmp_diagonal2(const Array<int,2>& diag, const Array<int,1>& dirspin, const Cix& cix)
  ////-> std::tuple<map<int,Array<int,2>>&>
  //  -> /*unordered_*/map<int,Array<int,2>>&
  //{
  //  diag2op(op,diag,dirspin);
  //  diagonal.clear();
  //  int Nd = 2*diag.extent(0);
  //  int total_size=0;
  //  for (int ib=0; ib<cix.Nblock; ++ib){
  //    int ib0 = ib;
  //    int _ib1_ =  cix.Fbc(op(0),ib0);
  //    int ii=0;
  //    for (; ii<Nd; ++ii){ // going around the ring and see if this block survives.
  //      int ib1 = cix.Fbc(op(ii),ib0);
  //      if (ib1<0) break;  // if we break, it did not survice
  //      ib0 = ib1;
  //    } 
  //    if (ib0>=0 and ii==Nd){    // we survived, hence this block contributes
  //      int current_size = cix.Block[ib].size();
  //      current(bl::Range(0,current_size-1)) = cix.Block[ib]; // all states in this block to start with
  //	//cout<<"start with current="<<current<<endl;
  //      for (int i=0; i<Nd; ++i){
  //	  int il=0;
  //	  for (int j=0; j<current_size; ++j){
  //	    int iq=current(j);
  //	    int ir = cix.Fi(op(i),iq);
  //	    if (ir!=0)
  //	      current(il++)=abs(ir); // only ir!=0 contribute, hence il is smaller than current_size at ech iteration
  //	  }
  //	  //cout<<"  current["<<i<<"]="<<current(Range(0,il-1))<<endl;
  //	  current_size=il; // reducing it for only those states that survive.
  //      }
  //	total_size += current_size;
  //      if (current_size>0){ // should survive
  //	  diagonal[ib].resize(current_size,Nd); // going through once more and now remember it
  //	  for (int i=0; i<Nd; ++i){
  //	    for (int j=0; j<current_size; ++j){
  //	      int _ab_ = cix.Fi(op(i),abs(current(j)));
  //	      int abs_ab = abs(_ab_);
  //	      int ab = cix.inBlock(abs_ab);
  //	      int blck = cix.iBlck(abs_ab)-1;
  //	      diagonal[ib](j,i)= cix.Sgm_symmetry[blck](ab,ab)*sign(_ab_);
  //	      current(j) = _ab_;
  //	    }
  //	  }
  //	  //cout<<"current_final="<<diagonal[ib]<<endl;
  //      }
  //    }
  //  }
  //  if (total_size==0){ // diagonal approximation vanishes, which is not good enough approximation, i.e., this part of the space would be completely arised.
  //    // We need to allow it by dinding at least one off-diagonal term (randomly).
  //    //deque<int> randomize;
  //    //for (int i=0; i<Nd-1; ++i)
  //    //if (op(i)==opn(i+1)) randomize.push_back(i+1);
  //    total_size=0;
  //    for (int ib=0; ib<cix.Nblock; ++ib){
  //	int ib0 = ib;
  //	int ii=0;
  //	for (; ii<Nd; ++ii){ // going around the ring and see if this block survives.
  //	  int ib1 = cix.Fbc(op(ii),ib0);
  //	  if (ib1<0) break;  // if we break, it did not survice
  //	  ib0 = ib1;
  //	}
  //	if (ib0>=0 and ii==Nd){    // we survived, hence this block contributes
  //	  int current_size = cix.Block[ib].size();
  //	  current(bl::Range(0,current_size-1)) = cix.Block[ib]; // all states in this block to start with
  //	  //cout<<"start with current="<<current<<endl;
  //	  cout<<"We know total_size=0, trying to get block ib0="<<ib0<<" to contribute. current="<<current(Range(0,current_size-1))<<endl;
  //	  int il=0, iq, ir;
  //	  int i=0;
  //	  for (; i<Nd; ++i){
  //	    il=0;
  //	    for (int j=0; j<current_size; ++j){
  //	      iq=current(j);
  //	      ir = cix.Fi(op(i),iq);
  //	      if (ir!=0){
  //		current(il)=abs(ir); // only ir!=0 contribute, hence il is smaller than current_size at ech iteration
  //		il+=1;
  //	      }
  //	    }
  //	    if (il==0){
  //	      cout<<"   iq="<<iq<<" and ir="<<ir<<" breaking out because il="<<il<<endl;
  //	      break; // no contribution, although we expected one, hence need off-diagonals
  //	    }
  //	    current_size = il;
  //	    cout<<"current["<<i<<"]="<<current(Range(0,current_size-1))<<endl;
  //	  }
  //	  Array<int,1> current_copy(current_size); current_copy = current(Range(0,current_size-1));
  //	  Array<int,2> current_ii(current_size,Nd);
  //	  int i_special=i;
  //	  Array<int,1> blck_special(current_size);
  //	  for (int j=0; j<current_size; ++j) blck_special(j) = cix.iBlck(current(j))-1;
  //	  
  //	  cout<<"il="<<il<<" i_special="<<i_special<<" current_copy="<<current_copy<<" blck_special="<<blck_special<<endl;
  //	  for (; i<Nd; ++i){
  //	    // repeat the same step, but select off-diagonal
  //	    il=0;
  //	    for (int j=0; j<current_size; ++j){
  //	      int iq=current(j);
  //	      int ir = cix.Fi(op(i),iq);
  //	      int a = cix.inBlock(iq);
  //	      int blck = cix.iBlck(iq)-1;
  //	      int Nblck = cix.Block[blck].size(); // size of the block
  //	      cout<<"iq="<<iq<<" ir="<<ir<<" a="<<a<<" blck="<<blck<<" Nblck="<<Nblck<<endl;
  //	      if (ir!=0){
  //		current_ii(j,i) = cix.Sgm_symmetry[blck](a,a);
  //		current(il)=abs(ir);
  //		cout<<" diagonal ["<<i<<","<<j<<"] is fine ii="<<cix.Sgm_symmetry[blck](a,a)<<endl;
  //		il+=1;
  //	      }else if (Nblck!=1){ // block of size=1 can not be changed.
  //		int b = a;  // before the propagator was diagonal, which made it vanish
  //		int ii;
  //		while (b==a or ii==0){
  //		  b = (drand48()*Nblck); // now b!=a, and we use G_{a,b} instead of G_{a,a}
  //		  ii = cix.Sgm_symmetry[blck](a,b);
  //		}
  //		current_ii(j,i) = ii;
  //		iq = cix.Block[blck](b); // now iq is corrected, because its off-diagonal
  //		int ir = cix.Fi(op(i),iq);
  //		cout<<" off-diagonal attempt ["<<i<<","<<j<<"] ii="<<ii<<" a="<<a<<" b="<<b<<" iq->"<<iq<<" ir="<<ir<<endl;
  //		if (ir!=0){
  //		  current(il)=abs(ir); // only ir!=0 contribute, hence il is smaller than current_size at ech iteration
  //		  il+=1;
  //		}
  //	      }
  //	    }
  //	    current_size=il; // reducing it for only those states that survive.
  //	    cout<<"  current["<<i<<"]="<<current(Range(0,il-1))<<endl;
  //	  }
  //	  cout<<"il="<<il<<" going backward from i_special="<<i_special<<" current_copy="<<current_copy<<endl;
  //	  for (int i=i_special-1; i>=0; --i){
  //	    int il=0;
  //	    for (int j=0; j<current_copy.size(); ++j){
  //	      int iq=current_copy(j);
  //	      int ds = op(i);
  //	      int op_opposite = 1-(ds & 1) + 2*(ds >> 1);  // operator opposite to op(i)
  //	      int ir = cix.Fi(op_opposite,iq);
  //	      cout<<"   i="<<i<<" j="<<j<<" iq="<<iq<<" op_opposite="<<op_opposite<<" ir="<<ir<<endl;
  //	      if (ir!=0){
  //		current_copy(il)=abs(ir); // only ir!=0 contribute, hence il is smaller than current_size at ech iteration
  //		int a = cix.inBlock(ir);
  //		int blck = cix.iBlck(iq)-1;
  //		ii = cix.Sgm_symmetry[blck](a,a);
  //		current_ii(j,i)=ii;
  //		il+=1;
  //		cout<<" ii="<<ii<<" a="<<a<<" blck="<<blck<<endl;
  //	      }
  //	    }
  //	    if (il!=current_copy.size()){
  //	      cout<<"ERROR"<<endl;
  //	    }
  //	  }
  //	  cout<<"current_ii="<<current_ii<<endl;
  //	  total_size += current_size;
  //	  if (current_size>0){ // should survive
  //	    diagonal[ib].resize(current_size,Nd); // going through once more and now remember it
  //	    diagonal[ib] = current_ii(Range(0,current_size-1),Range(0,Nd-1));
  //	    /*
  //	    for (int i=0; i<Nd; ++i){
  //	      for (int j=0; j<current_size; ++j){
  //		int _ab_ = cix.Fi(op(i),abs(current(j)));
  //		int abs_ab = abs(_ab_);
  //		int ab = cix.inBlock(abs_ab);
  //		int blck = cix.iBlck(abs_ab)-1;
  //		diagonal[ib](j,i)= cix.Sgm_symmetry[blck](ab,ab)*sign(_ab_);
  //		current(j) = _ab_;
  //	      }
  //	    }
  //	    */
  //	    //cout<<"current_final="<<diagonal[ib]<<endl;
  //	  }
  //	}
  //    }
  //  }
  //  //return std::tie(diagonal);
  //  return diagonal;
  //
  //}

  auto Cmp_diagonal3(const Array<int,2>& diag, const Array<int,1>& dirspin, const Cix& cix)
  //-> std::tuple<map<int,Array<int,2>>&>
    -> /*unordered_*/map<int,Array<int,2>>&
  {
    //cout<<"diag="<<diag<<" dirspin="<<dirspin<<endl;
    diag2op(op,diag,dirspin);
    diagonal.clear();
    int Nd = 2*diag.extent(0);
    int total_size=0;
    for (int ib=0; ib<cix.Nblock; ++ib){
      int ib0 = ib;
      int _ib1_ =  cix.Fbc(op(0),ib0);
      int ii=0;
      for (; ii<Nd; ++ii){ // going around the ring and see if this block survives.
        int ib1 = cix.Fbc(op(ii),ib0);
        if (ib1<0) break;  // if we break, it did not survice
        ib0 = ib1;
      } 
      if (ib0>=0 and ii==Nd){    // we survived, hence this block contributes
        int current_size = cix.Block[ib].size();
        current(Range(0,current_size-1)) = cix.Block[ib]; // all states in this block to start with
	starts(Range(0,current_size-1)) = cix.Block[ib];
	//cout<<"start with current="<<current<<" and op="<<op<<endl;
        for (int i=0; i<Nd; ++i){
	  int il=0;
	  for (int j=0; j<current_size; ++j){
	    int istart = starts(j);
	    int iq=current(j);
	    int ir = cix.Fi(op(i),iq);
	    if (ir!=0){
	      current(il)=abs(ir); // only ir!=0 contribute, hence il is smaller than current_size at ech iteration
	      starts(il)=istart;
	      il+=1;
	    }
	  }
	  //cout<<"  current["<<i<<"]="<<current(Range(0,il-1))<<" and starts="<<starts(Range(0,il-1))<<endl;
	  current_size=il; // reducing it for only those states that survive.
        }
	// now keep only closed loops: start == end
	int nloops = 0;
	for (int j = 0; j < current_size; ++j) {
	  if (starts(j) == current(j)) {
	    starts(nloops)  = starts(j);
	    current(nloops) = current(j);  // optional
	    ++nloops;
	  }
	}
	current_size = nloops;
	//if (current_size>0) cout<<"  current="<<current(Range(0,current_size-1))<<" and starts="<<starts(Range(0,current_size-1))<<endl;
	//else cout<<"  current=[] and starts=[]"<<endl;
	total_size += current_size;
        if (current_size>0){ // should survive
	  diagonal[ib].resize(current_size,Nd); // going through once more and now remember it
	  for (int i=0; i<Nd; ++i){
	    for (int j=0; j<current_size; ++j){
	      int _ab_ = cix.Fi(op(i),abs(current(j)));
	      if (_ab_==0){
		cout << "ERROR _ab_="<<_ab_<<" at current(j)="<<current(j)<<" and op("<<i<<")="<<op(i)<<endl;
		cout<<"current["<<j<<"]="<<current(j)<<" ib="<<ib<<",j="<<j<<",i="<<i<<" op(i)="<<op(i)<<endl;
	      }
	      int abs_ab = abs(_ab_);
	      int ab = cix.inBlock(abs_ab);
	      int blck = cix.iBlck(abs_ab)-1;
	      diagonal[ib](j,i)= cix.Sgm_symmetry[blck](ab,ab)*sign(_ab_);
	      current(j) = _ab_;
	      //cout<<"current["<<j<<"]="<<current(j)<<" diagonal[ib="<<ib<<"]("<<j<<","<<i<<")="<<diagonal[ib](j,i)<<" blck="<<blck<<" ab="<<ab<<" _ab_="<<_ab_<<endl;
	    }
	  }
	  //cout<<"current_final="<<diagonal[ib]<<endl;
        }
      }
    }
    if (total_size==0){
      // diagonal approximation vanishes, which is not good enough approximation.
      // Within such approximation this part of the space would be completely erased.
      // We need to allow it by finding at least one off-diagonal term (randomly).
      for (int ib=0; ib<cix.Nblock; ++ib){
	int ib0 = ib; // block index when vertex==0
	int ir=0;
	for (; ir<Nd; ++ir){// going around the ring and see if this block survives.
	  int ib1 = cix.Fbc(op(ir),ib0);
	  if (ib1<0) break;// if we break, it did not survice
	  ib0 = ib1;
	}
	if (ib0>=0 and ir==Nd){// we survived, hence this block should contribute
	  ib0 = ib; // back to block at vertex==0
	  vector<Array<int,2>> allowed_(Nd); // first we get the same quantity as needed for Cmp_matrix
	  int max_size=0;
	  for (int i=0; i<Nd; ++i){
	    int ib1 = cix.Fbc(op(i),ib0); // next block index
	    int il=0;
	    deque<TinyVector<int,2>> rr;  // rr={(iq,ir)} iq:before and ir:after applying |ir> = Fi(op(i))|iq>
	    for (int j=0; j<cix.Block[ib0].size(); ++j){
	      int iq = cix.Block[ib0](j);
	      int ir = cix.Fi(op(i),iq);
	      if (ir!=0){
		rr.push_back(TinyVector<int,2>(iq,ir));
		il+=1;
	      }
	    }
	    // storing deque rr into allowed_, which is a vector
	    allowed_[i].resize(rr.size(),2); // i goes around the ring. currently we have rr.size terms that survive
	    if (rr.size()>max_size) max_size = rr.size();
	    for (int il=0; il<rr.size(); ++il){
	      allowed_[i](il,0) = rr[il][0]; // |iq>
	      allowed_[i](il,1) = rr[il][1]; // |ir>=Fi(op(i))|iq>
	    }
	    ib0 = ib1;
	  }
	  // We just computed allowed_, which is output in Cmp_matrix. Now we will select one
	  // random path through all matrices.
	  //for (int i=0; i<Nd; ++i){
	  //  cout<<i<<" [";
	  //  for (int j=0; j<allowed_[i].extent(0); ++j)
	  //    cout<<allowed_[i](j,0)<<":"<<allowed_[i](j,1)<<",";
	  //  cout<<"] ";
	  //}
	  //cout<<endl;
	  diagonal[ib].resize(1,Nd); // only one possible path will be targeted. We don't need more than one.
	  int N0 = allowed_[0].extent(0);
	  int j0 = static_cast<int>(drand48()*N0); // first is one random index between 0...N0-1, which are possible at vertex=0
	  int j=j0;
	  for (int i=0; i<Nd; ++i){
	    int ip1 = (i+1)%Nd;
	    int _a_ = allowed_[i](j,1); // G_{a,b} and a=allowed_[](:,1)
	    int a = abs(_a_);
	    int ai = cix.inBlock(a);    // index inside the block
	    int blck = cix.iBlck(a)-1;  // and index of the block
	    int ii=0;
	    if (i<Nd-1){                // last vertex need to connect back to vertex 0, hence it is different
	      for (int l=0; l<allowed_[ip1].extent(0); ++l){ // checking if diagonal G_{a,a} exists. We want allowed[ip1](l,0)==|allowed[i](j,1)|
		if (allowed_[ip1](l,0)==a){
		  j = l; // yes, found diagonal propagator, and for the next vertex we set j
		  ii = cix.Sgm_symmetry[blck](ai,ai)*sign(_a_); // remember what is the index in G[ii]
		  //cout<<i<<") i-diagonal ("<<_a_<<","<<a<<") ii="<<ii<<" j="<<j<<" ai="<<ai<<" bi="<<ai<<" blck="<<blck<<endl;
		  break;
		}
	      }
	      int how_many_trials=0;
	      while (ii==0){ // if we found diagonal above ii!=0, and we skip.
		// if diagonal does not exist, we try to take a random off-diagonal term
		int l = static_cast<int>(allowed_[ip1].extent(0)*drand48()); // this is the random off-diagonal
		int b = allowed_[ip1](l,0); // b in G_{a,b}
		int bi = cix.inBlock(b);    // index inside the block for b
		ii = cix.Sgm_symmetry[blck](ai,bi)*sign(_a_); // index for G[ii]
		j = l; // j in the next vertex
		//cout<<i<<") o-diagonal ("<<_a_<<","<<b<<") ii="<<ii<<" j="<<j<<" ai="<<ai<<" bi="<<bi<<" blck="<<blck<<endl;
		how_many_trials+=1; // sometimes we just can't find finite contribution, and eventually we need to give up.
		if (how_many_trials>allowed_[ip1].extent(0)*2) break;
	      }
	      if (ii==0){// we tried maximum how_many_trials, and now we give up and erase this ib in diagonal.
		diagonal.erase(ib);
		break;
	      }
	    }else{ // last step to connect back to vertex==0. We have no freedom.
	      // We selected j0 for vertex==0, and now we need to tie back to it.
	      int b = allowed_[ip1](j0,0);
	      int bi = cix.inBlock(b);
	      ii = cix.Sgm_symmetry[blck](ai,bi)*sign(_a_);
	      j=j0;
	      //cout<<i<<") last diagonal ("<<_a_<<","<<b<<") ii="<<ii<<" j="<<j<<" ai="<<ai<<" bi="<<bi<<" blck="<<blck<<endl;
	    }
	    diagonal[ib](0,i)= ii;
	  }
	}
      }
    }
    //cout << "We end up with diagonal=";
    //for (auto const& [ib, alw] : diagonal)
    //  cout << "["<<ib<<"]="<< alw << ", ";
    //cout<<endl;
    return diagonal;
  }
};

struct DSU{
  // Disjoint-Set Union (DSU), also known as union-find.
  //
  std::vector<int> p; //parent pointer for each element (initially p[i]=i)
  std::vector<int> r; //rank (approximate depth) of each root (initially 0)
  DSU(int n): p(n), r(n,0) { for(int i=0;i<n;++i) p[i]=i; }
  // Returns an integer that labels the set containing x.
  // It uses path compression: while it climbs the parent pointers (p) to the root,
  // it rewires every visited node directly to that root, flattening the tree and making future queries faster.
  int find(int x){ return p[x]==x?x:p[x]=find(p[x]); }
  // Merges the two sets that contain a and b.
  // •	It first calls find on each element, giving the two roots a and b.
  // •	If the roots are the same the sets are already united.
  // •	Otherwise it links the shallower tree under the deeper one (union by rank).
  // •	r stores the rank (an upper bound on tree depth).
  // •	If the ranks are equal it arbitrarily chooses one root to be the new parent and increments its rank.
  void unite(int a,int b){
    a = find(a); b = find(b);
    if(a==b) return;
    if(r[a]<r[b]) p[a]=b;
    else if(r[b]<r[a]) p[b]=a;
    else { p[b]=a; ++r[a]; }
  }
};
/******************************************************************************
 * group_blocks
 * -----------
 * Given `indices`, a list of N sets of integers, return a vector of groups,
 * where each group is a vector<int> containing the indices of all blocks
 * that share at least one integer value.
 *
 * Two blocks belong to the same group if their sets overlap *directly* or
 * through a chain of overlaps (i.e. we take the transitive closure).
 *
 * Example
 *   indices[0] = {1, 2}
 *   indices[1] = {4}
 *   indices[2] = {2, 5}
 *   indices[3] = {8}
 *
 *   → groups = { {0,2}, {1}, {3} }
 *
 *  This routine is used to determine exponent a of mGexp class, which is a common 
 *  exponent shared in a block or a set of blocks that share common elements. Namely,
 *  if two blocks contain the same PP-function component, they must also share a common
 *  exponent, otherwise we can not consistently determine exponents for blocks and 
 *  also scalar list of functions, which stores only PP-unique functions.
 ******************************************************************************/
vector<vector<int>> group_blocks(const vector<set<int>>& indices){
  int N = indices.size();
  /*----------------------------------------------------------------------
    1.  Build a disjoint–set union (DSU) over block indices 0 … N-1.
        We will “union” two blocks whenever we discover that they share
        (at least) one common integer.
    ---------------------------------------------------------------------*/
  DSU dsu(N); // parent and rank initialised for all blocks
   /*----------------------------------------------------------------------
     2.  Map each integer → “first block that contained it”.
         As we iterate, when we encounter a value that has already been
         seen in block `j`, we union the *current* block with block `j`.
     ---------------------------------------------------------------------*/
  std::unordered_map<int,int> first; // value → block index
  first.reserve(1024);   // pre-allocate to avoid rehashes
  
  for(int ib=0; ib<N; ++ib){  // for every block
    for(int x: indices[ib]){  // for every value in it
      auto it = first.find(x);
      if(it == first.end()){
	// first time we see value `x`: remember its block
	first[x] = ib;
      } else {
	// `x` already appeared in block it->second → merge groups
	dsu.unite(ib, it->second);
      }
    }
  }
  /*----------------------------------------------------------------------
      3.  After the loop, `dsu` encodes all connected components.
          We now collect the blocks that share the same DSU root.
   ---------------------------------------------------------------------*/
  std::map<int, std::vector<int>> tmp;   // root → list of blocks
  //tmp.reserve(N);
  for(int ib=0; ib<N; ++ib){
    int root = dsu.find(ib);
    tmp[root].push_back(ib);
  }
  /*----------------------------------------------------------------------
    4.  Move the grouped vectors into the final result.
        (std::map keeps roots sorted; if you don’t care about order,
        an unordered_map would work, too.)
    ---------------------------------------------------------------------*/
  std::vector<std::vector<int>> groups;
  groups.reserve(tmp.size());
  for(auto &kv : tmp){
    groups.push_back(std::move(kv.second));  // take ownership
  }
  return groups;
}

struct Int128Hash {
    std::size_t operator()(const __int128& x) const noexcept {
        unsigned __int128 ux = (unsigned __int128)x;
        std::uint64_t low  = (std::uint64_t)ux;
        std::uint64_t high = (std::uint64_t)(ux >> 64);
        // Mix the two halves (this is a standard combination pattern)
        std::size_t h = std::hash<std::uint64_t>{}(low) 
                      ^ (std::hash<std::uint64_t>{}(high) + 0x9e3779b97f4a7c15ULL + (low<<6) + (low>>2));
        return h;
    }
};

// ======================================================================
//  assign_label
//
//  This function assigns a symmetry label (an integer) to the matrix
//  element Sgm[ib](i,j). The label is stored in Sgm_symmetry[ib](i,j).
//
//  The rules:
//
//   • If the value has appeared before, reuse its existing label.
//   • If the NEGATIVE of the value has appeared before, reuse the label
//     with opposite sign (because Sgm may have anti-symmetric patterns).
//   • Otherwise, this is a NEW unique value (up to sign), so assign a
//     new label (++inew) and record it in the hash table.
//
//  The “known_values” hash stores one representative location for each
//  discovered value, so we can copy the label from there.
//
//  Inputs:
//     Sgm           – vector of 2D __int128 arrays (self-energy matrices)
//     Sgm_symmetry  – vector of 2D int arrays (labels)
//     known_values  – maps (__int128 key) → (ib,i,j) representative
//     inew          – the next label to assign
//
// ======================================================================
void assign_label(std::vector< Array<int,2> >& Sgm_symmetry, std::unordered_map<__int128, std::tuple<int,int,int>, Int128Hash>& known_values, int& inew,
		  int ib, int i, int j, const std::vector< Array<__int128,2> >& Sgm)
{
    // Value of the self-energy element at this position
    __int128 key = Sgm[ib](i,j);

    // Zero entries get no label (remain 0)
    if (key == 0) {
        return;
    }

    // ------------------------------------------------------------------
    // Case 1: We have seen THIS exact value before
    // ------------------------------------------------------------------
    auto it = known_values.find(key);
    if (it != known_values.end()) {
        // Retrieve representative element and copy its label
        auto [ib1, i1, j1] = it->second;
        Sgm_symmetry[ib](i,j) = Sgm_symmetry[ib1](i1,j1);
        return;
    }

    // ------------------------------------------------------------------
    // Case 2: We have seen the NEGATIVE of this value before
    //         (common in antisymmetric or sign-related symmetries)
    // ------------------------------------------------------------------
    auto it2 = known_values.find(-key);
    if (it2 != known_values.end()) {
        auto [ib1, i1, j1] = it2->second;
        // Same label magnitude, but opposite sign
        Sgm_symmetry[ib](i,j) = -Sgm_symmetry[ib1](i1,j1);
        return;
    }

    // ------------------------------------------------------------------
    // Case 3: NEW unique value (up to sign).
    //         Assign a new label:  inew = 1,2,3,...
    // ------------------------------------------------------------------
    int new_label = ++inew;
    Sgm_symmetry[ib](i,j) = new_label;

    // Record this as the representative for this key
    known_values.emplace(key, std::make_tuple(ib, i, j));
}

void Determine_Symmetry(Cix& cix, Diag2Matrix& diag2matrix)
{
  bool DEBUG=false;
  vector<Array<__int128,2>> Sgm(cix.Nblock);
  vector<Array<long long,2>> Gm0(cix.Nblock);
  Array<double,1> work(3*cix.Nmaxsize);
  for (int ib=0; ib<cix.Nblock; ++ib){
    int N=cix.Block[ib].size();
    Sgm[ib].resize(N,N);
    Gm0[ib].resize(N,N);
    Array<double,2> H(cix.Ham[ib].copy()); // copy constructor. We don't want to modify cix.Ham. It will contain eigenvectors.
    Array<double,1> E(N), aEne(N); aEne = 1.0;
    if (cix.nmBlock[ib](0)==cix.nmBlock[ib](N-1)){ // we have one block, and can not block diagonalize Ham
      xsyev(H.extent(0), H.data(), H.extent(1), E.data(), work.data(), work.size(), "V");
      aEne = (sum(E)/N - cix.Emin)/(cix.Emax-cix.Emin);
    }else{
      // We will block-diagonalize instead of diagonalizing one large H
      int i=0;
      while (i < N) {
	int label = cix.nmBlock[ib](i);
	int j = i + 1;
	// find end of this contiguous block
	while (j < N && cix.nmBlock[ib](j) == label) ++j;
	int L = j - i;  // block size
	if (L == 1) {
	  // 1x1 block: eigenvalue is H(i,i), eigenvector is basis vector e_i
	  E(i) = H(i,i);
	  // Make column i of H be the unit vector e_i
	  H(Range::all(), i) = 0.0;
	  H(i, i) = 1.0;
	  aEne(i) = (E(i) - cix.Emin) / (cix.Emax - cix.Emin);
	} else {
	  // Extract LxL sub-block into a temporary
	  Array<double,2> Hb(L, L);
	  Hb(Range::all(),Range::all()) = H(Range(i,i+L-1),Range(i,i+L-1));
	  Array<double,1> Eb(L);
	  // Diagonalize this block
	  xsyev(L, Hb.data(), L, Eb.data(), work.data(), work.size(), "V");
	  // Copy eigenvectors back into the corresponding sub-block of H
	  H(Range(i,i+L-1),Range(i,i+L-1)) = Hb(Range::all(),Range::all());
	  // Copy eigenvalues into the global E array
	  E(Range(i,i+L-1)) = Eb(Range::all());
	  aEne(Range(i,i+L-1)) = (sum(Eb)/L - cix.Emin)/(cix.Emax-cix.Emin);
	}
	i = j; // move to next block
      }
    }
    for (int i=0; i<N; ++i){
      if (aEne(i) < 0.1) aEne(i)=0.1;
      if (aEne(i) > 2.0) aEne(i)=2.0;
    }
    aEne *= 10;
    
    // Note that diagonalization gives: Ham = H^+ * E * H
    // In our case
    // G = H^+ exp(-E*t) * H => G(i,j) = H^H(i,k) exp(-E(k)*t) H(k,j) = Ht^*(i,k) exp(-E(k)*t) Ht(j,k)
    double emin = bl::min(E);
    Array<double,1> expE(N);
    for (int k=0; k<N; ++k) expE(k) = exp(-(E(k)-emin));
    Array<double,2> gt(N,N); gt=0;
    for (int i=0; i<N; ++i){
      for (int j=0; j<N; ++j){
	for (int k=0; k<N; ++k) gt(i,j) += H(k,i)*(expE(k)+E(k)/7.732)*H(k,j); // if Ham is complex, we would need to conjugate first Ht(i,k)
	// This was added later to break the symmetry between blocks that correspond to very different valency, and can be accidentally degenerate.
	// The off diagonal terms tend to be very hard to break in high-symmetry cases. While diagonal terms are already shifted by E->E-mu*N, the off-diagonal
	// terms in highly symmetric situations remain degenerate for different valencies, which leads to wrong symmetry index. So, we here multiply
	// g0[i,j]*=aEne[i]*aEne[j] where aEne[i] = (<E>^{block}-Emin)/(Emax-Emin)*10. This just multiplies each block of constant N with a different number, to break the symmetry.
	gt(i,j) *= aEne(i)*aEne(j); 
      }
    }
    Array<double,1> tmp(N); tmp = (expE+E/7.732);
    //cout<<"g00["<<ib+1<<"]="<< tmp << endl;
    //cout<<"g0["<<ib+1<<"]="<<gt<<endl;
    //cout<<"H["<<ib+1<<"]="<<H<<endl;
    //cout<<"aEne="<<ib+1<<"]="<<aEne<<endl;
    
    for (int i=0; i < N; ++i){
      for (int j=0; j < N; ++j)
	Gm0[ib](i,j) = static_cast<long long>(gt(i,j)*1000);
    }
    //cout<<"ib="<<ib+1<<" expE="<<expE<<" Ek="<<E<<" H="<<H<<" gt="<<gt<<endl;
    Sgm[ib] = Gm0[ib];
  }
  if (DEBUG){
    for (int ib=0; ib<cix.Nblock; ++ib){
      cout<<"G0["<<ib+1<<"]="<<endl;
      for (int i=0; i < cix.Block[ib].size(); ++i){
	cout<<"    ";
	for (int j=0; j < cix.Block[ib].size(); ++j)
	  cout<<setw(10)<<Gm0[ib](i,j)<<" ";
	cout<<endl;
      }
    }
    cout<<"bath_sigind="<<cix.bath_sigind << endl;
    cout<<"fake_sigind=[";
    for (int iorb=0; iorb<cix.Nbath; ++iorb){
      int iib0 = cix.bath_sigind0(iorb); 
      double ib_sg = (iib0 < 0 ? -1 : 1);
      int iib   = std::abs(iib0);
      int fake_hybrid = (iib+1)*ib_sg; // hybridization is set to integer
      cout<<fake_hybrid;
      if (iorb<cix.Nbath-1) cout<<",";
    }
    cout<<"]"<<endl;
    cout<<"Fi="<<cix.Fi<<endl;
  }
  int Nitt=5;
  for (int itt=0; itt<Nitt; ++itt){// NCA type diagrams
    const int Nd=2;
    Array<int,2> diag(1,2); diag=0,1;
    Array<int,1> dirspin(1);
    for (int iorb=0; iorb<cix.Nbath; ++iorb){
      dirspin = 2*iorb; //
      int iib0 = cix.bath_sigind0(iorb); 
      double ib_sg = (iib0 < 0 ? -1 : 1);
      int iib   = std::abs(iib0);
      int fake_hybrid = (iib+1)*ib_sg; // hybridization is set to integer
      auto& allowed2 = diag2matrix.Cmp_matrix(diag,dirspin,cix);
      if (DEBUG) cout<<"iorb="<<iorb<<" dirspin="<<dirspin<<" fake_hybrid="<<fake_hybrid<<endl;
      for (const auto& [ib, alw] : allowed2){
	if (alw[0].extent(0)>0 and alw[1].extent(0)>0){
	  if (DEBUG){
	    for (int i=0; i<Nd; ++i){
	      int allowed_size = alw[i].extent(0);
	      if (allowed_size>0){
		cout<<"[";
		int j=0;
		for (; j<allowed_size-1; ++j) cout << alw[i](j,1) <<", ";
		cout << alw[i](allowed_size-1,1)<<"] x [";
		j=0;
		allowed_size = alw[(i+1)%Nd].extent(0);
		for (; j<allowed_size-1; ++j) cout << alw[(i+1)%Nd](j,0) << ", ";
		cout << alw[(i+1)%Nd](allowed_size-1,0) << "]" << endl;
	      }
	    }
	  }
	  TinyVector<int,2> ibs(cix.Fbc(diag2matrix.op(0),ib), ib);
	  TinyVector<int,Nd> js;
	  TinyVector<int,2> l_0, l_1;
	  for (js[0]=0; js[0]<alw[0].extent(0); ++js[0]){
	    for (js[1]=0; js[1]<alw[1].extent(0); ++js[1]){
	      //               js[0]
	      //            delta  alpha   where Fi[2*orb+1,delta]=alpha or psi_{orb}|delta>=|alpha>
	      //              -->----
	      //            /    |    \   block ibs[0]
	      //     G[jj] | D[orb,t]  |  G[ii]=G(alpha,beta)=G[ib](_alpha_,_beta_)
	      //           |     v     | 
	      //block ibs[1]\    |    /
	      //              ---<---
	      //            gamma  beta   where Fi[2*orb,beta]=gamma or psi^+_{orb}|beta>=|gamma>
	      //                js[1]
	      //                 t
	      int alpha = alw[0](js[0],1);
	      int beta  = alw[1](js[1],0);
	      int gamma = alw[1](js[1],1);
	      int delta = alw[0](js[0],0);
	      int _alpha_ = cix.inBlock(abs(alpha));
	      int _beta_ = cix.inBlock(abs(beta));
	      int _delta_ = cix.inBlock(abs(delta));
	      int _gamma_ = cix.inBlock(abs(gamma));
	      __int128 Phi = 3*fake_hybrid*Gm0[ibs[0]](_alpha_,_beta_)*Gm0[ibs[1]](_gamma_,_delta_)*sign(alpha)*sign(gamma);
	      {
		//int negative_time_sign = cix.bath_fermion(iorb) ? 1 : -1; // bosonic propagator does not cancel with permutation of psi operators
		bool SC0 = cix.nmBlock[ibs[0]](_alpha_)!=cix.nmBlock[ibs[0]](_beta_);
		bool SC1 = cix.nmBlock[ibs[1]](_gamma_)!=cix.nmBlock[ibs[1]](_delta_);
		double mm0 = sign(alpha)*sign(gamma);
		double mm1=mm0, mm2=mm0;
		// This is very subtle issue. We have the following expression for
		//   Sigma_{delta,gamma}(t2-t1) = mm0 * G_{alpha,beta}(t2-t1) * Delta(t2-t1)
		//   Sigma_{beta,alpha}(beta-t2+t1) = mm0 * G_{gamma,delta}(beta-t2+t1) * Delta(t2-t1)
		// which can also be written as
		//   Sigma_{beta,alpha}(t2-t1) = mm0 * G_{gamma,delta}(t2-t1) * Delta(beta-t2+t1)
		// but when Delta is superconducting, we have a symmetry Delta(beta-t)=-Delta(t), hence for SC bath, we have extra minus for Sigma_{beta,alpha}, i.e.,
		//   Sigma_{beta,alpha}(t2-t1) = -mm0 * G_{gamma,delta}(t2-t1) * Delta(t2-t1)
		// Now we can approximate Delta(t2-t1) with a number amd G_p(t2-t1) with another number, and the symmetry should give correct result
		if (cix.bath_SC(iorb)) mm1 *= -1;
		long long S_beta_alpha = 3*static_cast<__int128>(fake_hybrid*mm1)*Gm0[ibs[1]](_gamma_,_delta_); 
		Sgm[ibs[0]](_beta_,_alpha_) += S_beta_alpha;
		long long S_delta_gamma = 3*static_cast<__int128>(fake_hybrid*mm2)*Gm0[ibs[0]](_alpha_,_beta_);// here we have (+1) because G_{a,b}(t)*Delta(t)
		Sgm[ibs[1]](_delta_,_gamma_) += S_delta_gamma;
		if (DEBUG and (SC0 or SC1) and (S_beta_alpha!=0 or S_delta_gamma!=0)){
		  //if (mm<0 and (S_beta_alpha!=0 or S_delta_gamma!=0)){
		  cout << "  iorb="<<iorb<<" j0="<<js[0]<<" j1="<<js[1]<<" Phi="<<static_cast<long long>(Phi)<<" fh="<<3*fake_hybrid<<" sgn="<<sign(alpha)*sign(gamma)<<" mm_{ab}="<<mm1<<" m_{gd}="<<mm2;
		  cout << " G("<<alpha<<","<<beta<<")=G["<<ibs[0]+1<<"("<<_alpha_<<","<<_beta_<<")]="<<Gm0[ibs[0]](_alpha_,_beta_)<<"* G("<<gamma<<","<<delta<<")=G["<<ibs[1]+1<<"("<<_gamma_<<","<<_delta_<<")]="<<Gm0[ibs[1]](_gamma_,_delta_)<<endl;
		  cout << "  Sgm[("<<beta<<","<<alpha<<")="<<ibs[0]+1<<"]("<<_beta_<<","<<_alpha_<<")+="<< S_beta_alpha;
		  cout << ";   Sgm[("<<delta<<","<<gamma<<")="<<ibs[1]+1<<"]("<<_delta_<<","<<_gamma_<<")+="<<S_delta_gamma<<endl;
		}
	      }
	    }
	  }
	}
      }
    }
    // G0 <- G0 + Sigma_NCA
    for (int ib=0; ib<cix.Nblock; ++ib) Gm0[ib](Range::all(),Range::all()) = Sgm[ib](Range::all(),Range::all());

    if (DEBUG){
      cout<<"Itteration "<<itt<<" finished"<<endl;
      for (int ib=0; ib<cix.Nblock; ++ib){
	cout<<"S["<<ib+1<<"]="<<endl;
	for (int i=0; i<cix.Block[ib].size(); ++i){
	  cout<<"       ";
	  for (int j=0; j<cix.Block[ib].size(); ++j)
	    cout<<setw(4)<<Gm0[ib](i,j)<<", ";
	  cout<<endl;
	}
	cout<<endl;
      }
    }
    
  }

  bool second_order=false;
  if (second_order){
    Array<int,2> diag(2,2);
    diag = 0,2,1,3;
    Array<int,1> dirspin(2);
    vector<Array<long long,2>> Gs(4), Gr(4);
    const int Nd=4;
    for (int orb1=0; orb1<cix.Nbath; ++orb1){
      for (int orb2=0; orb2<cix.Nbath; ++orb2){
	dirspin = 2*orb1, 2*orb2;
	int iib0_1 = cix.bath_sigind0(orb1);
	int iib0_2 = cix.bath_sigind0(orb2);
	double ib_sg1 = (iib0_1 < 0 ? -1 : 1);
	double ib_sg2 = (iib0_2 < 0 ? -1 : 1);
	int iib1   = std::abs(iib0_1);
	int iib2   = std::abs(iib0_2);
	int fake_hybrid1 = (iib1+1)*ib_sg1; // hybridization is set to integer
	int fake_hybrid2 = (iib2+1)*ib_sg2; // hybridization is set to integer
	auto& allowed2  = diag2matrix.Cmp_matrix(diag,dirspin,cix);
	for (const auto& [ib, alw] : allowed2){
	  if (alw[0].extent(0)>0 and alw[1].extent(0)>0 and alw[2].extent(0)>0 and alw[3].extent(0)>0){
	    TinyVector<int,Nd> ibs;
	    ibs[0] = cix.Fbc(diag2matrix.op(0),ib);
	    ibs[1] = cix.Fbc(diag2matrix.op(1),ibs[0]);
	    ibs[2] = cix.Fbc(diag2matrix.op(2),ibs[1]);
	    ibs[3] = ib;
	    vector<vector<int>> l_0(Nd), l_1(Nd);
	    for (int i=0; i<Nd; ++i){
	      int ic = i, id=(i+1)%Nd; // we concentrate on propagator between operator i and operator i+1
	      int ibc = ibs[i];        // this propagator should be in block ibc
	      int dim1 = alw[ic].extent(0), dim2 = alw[id].extent(0); // dimensions of this block;
	      
	      l_0[i].resize(dim1);
	      for (int j0=0; j0<dim1; ++j0)
		l_0[i][j0] = cix.inBlock(abs(alw[ic](j0,1)));
	      
	      l_1[i].resize(dim2);
	      for (int j1=0; j1<dim2; ++j1)
		l_1[i][j1] = cix.inBlock(abs(alw[id](j1,0)));
	      
	      Gs[i].resize(dim1, dim2);
	      Gr[i].resize(dim1, dim2);
	      for (int j0=0; j0<dim1; ++j0){
		for (int j1=0; j1<dim2; ++j1){
		  Gr[i](j0,j1) = Gm0[ibc](l_0[i][j0],l_1[i][j1]);
		  if (sign(alw[id](j1,0))<0) cout<<"ERROR alw["<<id<<"]("<<j1<<","<<0<<")<0 = "<<alw[id](j1,0)<<endl;
		}
		Gs[i](j0,Range::all()) = Gr[i](j0,Range::all())*sign(alw[ic](j0,1));
	      }
	    }
	    TinyVector<int,Nd> js;
	    for (js[0]=0; js[0]<alw[0].extent(0); ++js[0]){
	      for (js[1]=0; js[1]<alw[1].extent(0); ++js[1]){
		for (js[2]=0; js[2]<alw[2].extent(0); ++js[2]){
		  for (js[3]=0; js[3]<alw[3].extent(0); ++js[3]){
		    //__int128 Phi = static_cast<__int128>(3*fake_hybrid1*fake_hybrid2)*static_cast<__int128>(Gs[0](js[0],js[1])*Gs[1](js[1],js[2]))*Gs[2](js[2],js[3])*Gs[3](js[3],js[0]);
		    __int128 Phi_common = static_cast<__int128>(3*fake_hybrid1*fake_hybrid2);
		    for (int i=0; i<Nd; ++i){
		      int j0 = js[i], j1 = js[(i+1)%Nd];
		      int bosonic_minus1 = 1;//cix.bath_fermion(orb1)==false and (0<=i and i<2) ? -1 : 1;
		      int bosonic_minus2 = 1;//cix.bath_fermion(orb2)==false and (1<=i and i<3) ? -1 : 1;
		      int Fi0 = sign(alw[i](js[i],1));
		      __int128 prod=1;
		      for (int j=0; j<Nd; ++j)
			if (j!=i) prod *= Gs[j](js[j],js[(j+1)%Nd]);
		      
		      Sgm[ibs[i]](l_1[i][j1],l_0[i][j0]) += Phi_common*Fi0*bosonic_minus1*bosonic_minus2*prod;
		      //Sgm[ibs[i]](l_1[i][j1],l_0[i][j0]) += Phi/Gr[i](j0,j1);
		    }
		    //}
		    //const __int128 K = 3LL * fake_hybrid1 * fake_hybrid2;
		    //const long long a0 = Gs[0](js[0], js[1]);  // includes sign
		    //const long long a1 = Gs[1](js[1], js[2]);
		    //const long long a2 = Gs[2](js[2], js[3]);
		    //const long long a3 = Gs[3](js[3], js[0]);
		    //// Phi == K*a0*a1*a2*a3;
		    //if (a1 && a2 && a3)
		    //  Sgm[ibs[0]](l_1[0][js[1]], l_0[0][js[0]]) += K * a1 * a2 * a3; // Phi/a0
		    //if (a0 && a2 && a3)
		    //  Sgm[ibs[1]](l_1[1][js[2]], l_0[1][js[1]]) += K * a0 * a2 * a3; // Phi/a1
		    //if (a0 && a1 && a3)
		    //  Sgm[ibs[2]](l_1[2][js[3]], l_0[2][js[2]]) += K * a0 * a1 * a3; // Phi/a2
		    //if (a0 && a1 && a2)
		    //  Sgm[ibs[3]](l_1[3][js[0]], l_0[3][js[3]]) += K * a0 * a1 * a2; // Phi/a3
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  // ======================================================================
  //   SYMMETRY LABELING FOR SELF–ENERGY MATRICES
  //
  //   Goal:
  //     • Assign consistent integer labels to matrix entries Sgm[ib](i,j)
  //     • Entries inside the normal state Hamiltonian block (intra-block) should
  //       get SMALL labels.
  //     • Entries connecting different blocks (exchange of Cooper pair, inter-block)
  //       should get LARGER labels, assigned only after all intra-block labels.
  //
  //   Why two phases?
  //     • We want block-internal structure to dominate: these represent
  //       the “intra-Hilbert-space” operators.
  //     • Entries between blocks represent anomalous propagators, and we want them
  //       to appear at the end.
  // ======================================================================
  // Resize and initialize the symmetry arrays
  cix.Sgm_symmetry.resize(cix.Nblock);
  for (int ib=0; ib<cix.Nblock; ++ib){
    int N = cix.Block[ib].size();
    cix.Sgm_symmetry[ib].resize(N,N);
    cix.Sgm_symmetry[ib]=0;
  }
  unordered_map<__int128, tuple<int,int,int>, Int128Hash> known_values;
  int inew=0;
  /*
  for (int ib=0; ib<cix.Nblock; ++ib){
    const Array<int,1>& indices = cix.Block[ib];
    for (int i=0; i<indices.size(); ++i){
      int j=i;
      // first fill the diagonal values
      __int128 key = Sgm[ib](i,j);
      if (key!=0){
	auto it = known_values.find(key);
	if (it != known_values.end()){
	  std::tuple<int,int,int>& tup = it->second;
	  auto [ib1, i1, j1] = tup;
	  cix.Sgm_symmetry[ib](i,j) = cix.Sgm_symmetry[ib1](i1,j1);
	}else{
	  auto it2 = known_values.find(-key);
	  if (it2 != known_values.end()){
	    std::tuple<int,int,int>& tup = it2->second;
	    auto [ib1, i1, j1] = tup;
	    cix.Sgm_symmetry[ib](i,j) = -cix.Sgm_symmetry[ib1](i1,j1);
	  }else{
	    cix.Sgm_symmetry[ib](i,j) = ++inew;
	    known_values[key] = std::tuple(ib,i,j);
	  }
	}
      }
      // than fill the off-diagonal, This is just for convenience.
      for (int j=0; j<indices.size(); ++j){
	if (j==i) continue;
	__int128 key = Sgm[ib](i,j);
	if (key!=0){
	  auto it = known_values.find(key);
	  if (it != known_values.end()){
	    std::tuple<int,int,int>& tup = it->second;
	    auto [ib1, i1, j1] = tup;
	    cix.Sgm_symmetry[ib](i,j) = cix.Sgm_symmetry[ib1](i1,j1);
	  }else{
	    auto it2 = known_values.find(-key);
	    if (it2 != known_values.end()){
	      std::tuple<int,int,int>& tup = it2->second;
	      auto [ib1, i1, j1] = tup;
	      cix.Sgm_symmetry[ib](i,j) = -cix.Sgm_symmetry[ib1](i1,j1);
	    }else{
	      cix.Sgm_symmetry[ib](i,j) = ++inew;
	      known_values[key] = std::tuple(ib,i,j);
	    }
	  }
	}
      }
    }
  }
  */
  // ----------------------------------------------------------------------
  //  PHASE 1: INTRA-BLOCK LABELING (assign *small* labels first)
  //  We process all (i,j) pairs where nmBlock(i) == nmBlock(j).
  //  These correspond to elements of the normal state blocks.
  //  Any new value discovered in this phase will get a low label.
  // ----------------------------------------------------------------------
  for (int ib = 0; ib < cix.Nblock; ++ib) {
    const Array<int,1>& nm = cix.nmBlock[ib];     // block label per index
    int N = cix.Block[ib].size();
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
	if (nm(i) == nm(j)) {  // both indices inside the same block
	  assign_label(cix.Sgm_symmetry, known_values, inew, ib, i, j, Sgm);
	}
      }
    }
  }
  cix.first_anomalous = inew+1;
  // ----------------------------------------------------------------------
  //  PHASE 2: INTER-BLOCK LABELING (assign *large* labels later)
  //  Now we process pairs (i,j) where nmBlock(i) != nmBlock(j).
  //  These represent anomalous propagators due to superconductivity.
  //  Any new values found here will get labels greater than all the
  //  intra-block ones, making the structure much clearer.
  // ----------------------------------------------------------------------
  for (int ib = 0; ib < cix.Nblock; ++ib) {
    const Array<int,1>& nm = cix.nmBlock[ib];
    int N = cix.Block[ib].size();
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
	if (nm(i) != nm(j)) {  // different block labels → inter-block element
	  assign_label(cix.Sgm_symmetry, known_values, inew, ib, i, j, Sgm);
	}
      }
    }
  }
  
  cix.Npseudo_unique = inew;
  cix.Sgm_indx.resize(cix.Npseudo_unique);
  for (int ib=0; ib<cix.Nblock; ++ib){
    for (int i=0; i<cix.Block[ib].size(); ++i){
      for (int j=0; j<cix.Block[ib].size(); ++j){
	int _ii_ = cix.Sgm_symmetry[ib](i,j);
	if (_ii_!=0) cix.Sgm_indx[abs(_ii_)-1].push_back(TinyVector<int,3>(ib*sign(_ii_),i,j));
      }
    }
  }
  
  cix.deg_pp.resize(cix.Npseudo_unique); cix.deg_pp=0; // degeneracy of each pseudoparticle
  for (int ib=0; ib<cix.Nblock; ++ib){
    int Ni = cix.Block[ib].size();
    for (int _i1_=0; _i1_<Ni; ++_i1_){
      for (int _i2_=0; _i2_<Ni; ++_i2_){
	int _ii_ = cix.Sgm_symmetry[ib](_i1_,_i2_);
	int _jj_ = cix.Sgm_symmetry[ib](_i2_,_i1_);
	if (_ii_ != _jj_){
	  cout<<"WARNING PPGreens function is not symmetric and current code does not work correctly. Need complex G' and G''!"<<endl;
	}
	if (_ii_==0) continue;
	cix.deg_pp(abs(_ii_)-1) += 1;
      }
    }
  }
  //cout<<"Nbath_unique="<<cix.Nbath_unique<<endl;
  cix.deg_bth.resize(cix.Nbath_unique); cix.deg_bth=0; // degeneracy of each bath
  for (int iorb=0; iorb<cix.Nbath; ++iorb) cix.deg_bth(cix.bath_sigind(iorb))+= 1;
  
  {
    vector<set<int>> indices_in_block(cix.Nblock);
    for (int ib=0; ib<cix.Nblock; ++ib){
      int Ni = cix.Block[ib].size();
      for (int a=0; a<Ni; ++a)
	for (int b=0; b<Ni; ++b)
	  if (cix.Sgm_symmetry[ib](a,b)!=0)
	    indices_in_block[ib].insert( abs(cix.Sgm_symmetry[ib](a,b)) );
      //cout<<"indices_in_block["<<ib<<"]="<<indices_in_block[ib]<<endl;
    }
    cix.unique_blocks = group_blocks(indices_in_block);
    //for (int i=0; i<unique_blocks.size(); ++i){
    //  cout<<"unique_blocks["<<i<<"]="<<unique_blocks[i]<<endl;
    //}
  }
  if (DEBUG){
    cout<<"Npseudo_unique="<<cix.Npseudo_unique<<endl;
    for (int ii=0; ii<cix.Npseudo_unique; ++ii){
      cout<<"Sgm_indx["<<ii<<"]: ";
      for (int j=0; j<cix.Sgm_indx[ii].size(); ++j)
	cout << cix.Sgm_indx[ii][j] << ", ";
      cout<<endl;
    }
    cout<<endl;
    for (int ib=0; ib<cix.Nblock; ++ib){
      cout<<"S["<<ib+1<<"]="<<endl;
      for (int i=0; i<cix.Block[ib].size(); ++i){
	cout<<"    ";
	for (int j=0; j<cix.Block[ib].size(); ++j)
	  cout<<setw(4)<<cix.Sgm_symmetry[ib](i,j)<<", ";
	cout<<"    ";
	for (int j=0; j<cix.Block[ib].size(); ++j)
	  cout<<setw(16)<< static_cast<long long>(Sgm[ib](i,j))<<", ";
	cout<<endl;
      }
      cout<<endl;
    }
  }
}

void FinalizeSindex(Cix& cix, int header_first_anomalous=-1)
{
  cix.Sgm_indx.clear();
  cix.Sgm_indx.clear();
  cix.unique_blocks.clear();

  int iimax=0;
  for (int i=0; i<cix.Sgm_symmetry.size(); ++i){
    int imax = max(cix.Sgm_symmetry[i]);
    if (imax>iimax) iimax=imax;
    //cout<<"S["<<i+1<<"]="<<cix.Sgm_symmetry[i]<<endl;
  }
  cix.Npseudo_unique = iimax;
  cix.Sgm_indx.resize(cix.Npseudo_unique);
  for (int ib=0; ib<cix.Nblock; ++ib){
    for (int i=0; i<cix.Block[ib].size(); ++i){
      for (int j=0; j<cix.Block[ib].size(); ++j){
	int _ii_ = cix.Sgm_symmetry[ib](i,j);
	if (_ii_!=0) cix.Sgm_indx[abs(_ii_)-1].push_back(TinyVector<int,3>(ib*sign(_ii_),i,j));
      }
    }
  }
  
  cix.deg_pp.resize(cix.Npseudo_unique); cix.deg_pp=0; // degeneracy of each pseudoparticle
  for (int ib=0; ib<cix.Nblock; ++ib){
    int Ni = cix.Block[ib].size();
    for (int _i1_=0; _i1_<Ni; ++_i1_){
      for (int _i2_=0; _i2_<Ni; ++_i2_){
	int _ii_ = cix.Sgm_symmetry[ib](_i1_,_i2_);
	int _jj_ = cix.Sgm_symmetry[ib](_i2_,_i1_);
	if (_ii_ != _jj_){
	  cout<<"WARNING PPGreens function is not symmetric and current code does not work correctly. Need complex G' and G''!"<<endl;
	}
	if (_ii_==0) continue;
	cix.deg_pp(abs(_ii_)-1) += 1;
      }
    }
  }
  cix.deg_bth.resize(cix.Nbath_unique); cix.deg_bth=0; // degeneracy of each bath
  for (int iorb=0; iorb<cix.Nbath; ++iorb) cix.deg_bth(cix.bath_sigind(iorb))+= 1;
  
  {
    vector<set<int>> indices_in_block(cix.Nblock);
    for (int ib=0; ib<cix.Nblock; ++ib){
      int Ni = cix.Block[ib].size();
      for (int a=0; a<Ni; ++a)
	for (int b=0; b<Ni; ++b)
	  if (cix.Sgm_symmetry[ib](a,b)!=0)
	    indices_in_block[ib].insert( abs(cix.Sgm_symmetry[ib](a,b)) );
      //cout<<"indices_in_block["<<ib<<"]="<<indices_in_block[ib]<<endl;
    }
    cix.unique_blocks = group_blocks(indices_in_block);
  }

  if (header_first_anomalous > 0){
    cix.first_anomalous = header_first_anomalous;
  }else{
    int first_anomalous = cix.Npseudo_unique+1;
    for (int ib=0; ib<cix.Nblock; ++ib){
      int Ni = cix.Block[ib].size();
      for (int a=0; a<Ni; ++a){
        for (int b=0; b<Ni; ++b){
          if (cix.nmBlock[ib](a) != cix.nmBlock[ib](b)){
            int lbl = abs(cix.Sgm_symmetry[ib](a,b));
            if (lbl != 0 && lbl < first_anomalous) first_anomalous = lbl;
          }
        }
      }
    }
    cix.first_anomalous = first_anomalous;
  }
}

void ReadSindex(string fname, Cix& cix)
{
  std::ifstream in(fname);
  if (!in)
    throw std::runtime_error("Cannot open file: " + fname);

  std::string line;
  std::vector<std::vector<int>> rows;
  int header_nblock = -1;
  int header_npseudo = -1;
  int header_first_anomalous = -1;
  std::vector<int> header_block_sizes;

  cix.Sgm_symmetry.clear();

  auto flush_block = [&]() {
    if (rows.empty()) return;
    const int nrow = rows.size();
    const int ncol = rows[0].size();
    for (int i = 1; i < nrow; ++i){
      if (rows[i].size() != ncol)
        throw std::runtime_error("Inconsistent row width while reading " + fname);
    }
    Array<int,2> A(nrow, ncol);
    for (int i = 0; i < nrow; ++i)
      for (int j = 0; j < ncol; ++j)
	A(i,j) = rows[i][j];
    cix.Sgm_symmetry.push_back(A);
    rows.clear();
  };

  while (std::getline(in, line)) {
    if (!line.empty() && line[0] == '#'){
      std::stringstream hs(line.substr(1));
      std::string key;
      hs >> key;
      if (key == "Nblock"){
        hs >> header_nblock;
      }else if (key == "Npseudo_unique"){
        hs >> header_npseudo;
      }else if (key == "first_anomalous"){
        hs >> header_first_anomalous;
      }else if (key == "block_sizes"){
        int bs;
        while (hs >> bs) header_block_sizes.push_back(bs);
      }
      continue;
    }
    if (line.find("S[") != std::string::npos) {
      flush_block();
      continue;
    }
    if (line.find_first_not_of(" \t\r") == std::string::npos)
      continue;
    for (char& ch : line) {
      if (ch == ',') ch = ' ';
    }
    std::stringstream ss(line);
    std::vector<int> row;
    int val;
    while (ss >> val)
      row.push_back(val);
    if (!row.empty()){
      rows.push_back(row);
    }
  }
  flush_block();

  if (cix.Sgm_symmetry.size() != cix.Nblock){
    throw std::runtime_error("Sindex block count mismatch in " + fname + ": expected " + to_string(cix.Nblock) + ", got " + to_string(cix.Sgm_symmetry.size()));
  }
  if (header_nblock >= 0 && header_nblock != cix.Nblock){
    throw std::runtime_error("Sindex header Nblock mismatch in " + fname);
  }
  if (!header_block_sizes.empty()){
    if (header_block_sizes.size() != cix.Nblock)
      throw std::runtime_error("Sindex header block_sizes mismatch in " + fname);
    for (int ib=0; ib<cix.Nblock; ++ib){
      if (header_block_sizes[ib] != cix.Block[ib].size())
        throw std::runtime_error("Sindex header block size mismatch in " + fname + " for block " + to_string(ib+1));
    }
  }
  for (int ib=0; ib<cix.Nblock; ++ib){
    if (cix.Sgm_symmetry[ib].extent(0) != cix.Block[ib].size() || cix.Sgm_symmetry[ib].extent(1) != cix.Block[ib].size()){
      throw std::runtime_error("Sindex matrix size mismatch in " + fname + " for block " + to_string(ib+1));
    }
  }

  FinalizeSindex(cix, header_first_anomalous);
  if (header_npseudo >= 0 && header_npseudo != cix.Npseudo_unique){
    throw std::runtime_error("Sindex header Npseudo_unique mismatch in " + fname);
  }
}


void Bcast_Sindex(Cix& cix, MPIcls* mpi)
{    
  int first_anomalous = cix.first_anomalous;
  int nblock = cix.Nblock;
  mpi->bcast_ints(&first_anomalous, 1);
  int total_size = 0;
  for (int ib=0; ib<nblock; ++ib) total_size += cix.Block[ib].size()*cix.Block[ib].size();
  Array<int,1> flat(total_size);
  if (mpi->rank==mpi->master){
    int pos = 0;
    for (int ib=0; ib<nblock; ++ib){
      int Ni = cix.Block[ib].size();
      for (int i=0; i<Ni; ++i){
	for (int j=0; j<Ni; ++j){
	  flat(pos++) = cix.Sgm_symmetry[ib](i,j);
	}
      }
    }
  }
  mpi->bcast_ints(flat.data(), total_size);
  if (mpi->rank!=mpi->master){
    cix.Sgm_symmetry.clear();
    int pos = 0;
    for (int ib=0; ib<nblock; ++ib){
      int Ni = cix.Block[ib].size();
      Array<int,2> A(Ni,Ni);
      for (int i=0; i<Ni; ++i){
	for (int j=0; j<Ni; ++j){
	  A(i,j) = flat(pos++);
	}
      }
      cix.Sgm_symmetry.push_back(A);
    }
    FinalizeSindex(cix, first_anomalous);
  }
}

/*
struct DisjointSet {
    std::vector<int> parent;
    std::vector<int> rank;

    DisjointSet(int n) : parent(n), rank(n, 0) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    int find(int x) {
        if (parent[x] != x)
            parent[x] = find(parent[x]);  // path compression
        return parent[x];
    }

    void unite(int x, int y) {
        int rx = find(x);
        int ry = find(y);
        if (rx == ry) return;
        if (rank[rx] < rank[ry]) {
            parent[rx] = ry;
        } else if (rank[rx] > rank[ry]) {
            parent[ry] = rx;
        } else {
            parent[ry] = rx;
            rank[rx]++;
        }
    }
};
std::vector<std::vector<int>> find_blocks(
    int Norb,
    const std::vector<std::pair<int,int>>& nonzeros)
{
    DisjointSet dsu(Norb);

    // Unite connected orbitals
    for (auto [i, j] : nonzeros) {
        // safety: ignore out-of-range
        if (i < 0 || j < 0 || i >= Norb || j >= Norb) continue;

        // Diagonal element (i,i) just ensures the orbital exists; unite(i,j)
        // is still harmless even if i == j.
        dsu.unite(i, j);
    }

    // Group orbitals by their representative root
    std::unordered_map<int, std::vector<int>> tmp;
    tmp.reserve(Norb);
    for (int i = 0; i < Norb; ++i) {
        int r = dsu.find(i);
        tmp[r].push_back(i);
    }

    // Flatten into a vector of blocks
    std::vector<std::vector<int>> blocks;
    blocks.reserve(tmp.size());
    for (auto &kv : tmp) {
        blocks.push_back(std::move(kv.second));
    }

    return blocks;
}
*/
/*
int main()
{
  Cix cix;
  cix.Read_cix("bctqmc.cix", 8.0);
  //cix.Read_cix("single_band.cix", 1.0);

  Diag2Matrix diag2matrix(6, cix);

  diag2matrix.Determine_Symmetry(cix);
  
  diag2matrix.PrintSymmetry(cix);

  return 0;
*/
  /*  
  auto [allowed1, diagonal2] = diag2matrix.Cmp_both(diag,dirspin,cix);
  auto [diagonal1]           = diag2matrix.Cmp_diagonal(diag,dirspin,cix);

  int Nd = diag.extent(0)*2;
  
  
  
  for (const auto& [ib, diagon] : diagonal2){
    cout<<" diagonal["<<setw(2)<<setw(2)<<ib<<"]=";
    for (int i=0; i<Nd; ++i){
      for (int il=0; il<diagon.extent(1); ++il){
	cout<<setw(2)<<diagon(i,il)<<",";
      }
      cout<<"; ";
    }
    cout<<endl;
  }



  Array<int,1> iN_blc(cix.Nblock);
  for (int ib=0; ib<cix.Nblock; ++ib){
    iN_blc(ib) = cix.iN(cix.Block[ib](0)-1);
    //cout<<"ib="<<ib<<" block="<<cix.Block[ib]<<" "<<cix.Block[ib](0)-1<<" N="<<cix.iN(cix.Block[ib](0)-1)<<endl;
  }
  //cout<<"iN_blc="<<iN_blc<<endl;

  vector<Array<long long,2>> Sgm(cix.Nblock), Gm0(cix.Nblock);
  for (int ib=0; ib<cix.Nblock; ++ib){
    Sgm[ib].resize(cix.Block[ib].size(),cix.Block[ib].size());
    Gm0[ib].resize(cix.Block[ib].size(),cix.Block[ib].size());
    for (int i=0; i < cix.Block[ib].size(); ++i){
      for (int j=0; j < cix.Block[ib].size(); ++j)
	Gm0[ib](i,j) = static_cast<long long>(cix.Ham[ib](i,j)*1000);
      //Gm0[ib](i,i) += iN_blc(ib)*7;
    }
    for (int i=0; i<Gm0[ib].extent(0); i++)
      if (fabs(Gm0[ib](i,i))==0) Gm0[ib](i,i) += 100;
    Sgm[ib] = Gm0[ib];
  }

  
  for (int ii=0; ii<1; ++ii){
    const int Nd=2;
    for (int iorb=0; iorb<cix.Nbath; ++iorb){
      dirspin = 2*iorb;
      int fake_hybrid = (cix.bath_sigind(iorb)+1);
      auto [allowed2] = diag2matrix.Cmp_matrix(diag,dirspin,cix);
      for (const auto& [ib, alw] : diag2matrix.allowed2){
	if (alw[0].extent(0)>0 and alw[1].extent(0)>0){
	  TinyVector<int,2> ibs(cix.Fbc(diag2matrix.op(0),ib), ib);
	  //cout<<"ibs="<<ibs<<endl;
	  //cout<<"Block["<<ibs[0]<<"]="<<cix.Block[ibs[0]]<<" Block["<<ibs[1]<<"]="<<cix.Block[ibs[1]]<<endl;
	  //cout<<"alw[0]="<<alw[0]<<" alw[1]="<<alw[1]<<endl;
	  // just printing!
	  TinyVector<int,Nd> js;
	  TinyVector<int,2> l_0, l_1;
	  for (js[0]=0; js[0]<alw[0].extent(0); ++js[0]){
	    for (js[1]=0; js[1]<alw[1].extent(0); ++js[1]){
	      l_0[0] = cix.inBlock(abs(alw[0](js[0],1)));
	      l_1[1] = cix.inBlock(abs(alw[0](js[0],0)));
	      l_1[0] = cix.inBlock(abs(alw[1](js[1],0)));
	      l_0[1] = cix.inBlock(abs(alw[1](js[1],1)));
	      //cout<<"l_0="<<l_0<<" l_1="<<l_1<<endl;
	      long long Phi = 3*fake_hybrid*Gm0[ibs[0]](l_0[0],l_1[0])*Gm0[ibs[1]](l_0[1],l_1[1])*sign(alw[0](js[0],1))*sign(alw[1](js[1],1));
	      //cout<<"Phi="<<Phi<<" G1="<<Gm0[ibs[1]](l_0[1],l_1[1])<<" G0="<<Gm0[ibs[0]](l_0[0],l_1[0])<<endl;
	      if (Phi!=0){
		Sgm[ibs[0]](l_1[0],l_0[0]) += Phi/Gm0[ibs[0]](l_0[0],l_1[0]); //3*fake_hybrid*Gm0[ibs[1]](l_0[1],l_1[1])*sign(k0_0)*sign(k1_0);
		Sgm[ibs[1]](l_1[1],l_0[1]) += Phi/Gm0[ibs[1]](l_0[1],l_1[1]); //3*fake_hybrid*Gm0[ibs[0]](l_0[0],l_1[0])*sign(k0_0)*sign(k1_0);
	      }
	    }
	  }
	}
      }
    }
    for (int ib=0; ib<cix.Nblock; ++ib) Gm0[ib](Range::all(),Range::all()) = Sgm[ib](Range::all(),Range::all());
  }

  bool second_order=true;
  if (second_order){
    Array<int,2> diag(2,2);
    diag = 0,2,1,3;
    Array<int,1> dirspin(2);
    cout<<"Second order"<<endl;
    vector<Array<long long,2>> Gs(4), Gr(4);
    const int Nd=4;
    for (int orb1=0; orb1<cix.Nbath; ++orb1){
      for (int orb2=0; orb2<cix.Nbath; ++orb2){
	dirspin = 2*orb1, 2*orb2;
	int fake_hybrid1 = (cix.bath_sigind(orb1)+1);
	int fake_hybrid2 = (cix.bath_sigind(orb2)+1);
	//cout<<"orb1="<<orb1<<" orb2="<<orb2<<" dirspin="<<dirspin<<" hyb1="<<fake_hybrid1<<" hyb2="<<fake_hybrid2<<endl;
	auto [allowed2]  = diag2matrix.Cmp_matrix(diag,dirspin,cix);
	for (const auto& [ib, alw] : diag2matrix.allowed2){
	  if (alw[0].extent(0)>0 and alw[1].extent(0)>0 and alw[2].extent(0)>0 and alw[3].extent(0)>0){
	    TinyVector<int,Nd> ibs;
	    ibs[0] = cix.Fbc(diag2matrix.op(0),ib);
	    ibs[1] = cix.Fbc(diag2matrix.op(1),ibs[0]);
	    ibs[2] = cix.Fbc(diag2matrix.op(2),ibs[1]);
	    ibs[3] = ib;
	    //cout<<"ibs="<<ibs<<endl;
	    //cout<<"alw[0]="<<alw[0]<<" alw[1]="<<alw[1]<<" alw[2]="<<alw[2]<<" alw[3]="<<alw[3]<<endl;
	    // just printing!
	    vector<vector<int>> l_0(Nd), l_1(Nd);
	    for (int i=0; i<Nd; ++i){
	      int ic = i, id=(i+1)%Nd; // we concentrate on propagator between operator i and operator i+1
	      int ibc = ibs[i];        // this propagator should be in block ibc
	      int dim1 = alw[ic].extent(0), dim2 = alw[id].extent(0); // dimensions of this block;
	      
	      l_0[i].resize(dim1);
	      for (int j0=0; j0<dim1; ++j0)
		l_0[i][j0] = cix.inBlock(abs(alw[ic](j0,1)));
	      
	      l_1[i].resize(dim2);
	      for (int j1=0; j1<dim2; ++j1)
		l_1[i][j1] = cix.inBlock(abs(alw[id](j1,0)));
	      
	      Gs[i].resize(dim1, dim2);
	      Gr[i].resize(dim1, dim2);
	      for (int j0=0; j0<dim1; ++j0){
		for (int j1=0; j1<dim2; ++j1){
		  Gr[i](j0,j1) = Gm0[ibc](l_0[i][j0],l_1[i][j1]);
		  if (sign(alw[id](j1,0))<0) cout<<"ERROR alw["<<id<<"]("<<j1<<","<<0<<")<0 = "<<alw[id](j1,0)<<endl;
		}
		Gs[i](j0,Range::all()) = Gr[i](j0,Range::all())*sign(alw[ic](j0,1));
	      }
	    }
	    TinyVector<int,Nd> js;
	    for (js[0]=0; js[0]<alw[0].extent(0); ++js[0]){
	      for (js[1]=0; js[1]<alw[1].extent(0); ++js[1]){
		for (js[2]=0; js[2]<alw[2].extent(0); ++js[2]){
		  for (js[3]=0; js[3]<alw[3].extent(0); ++js[3]){
		    long long Phi = 3*fake_hybrid1*fake_hybrid2*Gs[0](js[0],js[1])*Gs[1](js[1],js[2])*Gs[2](js[2],js[3])*Gs[3](js[3],js[0]);
		    if (Phi!=0){
		      for (int i=0; i<Nd; ++i){
			int j0 = js[i], j1 = js[(i+1)%Nd];
			Sgm[ibs[i]](l_1[i][j1],l_0[i][j0]) += Phi/Gr[i](j0,j1);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  vector<Array<double,2>> Sgm_symmetry(cix.Nblock);
  for (int ib=0; ib<cix.Nblock; ++ib){
    Sgm_symmetry[ib].resize(cix.Block[ib].size(),cix.Block[ib].size());
    Sgm_symmetry[ib]=0;
  }
  unordered_map<long long,tuple<int,int,int>> known_values;
  int inew=0;
  for (int ib=0; ib<cix.Nblock; ++ib){
    Array<int,1>& indices = cix.Block[ib];
    for (int i=0; i<indices.size(); ++i){
      for (int j=0; j<indices.size(); ++j){
	long long key = Sgm[ib](i,j);
	if (key!=0){
	  auto it = known_values.find(key);
	  if (it != known_values.end()){
	    std::tuple<int,int,int>& tup = it->second;
	    auto [ib1, i1, j1] = tup;
	    //if (iN_blc(ib1)==iN_blc(ib)){
	    Sgm_symmetry[ib](i,j) = Sgm_symmetry[ib1](i1,j1);
	    //}else{
	    //Sgm_symmetry[ib](i,j) = ++inew;
	    //known_values[key] = std::tuple(ib,i,j);
	    //}
	  }else{
	    Sgm_symmetry[ib](i,j) = ++inew;
	    known_values[key] = std::tuple(ib,i,j);
	  }
	}
      }
    }
  }
  
  
  for (int ib=0; ib<cix.Nblock; ++ib){
    cout<<"S["<<ib+1<<"]="<<endl;
    for (int i=0; i<cix.Block[ib].size(); ++i){
      cout<<"       ";
      for (int j=0; j<cix.Block[ib].size(); ++j)
	cout<<setw(4)<<Sgm_symmetry[ib](i,j)<<", ";
      
      cout<<endl;
    }
    cout<<endl;
  }
  return 0;


  
  int all_nn = 0;
  for (const auto& [ib, alw] : diag2matrix.allowed2){
    cout<<"iblock="<<ib<<endl;
    int ann=1, bnn=0;
    for (int i=0; i<Nd; ++i){
      int allowed_size = alw[i].extent(0);
      if (allowed_size>0){
	cout<<"[";
	int j=0;
	for (; j<allowed_size-1; ++j) cout << alw[i](j,1) <<", ";
	cout << alw[i](allowed_size-1,1)<<"] x [";
	j=0;
	allowed_size = alw[(i+1)%Nd].extent(0);
	for (; j<allowed_size-1; ++j) cout << alw[(i+1)%Nd](j,0) << ", ";
	cout << alw[(i+1)%Nd](allowed_size-1,0) << "]" << endl;

	int nn = alw[i].extent(0)*alw[(i+1)%Nd].extent(0);
	ann *= alw[i].extent(0);
	bnn += nn;
      }
    }
    all_nn += ann + bnn;
    cout<<"ann="<<ann<<" all_nn="<<all_nn<<endl;
  }
  return 0;
  */
  //}

#endif //_CIX_
