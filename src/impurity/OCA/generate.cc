// @Copyright 2007 Kristjan Haule
// 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <algorithm>

using namespace std;
typedef vector<int>::size_type vint;

// Crucial function that recursively generates the local base for any given number of baths ans its degeneracy
void generateBase(int b, vector<int>& state, vector<vector<int> >& base, int baths, const vector<int>& Ns) 
{ 
  for(int i = 0; i<=Ns[b]; i++){
    state[b] = i;
    //    state[baths-b-1] = i;
    if (b<baths-1) 
      generateBase(b+1, state,base,baths,Ns);
    else{      
      base.push_back(state);
    }
  }
}

// This comparisson is needed to sort the base by the total number of electrons in the
// local state. Inside the subspace of the same number of particles the states are sorted
// such that firts bit is first filled. This is not necessary but just more convenient
// when there is only one particle per bath.
class cmp{
  const vector<int>& Ns;
public:
  cmp(const vector<int>& Ns_) : Ns(Ns_){};
  bool operator()(const vector<int>& a, const vector<int>& b){
    int tota = 0; int totb = 0;
    for (vint i=0; i<a.size(); i++) tota += a[i];
    for (vint i=0; i<b.size(); i++) totb += b[i];
    if (tota!=totb) return tota<totb;
    int t=1;
    int linda=0, lindb=0;
    for (vint i=0; i<a.size(); i++) {
      linda += a[i]*t;
      lindb += b[i]*t;
      t*=Ns[i];
    }
    return linda<lindb;
  }
};

// Just old good !
int Binomial(int n, int m)
{
  double Mf = 1;
  for (int i=2; i<=m; i++) Mf *= i;
  double r = 1;
  for (int i=n; i>=n-m+1; i--) r*=i;
  return static_cast<int>(r/Mf);
}

// For parsing an entity in the parsInputChoice.
// Given a string like 10-12 return par of numbers
// 10 and 12.
pair<int,int> parsRange(const string& str){
  string::size_type pos;
  if ((pos=str.find('-'))<string::npos){
    int start = atoi(str.substr(0,pos).c_str());
    int end   = atoi(str.substr(pos+1).c_str());
    return make_pair(start,end);
  } else{
    int start = atoi(str.c_str());
    return make_pair(start,start);
  }
}

// Used to parse the input choice when the user cuts the base.
// Gets a string from the input (like 1-10,12,13) and returns
// a list of int numbers (in this 1,2,....10,12,13).
void parsInputChoice(string& str, list<int>& small, int base_size)
{
  list<pair<int,int> > keep;
  string::size_type pos;
  while((pos=str.find(','))<string::npos){
    keep.push_back(parsRange(str.substr(0,pos)));
    str.erase(0,pos+1);
  }
  
  keep.push_back(parsRange(str));
  for (list<pair<int,int> >::iterator i=keep.begin(); i!=keep.end(); i++)
    for (int j=(i->first); j<=(i->second); j++)
      if (j<base_size && j>=0) small.push_back(j);   
}

// Given the base and bath degeneracy calculates index to the base, total number of particles
// in the state, degeneracy of each state and NCA diagrams.
void SetUp(const vector<vector<int> >& base, const vector<int>& Ns, map<vector<int>,int>& index,
	   vector<int>& Mtot, vector<int>& degeneracy){
  for(vint i=0; i<base.size(); i++){
    int nt=0; int dg=1;
    for (vint j=0; j<base[i].size(); j++){
      nt += base[i][j];
      dg *= Binomial(Ns[j],base[i][j]);
    }
    Mtot[i] = nt;
    degeneracy[i] = dg;
    index[base[i]]=i+1;
  }

}

void CalcNCADiag(int baths, const vector<vector<int> >& base, const vector<int>& Ns, map<vector<int>,int>& index,
		 const vector<int>& degeneracy, vector<vector<int> >& ncab, vector<vector<int> >& ncaf)
{
  for (vint i=0; i<base.size(); i++){
    for (int b=0; b<baths; b++){
      vector<int> st = base[i];
      st[b]++;
      int indb = index[st]-1;
      ncab[i].push_back(indb);
      st[b]--; st[b]--;
      int indf = index[st]-1;
      ncaf[i].push_back(indf);
    }
  }
}

// void CalcUNCADiag(int baths, const vector<vector<int> >& base, const vector<int>& Ns, map<vector<int>,int>& index,
// 		  const vector<int>& degeneracy, vector<vector<vector<int> > >& sobb, vector<vector<vector<int> > >& sobf,
// 		  vector<vector<vector<int> > >& sofb, vector<vector<vector<int> > >& soff)
// {
//   for (vint i=0; i<base.size(); i++){
//     for (int b1=0; b1<baths; b1++){
//       for (int b2=0; b2<baths; b2++){
// 	vector<int> st = base[i];
// 	vector<int> diagbb(3), diagff(3), diagbf(3), diagfb(3);
// 	st[b1]++;
// 	diagbb[0] = index[st]-1;
// 	st[b2]++;
// 	diagbb[1] = index[st]-1;
// 	st[b1]--;
// 	diagbb[2] = index[st]-1;
// 	st[b2]--;
// 	sobb[i].push_back(diagbb);
	
// 	st[b1]--;
// 	diagff[0] = index[st]-1;
// 	st[b2]--;
// 	diagff[1] = index[st]-1;
// 	st[b1]++;
// 	diagff[2] = index[st]-1;
// 	st[b2]++;
// 	soff[i].push_back(diagff);

// 	st[b1]++;
// 	diagbf[0] = index[st]-1;
// 	st[b2]--;
// 	diagbf[1] = index[st]-1;
// 	st[b1]--;
// 	diagbf[2] = index[st]-1;
// 	st[b2]++;
// 	sobf[i].push_back(diagbf);

// 	st[b1]--;
// 	diagfb[0] = index[st]-1;
// 	st[b2]++;
// 	diagfb[1] = index[st]-1;
// 	st[b1]++;
// 	diagfb[2] = index[st]-1;
// 	st[b2]--;
// 	sofb[i].push_back(diagfb);
//       }
//     }
//   }
// }

void CalcUNCADiag_(int baths, const vector<vector<int> >& base, const vector<int>& Ns, map<vector<int>,int>& index,
		   const vector<int>& degeneracy,
		   vector<vector<vector<int> > >& sobb, vector<vector<vector<int> > >& sobf,
		   vector<vector<vector<int> > >& sofb, vector<vector<vector<int> > >& soff,
		   vector<vector<vector<int> > >& sogl,
		   vector<vector<int> >& prbb, vector<vector<int> >& prbf,
		   vector<vector<int> >& prfb, vector<vector<int> >& prff, vector<vector<int> >& prgl)
{
  for (vint i=0; i<base.size(); i++){
    for (int b1=0; b1<baths; b1++){
      for (int b2=0; b2<baths; b2++){
	vector<int> st = base[i];
	vector<int> diagbb(3), diagff(3), diagbf(3), diagfb(3);
	// Only half of diagrams need to be calculated
	// The second half is equivalent due to symmetry
	// Only factor of 2 is added to the prefactor
	int prefact = (b1!=b2) ? 2*(Ns[b1]-st[b1])*(Ns[b2]-st[b2]) : (Ns[b1]-st[b1])*(Ns[b2]-st[b2]-1);
	if (b1<=b2){
	  st[b1]++;
	  diagbb[0] = index[st]-1;
	  st[b2]++;
	  diagbb[1] = index[st]-1;
	  st[b1]--;
	  diagbb[2] = index[st]-1;
	  st[b2]--;
	} else {
	  diagbb[0] = -1;
	  diagbb[1] = -1;
	  diagbb[2] = -1;
	}
	sobb[i].push_back(diagbb);
	prbb[i].push_back(prefact);
      }
    }
  }
  
  for (vint i=0; i<base.size(); i++){
    for (int b1=0; b1<baths; b1++){
      for (int b2=0; b2<baths; b2++){
	vector<int> st = base[i];
	vector<int> diagbb(3), diagff(3), diagbf(3), diagfb(3);
	// Only half of diagrams need to be calculated
	// The second half is equivalent due to symmetry
	// Only factor of 2 is added to the prefactor
	int prefact = (b1!=b2) ? 2*st[b1]*st[b2] : st[b1]*(st[b2]-1);
	if (b1<=b2){
	  st[b1]--;
	  diagff[0] = index[st]-1;
	  st[b2]--;
	  diagff[1] = index[st]-1;
	  st[b1]++;
	  diagff[2] = index[st]-1;
	  st[b2]++;
	} else {
	  diagff[0] = -1;
	  diagff[1] = -1;
	  diagff[2] = -1;
	}
	soff[i].push_back(diagff);
	prff[i].push_back(prefact);
      }
    }
  }

  for (vint i=0; i<base.size(); i++){
    for (int b1=0; b1<baths; b1++){
      for (int b2=0; b2<baths; b2++){
	vector<int> st = base[i];
	vector<int> diagbb(3), diagff(3), diagbf(3), diagfb(3);
	// Only half of diagrams need to be calculated
	// The second half is equivalent due to symmetry
	// Only factor of 2 is added to the prefactor
	int prefact = 2*st[b2]*(Ns[b1]-st[b1]);
	st[b1]++;
	diagbf[0] = index[st]-1;
	st[b2]--;
	diagbf[1] = index[st]-1;
	st[b1]--;
	diagbf[2] = index[st]-1;
	st[b2]++;
	sobf[i].push_back(diagbf);
	prbf[i].push_back(prefact);
	
	prefact = 2*st[b1]*(Ns[b2]-st[b2]);
	st[b1]--;
	diagfb[0] = index[st]-1;
	st[b2]++;
	diagfb[1] = index[st]-1;
	st[b1]++;
	diagfb[2] = index[st]-1;
	st[b2]--;
	sofb[i].push_back(diagfb);
	prfb[i].push_back(prefact);
      }
    }
  }

  for (vint i=0; i<base.size(); i++){
    for (int b1=0; b1<baths; b1++){
      for (int b2=0; b2<baths; b2++){
	vector<int> st = base[i];
	vector<int> diagbb(3), diagff(3), diagbf(3), diagfb(3);
	int prefact = (b1!=b2) ? (Ns[b1]-st[b1])*(Ns[b2]-st[b2]) : (Ns[b1]-st[b1])*(Ns[b2]-st[b2]-1);
	prefact *= 2*degeneracy[i];
	prefact /= Ns[b1];
	st[b1]++;
	diagbb[0] = index[st]-1;
	st[b2]++;
	diagbb[1] = index[st]-1;
	st[b1]--;
	diagbb[2] = index[st]-1;
	st[b2]--;
	sogl[i].push_back(diagbb);
	prgl[i].push_back(prefact);
      }
    }
  }

}

// Prints out the information calculated above
void Print(ostream& stream, int baths, const vector<vector<int> >& base, const vector<int>& Mtot, const vector<int>& degeneracy)
{
  stream<<"#"<<setw(2)<<"i"<<"  "<<setw(baths*2)<<"state"<<"   "<<setw(6)<<"Mtot"<<setw(6)<<"Deg"<<endl;
  for(vint i=0; i<base.size(); i++){
    stream<<setw(3)<<i<<"  ";
    for (vint j=0; j<base[i].size(); j++) stream<<setw(2)<<base[i][j];
    stream<<"   "<<setw(4)<<Mtot[i]<<setw(7)<<degeneracy[i]<<endl;
  }
}

void Print2O(ostream& stream, int baths, int base_size, const vector<vector<vector<int> > >& sodiag){
  for(int i=0; i<base_size; i++){
    stream<<setw(3)<<i<<"  ";
    for (int b1=0; b1<baths; b1++)
      for (int b2=0; b2<baths; b2++){
	vector<int> diag = sodiag[i][b1*baths+b2];
	bool diag_exists = (diag[0]>=0) && (diag[1]>=0) && (diag[2]>=0);
	for (int j=0; j<3; j++) stream<<setw(2)<<(diag_exists ? diag[j] : -1)<<" ";
	stream<<"  ";
      }
    stream<<endl;
  }
}
void Print2O(ostream& stream, int baths, int base_size, const vector<vector<vector<int> > >& sodiag, const vector<vector<int> >& pref){
  for(int i=0; i<base_size; i++){
    stream<<setw(3)<<i<<"  ";
    for (int b1=0; b1<baths; b1++)
      for (int b2=0; b2<baths; b2++){
	int prefact = pref[i][b1*baths+b2];
	vector<int> diag = sodiag[i][b1*baths+b2];
	bool diag_exists = (diag[0]>=0) && (diag[1]>=0) && (diag[2]>=0);
	stream<<setw(8)<<(diag_exists? prefact : 0)<<" ";
	stream<<"  ";
      }
    stream<<endl;
  }
}

void Print(ostream& stream, int baths, const vector<vector<int> >& base, const vector<int>& Ns, const vector<int>& Mtot,
	   const vector<int>& degeneracy, const vector<vector<int> >& ncab, const vector<vector<int> >& ncaf,
	   const vector<vector<vector<int> > >& sobb, const vector<vector<vector<int> > >& sobf,
	   const vector<vector<vector<int> > >& sofb, const vector<vector<vector<int> > >& soff,
	   const vector<vector<vector<int> > >& sogl){
  stream<<"#"<<setw(2)<<"i"<<"  "<<setw(baths*2)<<"state"<<"  "<<setw(4)<<"Mtot"<<setw(7)<<"Deg"<<"  ";
  stream<<setw(5*baths)<<"ncab"<<setw(5*baths)<<"prefactb"<<setw(5*baths)<<"ncaf"<<setw(5*baths)<<"prefactf";
  stream<<setw(5*baths)<<"prefactG_loc"<<endl;
  for(vint i=0; i<base.size(); i++){
    stream<<setw(3)<<i<<"  ";
    for (int j=0; j<baths; j++) stream<<setw(2)<<base[i][j];
    stream<<"  "<<setw(4)<<Mtot[i]<<setw(7)<<degeneracy[i]<<"  ";
    for (int j=0; j<baths; j++) stream<<setw(5)<<ncab[i][j];
    for (int j=0; j<baths; j++) stream<<setw(5)<<Ns[j]-base[i][j];
    for (int j=0; j<baths; j++) stream<<setw(5)<<ncaf[i][j];
    for (int j=0; j<baths; j++) stream<<setw(5)<<base[i][j];
    for (int j=0; j<baths; j++) stream<<setw(5)<<(Ns[j]>0 ? degeneracy[i]*(Ns[j]-base[i][j])/Ns[j]: 0);
    stream<<endl;
  }
  stream<<"# 2-nd order back-back"<<endl;
  Print2O(stream, baths, base.size(), sobb);
  stream<<"# 2-nd order back-forward"<<endl;
  Print2O(stream, baths, base.size(), sobf);
//   stream<<"# 2-nd order forward-back"<<endl;
//   Print2O(stream, baths, base.size(), sofb);
  stream<<"# 2-nd order forward-forward"<<endl;
  Print2O(stream, baths, base.size(), soff);
  stream<<"# 2-nd order local-gf"<<endl;
  Print2O(stream, baths, base.size(), sogl);
}

void Print(ostream& out, int baths, int base_size,
	   const vector<vector<vector<int> > >& sobb, const vector<vector<vector<int> > >& sobf,
	   const vector<vector<vector<int> > >& sofb, const vector<vector<vector<int> > >& soff, const vector<vector<vector<int> > >& sogl,
	   const vector<vector<int> >& prbb, const vector<vector<int> >& prbf,
	   const vector<vector<int> >& prfb, const vector<vector<int> >& prff, const vector<vector<int> >& prgl)
{
  
  out<<"# 2-nd order back-back"<<endl;
  Print2O(out, baths, base_size, sobb, prbb);
  out<<"# 2-nd order back-forward"<<endl;
  Print2O(out, baths, base_size, sobf, prbf);
//   out<<"# 2-nd order forward-back"<<endl;
//   Print2O(out, baths, base_size, sofb, prfb);
  out<<"# 2-nd order forward-forward"<<endl;
  Print2O(out, baths, base_size, soff, prff);
  out<<"# 2-nd order local-gf"<<endl;
  Print2O(out, baths, base_size, sogl, prgl);
}

int main()
{
  int baths;
  cout<<" Give number of baths ... "<<endl;
  cin>>baths;
  if (!cin) {cerr<<"The value entered for number of baths is not valid\n";exit(1);}
    
  vector<int> Ns(baths);

  for (int i=0; i<baths; i++){
    cout<<" Give degeneracy of bath number "<<i<<" ... "<<endl;
    cin>>Ns[i];
    if (!cin) {cerr<<"The value entered for degeneracy of bath is not valid\n";exit(1);}
  }
  
  vector<vector<int> > base;
  vector<int> state(baths);
  generateBase(0,state,base,baths,Ns);
  cmp lessthan(Ns);
  sort(base.begin(),base.end(),lessthan);
  
  map<vector<int>,int> index;
  vector<int> Mtot(base.size());
  vector<int> degeneracy(base.size());
  
  SetUp(base, Ns, index, Mtot, degeneracy);
  Print(cout, baths, base, Mtot, degeneracy);
  
  cout<<" Printed are all possible states in this case. Which states would you like to keep?"<<endl;
  cout<<"  Enter your choice in format [xx-yy]+[,xx]+ (Example: 0-12,15,20)"<<endl;
  string str;
  cin>>str;

  list<int> small;
  parsInputChoice(str, small, base.size());
  
  vector<vector<int> > small_base(small.size());
  map<vector<int>,int> small_index;
  vector<int> small_Mtot(small_base.size());
  vector<int> small_degeneracy(small_base.size());
  vector<vector<int> > small_ncab(small_base.size());
  vector<vector<int> > small_ncaf(small_base.size());

  vector<vector<vector<int> > > sobb(small_base.size());
  vector<vector<vector<int> > > sobf(small_base.size());
  vector<vector<vector<int> > > sofb(small_base.size());
  vector<vector<vector<int> > > soff(small_base.size());
  vector<vector<vector<int> > > sogl(small_base.size());
  
  vector<vector<int> > prbb(small_base.size());
  vector<vector<int> > prbf(small_base.size());
  vector<vector<int> > prfb(small_base.size());
  vector<vector<int> > prff(small_base.size());
  vector<vector<int> > prgl(small_base.size());
  
  int l=0;
  for (list<int>::iterator i=small.begin(); i!=small.end(); i++,l++) small_base[l] = base[*i];
  
  SetUp(small_base, Ns, small_index, small_Mtot, small_degeneracy);
  CalcNCADiag(baths, small_base, Ns, small_index, small_degeneracy, small_ncab, small_ncaf);
  CalcUNCADiag_(baths, small_base, Ns, small_index, small_degeneracy, sobb, sobf, sofb, soff, sogl, prbb, prbf, prfb, prff, prgl);
  
  cout<<"Give output filename?"<<endl;
  string filename;
  cin>>filename;
  ofstream out(filename.c_str());
  out<<"# Input file for UNCA impurity solver. Do not change this file if you are not absolutely aware of what you are doing!"<<endl;
  out<<baths<<" ";
  for (int i=0; i<baths; i++) out<<Ns[i]<<" ";
  out<<small_base.size()<<"  # Number of baths it's degeneracy an number of all local states"<<endl;
  Print(out, baths, small_base, Ns, small_Mtot, small_degeneracy, small_ncab,  small_ncaf, sobb, sobf, sofb, soff, sogl);
  Print(out, baths, small_base.size(), sobb, sobf, sofb, soff, sogl, prbb, prbf, prfb, prff, prgl);

  cout<<"Output written on "<<filename<<endl;
  return 0;
}
