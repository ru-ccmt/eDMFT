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
#include <sstream>
#include <cmath>

using namespace std;
typedef vector<int>::size_type vint;

bool spin_orbit = 1.0;

// Just old good !
int Binomial(int n, int m)
{
  double Mf = 1;
  for (int i=2; i<=m; i++) Mf *= i;
  double r = 1;
  for (int i=n; i>=n-m+1; i--) r*=i;
  return static_cast<int>(r/Mf);
}
// Crucial function that recursively generates the local base for any given number of baths ans its degeneracy
void generateBase(int b, vector<int>& state, vector<vector<int> >& base, int baths, const vector<int>& Ns) 
{ 
  for(int i = 0; i<=Ns[b]; i++){
    state[b] = i;
    if (b<baths-1) 
      generateBase(b+1, state,base,baths,Ns);
    else{      
      base.push_back(state);
    }
  }
}

void GenerateBase(vector<vector<int> >& base, int baths)
{
  vector<int> state(baths);
  vector<int> nn(baths);
  for (int i=0; i<baths; i++) nn[i]=3;
  generateBase(0,state,base,baths,nn);
  for (vint i=0; i<base.size(); i++)
    for (vint j=0; j<baths; j++) base[i][j]--;
}

// This comparisson is needed to sort the base by the total number of electrons in the
// local state. Inside the subspace of the same number of particles the states are sorted
// such that firts bit is first filled. This is not necessary but just more convenient
// when there is only one particle per bath.
class cmp{
  const vector<int>& Mtot;
  const vector<vector<int> >& Nv;
  const vector<double>& Jc;
  const vector<vector<int> >& dblo;
public:
  cmp(const vector<int>& Mtot_, const vector<vector<int> >& Nv_, const vector<double>& Jc_,
      const vector<vector<int> >& dblo_) :
    Mtot(Mtot_), Nv(Nv_), Jc(Jc_), dblo(dblo_) {};
  bool operator()(int a, int b){
    if (Mtot[a]!=Mtot[b]) return Mtot[a]<Mtot[b];
    for (int k=0; k<Nv.size(); k++)
      if (Nv[k][a]!=Nv[k][b]) return Nv[k][a]>Nv[k][b];
    if (fabs(Jc[a]-Jc[b])>1e-6) return Jc[a]<Jc[b];
    for (int k=0; k<dblo.size(); k++)
      if (dblo[k][a]!=dblo[k][b]) return dblo[k][a]>dblo[k][b];
    return 1;
  }
};

// Given the base and bath degeneracy calculates index to the base, total number of particles
// in the state, degeneracy of each state and NCA diagrams.
void SetUp(const vector<vector<int> >& base, const vector<int>& Ns, map<vector<int>,int>& index,
	   vector<int>& Mtot, vector<int>& degeneracy){
  for(vint i=0; i<base.size(); i++){
    int nt=0; int dg=1;
    for (vint j=0; j<base[i].size(); j++){
      nt += abs(base[i][j]);
      //      dg *= Binomial(Ns[j],base[i][j]);
      dg=1;
    }
    Mtot[i] = nt;
    degeneracy[i] = dg;
    index[base[i]]=i+1;
  }
}
// Given the base and bath degeneracy calculates index to the base, total number of particles
// in the state, degeneracy of each state and NCA diagrams.
void SetUp1(const vector<vector<int> >& base, int baths, const vector<int>& Ns, vector<vector<int> >& Nv, vector<int>& Mtot,
	    vector<double>& Jc, vector<vector<int> >& dblo){
  bool first_time=true;
  for (vint i=0; i<base.size(); i++){
    
    int nst=0;
    int dbl=0;
    Mtot[i]=0;
    for (int k=0; k<Ns.size(); k++){
      int dbl1=0;
      int b1=0;
      for (int j=nst; j<nst+Ns[k]; j++){
	b1 += abs(base[i][j]);
	dbl1 += (base[i][j]==2) ? 1 : 0;
      }
      Nv[k][i] = b1;
      Mtot[i] += b1;
      nst+=Ns[k];
      dblo[k][i] = dbl1;
      dbl += dbl1;
    }

    double cJ1 = 1.5*dbl;
    
    double cJ2=0;
    for (int j1=0; j1<baths; j1++){
      int sz1 = abs(base[i][j1]<2) ? base[i][j1] : 0;
      for (int j2=j1+1; j2<baths; j2++){
	int sz2 = abs(base[i][j2]<2) ? base[i][j2] : 0;
	cJ2 -= 0.5*sz1*sz2;
      }
    }

    double cJ3=0;
    if (Mtot[i]==6 && cJ1+cJ2==-7.5 && first_time){cJ3=-spin_orbit; first_time=false;}
      
    Jc[i] = cJ1+cJ2+cJ3;
  }
}

void Resort(vector<vector<int> >& base, vector<int>& Mtot, vector<double>& Jc, vector<vector<int> >& Nv,
	    vector<vector<int> >& dblo)
{
  // Sort
  vector<int> bindex(base.size());
  for (int i=0; i<base.size(); i++) bindex[i]=i;
  cmp lessthan(Mtot,Nv,Jc,dblo);
  sort(bindex.begin(),bindex.end(),lessthan);
  // Rearrange
  vector<vector<int> > tmpbase(base.size());
  for (int i=0; i<base.size(); i++) tmpbase[i] = base[bindex[i]];
  base = tmpbase;
  vector<int> tmpMtot(base.size());
  for (int i=0; i<base.size(); i++) tmpMtot[i] = Mtot[bindex[i]];
  Mtot = tmpMtot;
  vector<double> tmpJc(base.size());
  for (int i=0; i<base.size(); i++) tmpJc[i] = Jc[bindex[i]];
  Jc = tmpJc;
  for (int k=0; k<Nv.size(); k++){
    for (int i=0; i<base.size(); i++) tmpMtot[i] = Nv[k][bindex[i]];
    Nv[k] = tmpMtot;
  }
}

// Prints out the information calculated above
void Print(ostream& stream, int baths, const vector<vector<int> >& base, const vector<int>& Mtot,
	   const vector<double>& Jc, const vector<vector<int> >& Nv, const vector<int>& sindex,
	   const vector<vector<int> >& ncab, vector<vector<int> >& ncaf)
{
  stream<<"#"<<setw(3)<<"i"<<"  "<<setw(3)<<"vi"<<setw(baths*3)<<"state"<<"   "<<setw(5);
  for (int k=0; k<Nv.size(); k++) stream<<"N"<<k+1<<setw(3);
  stream<<setw(5)<<"Mtot"<<setw(7)<<"Jc"<<endl;
  
  for(vint i=0; i<base.size(); i++){
    stream<<setw(4)<<i<<"  "<<setw(3)<<sindex[i]<<" ";
    for (vint j=0; j<base[i].size(); j++) stream<<setw(3)<<base[i][j];
    stream<<"   ";
    for (int k=0; k<Nv.size(); k++) stream<<setw(4)<<Nv[k][i];
    stream<<setw(5)<<Mtot[i]<<setw(7)<<Jc[i]<<"  ";
    for (int j=0; j<ncab[i].size(); j++) stream<<setw(6)<<ncab[i][j];
    for (int j=0; j<ncaf[i].size(); j++) stream<<setw(6)<<ncaf[i][j];
    stream<<endl;
  }
}


void Determin_Nonequivalent_States(int Ns_size, int base_size, const vector<vector<int> >& Nv,
				   const vector<double>& Jc,
				   const vector<int>& Mtot, list<int>& sMtot,
				   vector<int>& sindex, int& sbase_size)
{
  // Determins to which group of states each basis state corresponds
  vector<int> iNv(Ns_size);
  for (int k=0; k<Ns_size; k++) iNv[k]=Nv[k][0];
  double iJc = Jc[0];
  sMtot.push_back(Mtot[0]);
  int ll=0;
  for (int j=0; j<base_size; j++){
    bool equal=true;
    if (abs(iJc-Jc[j])>1e-6) equal=false;
    for (int k=0; k<Ns_size; k++) if (iNv[k]!=Nv[k][j]) equal=false;
    if (!equal){
      ll++;
      for (int k=0; k<Ns_size; k++) iNv[k]=Nv[k][j];
      iJc=Jc[j];
      // poskus!
      sMtot.push_back(Mtot[j]);
    }
    sindex[j]=ll;
  }
  sbase_size = ll+1;
}

void  Determin_Quantities_in_New_base(int base_size, int sbase_size, const vector<int>& sindex,
				      const vector<int>& Mtot, const vector<double>& Jc, const vector<vector<int> >& Nv,
				      vector<int>& sMtot, vector<int>& sdegen, vector<double>& sJc, vector<vector<int> >& sNv)
{
  for (int i=0; i<sbase_size; i++) sdegen[i]=0;
  for (int k=0; k<sNv.size(); k++) sNv[k].resize(sbase_size);
  for (int i=0; i<base_size; i++){
    if (sindex[i]<0 || sindex[i]>sbase_size) {cerr<<"sindex has wrong value!"<<endl; exit(1);}
    sMtot[sindex[i]] = Mtot[i];
    sJc[sindex[i]] = Jc[i];
    for (int k=0; k<Nv.size(); k++) sNv[k][sindex[i]] = Nv[k][i];
    sdegen[sindex[i]]++;
  }
}

bool add_electron(int& st, int ud){
  if (st==2 || st==ud) return false;
  if (st==0) {st=ud; return true;}
  st=2; return true;
}
bool rem_electron(int& st, int ud){
  if (st==ud) {st=0; return true;}
  if (st==2)  {st=-ud; return true;}
  return false;
}

void Calc_NCA_Diag(const vector<vector<int> >& base, int sbase_size, int Ns_size, int baths, const vector<int>& nsin,
		   const vector<int>& sindex, map<vector<int>,int>& child_index,
		   vector<vector<int> >& ncab, vector<vector<int> >& ncaf,
		   vector<vector<map<int,int> > >& sncab, vector<vector<map<int,int> > >& sncaf)
{
  for (vint i=0; i<sbase_size; i++) sncab[i].resize(Ns_size);
  for (vint i=0; i<sbase_size; i++) sncaf[i].resize(Ns_size);
  
  for (vint i=0; i<base.size(); i++){
    int sind = sindex[i];
    for (int b=0; b<baths; b++){
      int ib = nsin[b];
      vector<int> state, state0 = base[i];

      child_index[state0];
      
      for (int ud=-1; ud<=1; ud+=2){
	state = state0;
	int indb=-1;
	if (add_electron(state[b], ud)){
	  indb = child_index[state]-1;
	  sncab[sind][ib][indb]++;
	}
	ncab[i].push_back(indb);
      }

      for (int ud=-1; ud<=1; ud+=2){
	state=state0;
	int indf=-1;
	if (rem_electron(state[b],ud)){
	  indf = child_index[state]-1;
	  sncaf[sind][ib][indf]++;
	}
	ncaf[i].push_back(indf);
      }
    }
  }
}

void Print_NCA_Diag(ostream& out, int sbase_size, int valence_size, const vector<int>& Ns, const vector<int>& sdegen, const vector<vector<int> >& sNv,
		    const vector<int>& sMtot, const vector<double>& sJc,
		    const vector<vector<map<int,int> > >& sncab, const vector<vector<map<int,int> > >& sncaf)
{
  out<<"# Input file for NCA impurity solver. Do not change this file if you are not absolutely aware of what you are doing!"<<endl;
  out<<Ns.size()<<" ";
  for (int i=0; i<Ns.size(); i++) out<<2*Ns[i]<<" ";
  out<<valence_size<<" "<<sbase_size-valence_size<<" ";
  out<<"# Number of baths it's degeneracy an number of local valence and local core states"<<endl;
  out<<setw(4)<<"#"<<" ";
  for (int k=0; k<sNv.size(); k++) { stringstream Nk; Nk<<"N"<<k; out<<setw(3)<<Nk.str();}
  out<<setw(4)<<"Mtot"<<setw(5)<<"deg"<<setw(6)<<"Jc"<<"  ";;
  for (int p=0; p<Ns.size(); p++) out<<setw(2)<<"#b"<<" ";
  for (int p=0; p<Ns.size(); p++) out<<setw(2)<<"#f"<<" ";
  out<<endl;
  
  for (int i=0; i<sbase_size; i++){
    out<<setw(4)<<i<<" ";
    for (int k=0; k<sNv.size(); k++) out<<setw(3)<<sNv[k][i];
    out<<setw(4)<<sMtot[i]<<setw(5)<<sdegen[i]<<setw(6)<<sJc[i]<<"  ";
    
    double deg=sdegen[i];
    
    for (int p=0; p<Ns.size(); p++) out<<setw(2)<<sncab[i][p].size()<<" ";
    for (int p=0; p<Ns.size(); p++) out<<setw(2)<<sncaf[i][p].size()<<" ";
    
    for (int p=0; p<Ns.size(); p++)
      for (map<int,int>::const_iterator l=sncab[i][p].begin(); l!=sncab[i][p].end(); l++)
	out<<setw(6)<<l->second<<" x "<<setw(4)<<left<<l->first<<right;
    //	out<<setw(6)<<(double)(l->second)/(2*Ns[p])<<","<<setw(6)<<l->second/deg<<" x "<<setw(4)<<left<<l->first<<right;
	
    for (int p=0; p<Ns.size(); p++)
      for (map<int,int>::const_iterator l=sncaf[i][p].begin(); l!=sncaf[i][p].end(); l++)
	out<<setw(6)<<l->second<<" x "<<setw(4)<<left<<l->first<<right;
    //	out<<setw(6)<<l->second/deg<<" x "<<setw(4)<<left<<l->first<<right;
    
    out<<endl;
  }
}


void SimplifyDiagrams(const vector<int>& N0, int sbase_size, const vector<int>& Ns, const vector<vector<int> >& sNv, const vector<int>& sdegen,
		      const vector<int>& sMtot, const vector<double>& sJc, const vector<vector<map<int,int> > >& sncab, const vector<vector<map<int,int> > >& sncaf,
		      int& ssbase_size, int& svalence_size, int valence_start, int valence_end, vector<vector<int> >& ssNv, vector<int>& ssdegen, vector<int>& ssMtot,
		      vector<double>& ssJc, vector<vector<map<int,int> > >& ssncab, vector<vector<map<int,int> > >& ssncaf)
{
  vector<int> iNv(Ns.size());
  vector<int> ssindex(sbase_size);
  for (int k=0; k<Ns.size(); k++) iNv[k]=sNv[k][0];
  int ll=0;
  for (int j=0; j<sbase_size; j++){
    bool equal=true;
    for (int k=0; k<Ns.size(); k++) if (iNv[k]!=sNv[k][j]) equal=false;
    bool sMtot_in_N0=false;
    for (vector<int>::const_iterator l=N0.begin(); l!=N0.end(); l++) if (*l==sMtot[j]) sMtot_in_N0=true;
    if (!equal || sMtot_in_N0){
      ll++;
      for (int k=0; k<Ns.size(); k++) iNv[k]=sNv[k][j];
    }
    ssindex[j]=ll;
  }
  ssbase_size = ll+1;

  ssdegen.resize(ssbase_size);
  ssMtot.resize(ssbase_size);
  ssJc.resize(ssbase_size);
  ssncab.resize(ssbase_size);
  ssncaf.resize(ssbase_size);
  for (int i=0; i<ssbase_size; i++) {ssdegen[i]=0; ssJc[i]=0;}
  for (int k=0; k<Ns.size(); k++) ssNv[k].resize(ssbase_size);
  for (int i=0; i<ssbase_size; i++) {
    ssncab[i].resize(Ns.size());
    ssncaf[i].resize(Ns.size());
  }
  for (int j=0; j<sbase_size; j++){
    int nst = ssindex[j];
    ssdegen[nst] += sdegen[j];
    ssMtot[nst] = sMtot[j];
    for (int k=0; k<Ns.size(); k++) ssNv[k][nst] = sNv[k][j];

    for (int p=0; p<Ns.size(); p++){
      for (map<int,int>::const_iterator l=sncab[j][p].begin(); l!=sncab[j][p].end(); l++){
	int ncab = ssindex[l->first];
	ssncab[nst][p][ncab] += l->second;
      }
      for (map<int,int>::const_iterator l=sncaf[j][p].begin(); l!=sncaf[j][p].end(); l++){
	int ncaf = ssindex[l->first];
	ssncaf[nst][p][ncaf] += l->second;
      }
    }
  }
  for (int j=0; j<sbase_size; j++)
    if (ssdegen[ssindex[j]]==sdegen[j]) ssJc[ssindex[j]]=sJc[j];

  svalence_size=0;
  for (svalence_size=0; svalence_size<ssbase_size; svalence_size++) if (ssMtot[svalence_size]>valence_end || ssMtot[svalence_size]<valence_start) break;
}

void debugSimple(int sbase_size, const vector<int>& Ns, const vector<vector<int> >& sNv, const vector<int>& sdegen, const vector<int>& sMtot,
		 const vector<vector<map<int,int> > >& sncab, const vector<vector<map<int,int> > >& sncaf)
/*** Just for debugging *****************************************/
{
  vector<int> iNv(Ns.size());
  vector<int> ssindex(sbase_size);
  for (int k=0; k<Ns.size(); k++) iNv[k]=sNv[k][0];
  int ll=0;
  for (int j=0; j<sbase_size; j++){
    bool equal=true;
    for (int k=0; k<Ns.size(); k++) if (iNv[k]!=sNv[k][j]) equal=false;
    if (!equal){
      ll++;
      for (int k=0; k<Ns.size(); k++) iNv[k]=sNv[k][j];
    }
    ssindex[j]=ll;
  }
  int ssbase_size = ll+1;

  vector<int> ssdegen(ssbase_size);
  vector<int> ssMtot(ssbase_size);
  vector<vector<int> > ssNv(Ns.size());
  for (int i=0; i<ssbase_size; i++) ssdegen[i]=0;
  for (int k=0; k<Ns.size(); k++) ssNv[k].resize(ssbase_size);
  vector<vector<int> > ssncab(ssbase_size);
  vector<vector<double> > ppncab(ssbase_size);
  vector<vector<int> > ssncaf(ssbase_size);
  vector<vector<double> > ppncaf(ssbase_size);
  vector<vector<double> > ppncaloc(ssbase_size);
  for (int i=0; i<ssbase_size; i++) {
    ssncab[i].resize(Ns.size());
    ppncab[i].resize(Ns.size());
    ssncaf[i].resize(Ns.size());
    ppncaf[i].resize(Ns.size());
    ppncaloc[i].resize(Ns.size());
    for (int k=0; k<Ns.size(); k++) ppncaloc[i][k]=0;
  }
  for (int j=0; j<sbase_size; j++){
    ssdegen[ssindex[j]] += sdegen[j];
    ssMtot[ssindex[j]] = sMtot[j];
    for (int k=0; k<Ns.size(); k++) ssNv[k][ssindex[j]] = sNv[k][j];
    
    for (int p=0; p<Ns.size(); p++){
      int ncab=-1;
      int prencab=0;
      for (map<int,int>::const_iterator l=sncab[j][p].begin(); l!=sncab[j][p].end(); l++){
	ncab = ssindex[l->first];
	prencab += l->second;
      }
      ssncab[ssindex[j]][p]=ncab;
      ppncab[ssindex[j]][p]=(double)prencab/sdegen[j];
      ppncaloc[ssindex[j]][p]+=(double)prencab/(2*Ns[p]);

      int ncaf=-1;
      int prencaf=0;
      for (map<int,int>::const_iterator l=sncaf[j][p].begin(); l!=sncaf[j][p].end(); l++){
	ncaf = ssindex[l->first];
	prencaf += l->second;
      }
      ssncaf[ssindex[j]][p]=ncaf;
      ppncaf[ssindex[j]][p]=(double)prencaf/sdegen[j];
    }
  }

  for (int i=0; i<ssbase_size; i++){
    cout<<setw(3)<<i<<"  ";
    for (int k=0; k<Ns.size(); k++) cout<<setw(2)<<ssNv[k][i];
    cout<<"  "<<setw(4)<<ssMtot[i]<<setw(7)<<ssdegen[i]<<"  ";
    for (int p=0; p<Ns.size(); p++) cout<<setw(5)<<ssncab[i][p];
    for (int p=0; p<Ns.size(); p++) cout<<setw(5)<<ppncab[i][p];
    for (int p=0; p<Ns.size(); p++) cout<<setw(5)<<ssncaf[i][p];
    for (int p=0; p<Ns.size(); p++) cout<<setw(5)<<ppncaf[i][p];
    for (int p=0; p<Ns.size(); p++) cout<<setw(5)<<ppncaloc[i][p];
    cout<<endl;
  }
}

void parsNs(const string& str, vector<int>& Ns)
{
  
  stringstream inp(str);

  list<int> lNs;
  int a; char ch;
  while (inp){
    inp>>a;
    if (inp) lNs.push_back(a);
    else break;
    inp>>ch;
  }
  Ns.resize(lNs.size());
  int k=0;
  for (list<int>::iterator l=lNs.begin(); l!=lNs.end(); l++,k++) Ns[k] = *l;
}

int main()
{
  
  int baths;
  cout<<" Give number of bands (for d electrons 5, for f 7)... "<<endl;
  cin>>baths;
  if (!cin) {cerr<<"The value entered for number of baths is not valid\n";exit(1);}

  cout<<" How many multiplets form these bands (for example 7=3+3+1 for f bans in QS)... "<<endl;
  cout<<baths<<"=";
  string str;
  cin>>str;

  vector<int> Ns;
  parsNs(str,Ns);

  int tNs=0;
  for (int i=0; i<Ns.size(); i++) tNs += Ns[i];
  if (tNs!=baths) {cerr<<"The value entered for # of multiplets is not valid. The sum of these numbers must be equal to number of bands!\n";exit(1);}

  int valence_start, valence_end;
  cout<<"Which occupany of local level should be treated dynamically (as valence state)? "<<endl;
  cout<<"Minimum value is 0 and maximum value is "<<baths*2<<"! Specify it like: 0-"<<baths*2<<endl;
  cin>>str;
  vector<int> valenceSE;
  parsNs(str,valenceSE);
  if (valenceSE.size()!=2) {cerr<<"Format not recongized! Did you entered two integer numbers and - inbetween?"<<endl;}
  valence_start = valenceSE[0];
  valence_end = valenceSE[1];
  if (valence_start<0) valence_start=0;
  if (valence_end>2*baths) valence_end = 2*baths;
  
  
  clog<<"Number of bans is "<<baths<<endl;
  for (int i=0; i<Ns.size(); i++){
    clog<<"Multiplet number "<<i+1<<" contains "<<Ns[i]<<" bands"<<endl;
  }
  clog<<"Press any key to continue...."<<endl;
  char ch;
  cin.get(ch);
  clog<<".....calculating WAIT...."<<endl;


  /*
  int baths=7;
  int nNs=3;
  vector<int> Ns(nNs);
  Ns[0]=3;
  Ns[1]=3;
  Ns[2]=1;
  */

  /*
  int baths=7;
  int nNs=2;
  vector<int> Ns(nNs);
  Ns[0]=4;
  Ns[1]=3;

  int valence_start = 5;
  int valence_end = 7;
  */
  
  /*
  int baths=7;
  int nNs=1;
  vector<int> Ns(nNs);
  Ns[0]=7;
  */
  
  vector<int> nsin(baths);
  int li=0;
  for (int i=0; i<Ns.size(); i++){
    for (int j=0; j<Ns[i]; j++,li++) nsin[li]=i;
  }

  vector<vector<int> > base;
  GenerateBase(base,baths);
  
  vector<int> Mtot(base.size());
  vector<double> Jc(base.size());
  vector<vector<int> > Nv(Ns.size());
  for (int k=0; k<Ns.size(); k++) Nv[k].resize(base.size());
  vector<vector<int> > dblo(Ns.size());
  for (int k=0; k<Ns.size(); k++) dblo[k].resize(base.size());
  SetUp1(base, baths, Ns, Nv, Mtot, Jc, dblo);

  Resort(base, Mtot, Jc, Nv, dblo);

  vector<int> sindex(base.size());
  int sbase_size;
  list<int> sMtot0;
  Determin_Nonequivalent_States(Ns.size(), base.size(), Nv, Jc, Mtot, sMtot0, sindex, sbase_size);


  vector<int> valence_sindex(sbase_size);
  int ii=0;
  int jj=0;
  for (list<int>::const_iterator i=sMtot0.begin(); i!=sMtot0.end(); ++i,++jj)
    if ((*i)>=valence_start && (*i)<=valence_end) valence_sindex[jj] = ii++;

  int valence_size = ii;
  jj=0;
  for (list<int>::const_iterator i=sMtot0.begin(); i!=sMtot0.end(); ++i,++jj)
    if (!((*i)>=valence_start && (*i)<=valence_end)) valence_sindex[jj] = ii++;

  vector<int> tmp_sindex(base.size());
  for (int i=0; i<base.size(); i++) tmp_sindex[i] = valence_sindex[sindex[i]];
  sindex = tmp_sindex;

  
  vector<int> sMtot(sbase_size);
  vector<double> sJc(sbase_size);
  vector<vector<int> > sNv(Nv.size());
  vector<int> sdegen(sbase_size);
  Determin_Quantities_in_New_base(base.size(), sbase_size, sindex, Mtot, Jc, Nv, sMtot, sdegen, sJc, sNv);

  map<vector<int>,int> child_index;
  for (int i=0; i<base.size(); i++) child_index[base[i]]=sindex[i]+1;

  vector<vector<int> > ncab(base.size());
  vector<vector<int> > ncaf(base.size());
  vector<vector<map<int,int> > > sncab(sbase_size);
  vector<vector<map<int,int> > > sncaf(sbase_size);

  Calc_NCA_Diag(base, sbase_size, Ns.size(), baths, nsin, sindex, child_index, ncab, ncaf, sncab, sncaf);

  // Here is the complicated total output!!
  //  Print(cout, baths, base, Mtot, Jc, Nv, sindex, ncab, ncaf);
  //  cout<<endl;

  ofstream out1("out1.cix");
  Print_NCA_Diag(cout, sbase_size, valence_size, Ns, sdegen, sNv, sMtot, sJc, sncab, sncaf);
  Print_NCA_Diag(out1, sbase_size, valence_size, Ns, sdegen, sNv, sMtot, sJc, sncab, sncaf);

  vector<int> N0(1);
  N0[0]=6;

  vector<int> ssdegen;
  vector<int> ssMtot;
  vector<double> ssJc;
  vector<vector<int> > ssNv(Ns.size());
  vector<vector<map<int,int> > > ssncab;
  vector<vector<map<int,int> > > ssncaf;
  int ssbase_size, svalence_size;
  SimplifyDiagrams(N0,sbase_size,Ns,sNv,sdegen,sMtot,sJc,sncab,sncaf,ssbase_size,svalence_size,valence_start,valence_end,ssNv,ssdegen,ssMtot,ssJc,ssncab,ssncaf);
  ofstream out2("out2.cix");
  Print_NCA_Diag(cout, ssbase_size, svalence_size, Ns, ssdegen, ssNv, ssMtot, ssJc, ssncab, ssncaf);
  Print_NCA_Diag(out2, ssbase_size, svalence_size, Ns, ssdegen, ssNv, ssMtot, ssJc, ssncab, ssncaf);
  
  return 0;
}
