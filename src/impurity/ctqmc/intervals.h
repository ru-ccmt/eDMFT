// @Copyright 2007 Kristjan Haule
// 
#include <map>
#include <deque>

using namespace std;

typedef multimap<double,int> Tinterval;

class IntervalIndex{
public: 
  int ifl;
  int type;
  int in;
  IntervalIndex(int ifl_, int type_, int in_) : ifl(ifl_), type(type_), in(in_){};
  IntervalIndex(int x=-1) : ifl(x), type(x), in(x){}; // This is for conversion from integer to default state of the class
  // This state of the index should never be used. The last slot in Operators (at t=beta) does not contain this index. This
  // default value -1 will be set in this case.
};
inline std::ostream& operator<< (std::ostream& stream, const IntervalIndex& r){
  int width = stream.width();
  stream << std::setw(width)<< r.ifl << " " << std::setw(width) << r.type << " " << std::setw(width) << r.in;
  return stream;
}


class nIntervals{
  function2D<dcomplex> _exp_e, _exp_s;
  double beta;
  function1D<double> _time_e, _time_s;
  function1D<int> index_e, index_s;
  function2D<dcomplex> b_exp_e, b_exp_s; // bosonic frequencies in case needed
  function1D<int> _btype_e, _btype_s;
  deque<int> Empty_e, Empty_s;
  const mesh1D* piom;  // fermionic frequencies
  const mesh1D* biom;  // bosonic frequencies
public:// members
  function1D<int> Nbtyp_s, Nbtyp_e;
  function1D<int> index_s_1;
  static const int cd = 0;
  static const int c = 1;
public:// methods  
  nIntervals() {};
  void SetUp(int Nmax, const mesh1D& iom, const mesh1D& iomb, double beta_, int dim);
  // inlined short methods to access class members
  int size() const {return _time_e.size()+_time_s.size();}
  int fullsize() const {return _time_e.fullsize()+_time_s.fullsize();}
  const funProxy<dcomplex>& exp_e(int i) const {return _exp_e[index_e[i]];}
  const funProxy<dcomplex>& exp_s(int i) const {return _exp_s[index_s[i]];}
  const funProxy<dcomplex>& exp(int type, int i) const {return type==cd ? _exp_s[index_s[i]] : _exp_e[index_e[i]];}
  double time_s(int i) const {return _time_s[index_s[i]];}
  double time_e(int i) const {return _time_e[index_e[i]];}
  double time(int type, int i) const {return type==cd ? _time_s[index_s[i]] : _time_e[index_e[i]];}
  
  double time_direct(int type, int i) const {return type==cd ? _time_s[i] : _time_e[i];}
  
  template <int boson_fermion>
  const funProxy<dcomplex>& exp_direct(int type, int i) const {cerr<<"Should not happen"<<endl; }//return NULL;}
  
  int btype_s(int i) const {return _btype_s[index_s[i]];}
  int btype_e(int i) const {return _btype_e[index_e[i]];}
  int btype(int type, int i) const {return type==cd ? _btype_s[index_s[i]] : _btype_e[index_e[i]];}
  int Nbtype(int type, int b) const {return type==cd ? Nbtyp_s[b] : Nbtyp_e[b];}
  void Find_is_ie(double t_start, int& is, double t_end, int& ie);
  pair<int,int> InsertExponents(double t_start, int is, int btyp_s, double t_end, int ie, int btyp_e);
  void RemoveExponents(int is, int ie);
  void MoveExponents(int type, double t_old, int i_old, double t_new, int i_new);
  int FindIndex(double time, double t_old, int type);
  int FindIndex_nearby(double t_new, int type, int i0);
  void print(ostream& stream);
  void TestOrder();
  int FindSuccessive(int type, int ie, int bfle);
  int FindWhichOneIs(int type, int bflx, int iop);
  double PreviousTime(int type, int to_move, bool add_beta) const;
  double NextTime(int type, int to_move, bool add_beta) const;
  void copy_data(const nIntervals& source);
  void check_index_1();
  pair<double,double> Closest_Start_Times(int bfl, double t_s) const;
  pair<double,double> Closest_End_Times(int bfl, double t_e) const;
  pair<double,double> Closest_Times(int type, int bfl, double t_current) const;
  void Closest_Times(pair<double,double>& ts, pair<int,int>& is, int type, int bfl, double t_current);
  bool CheckOnSegmentHere(int bfl, double t_start, double& next_t_start, double& next_t_end, double& previous_t_start, double& previous_t_end, long long istep);
  pair<int,int> NextTimeIndex_(int type, double t_current, int bfl);
  void Resize(int newNmax);
  friend void IntervalsExchangeExponents(int typea, nIntervals& intervala, int ia, int ia_new,int typeb, nIntervals& intervalb, int ib, int ib_new);
  void ReleaseSomeTempMemory();
};

template <>
const funProxy<dcomplex>& nIntervals::exp_direct<0>(int type, int i) const
{return type==cd ? b_exp_s[i] : b_exp_e[i];}

template <>
const funProxy<dcomplex>& nIntervals::exp_direct<1>(int type, int i)  const
{return type==cd ? _exp_s[i] : _exp_e[i];}

  

double nIntervals::PreviousTime(int type, int to_move, bool add_beta=false) const
{
  double dbt = add_beta ? common::beta : 0.0;// When it jumps around: should we add/subtract beta or not.
  if (type==c){
    int bfl = _btype_e[index_e[to_move]];
    for (int ii=to_move-1; ii>=0; ii--)
      if (_btype_e[index_e[ii]]==bfl) return _time_e[index_e[ii]];
    for (int ii=_time_e.size()-1; ii>=to_move; ii--)
      if (_btype_e[index_e[ii]]==bfl) return _time_e[index_e[ii]]-dbt;
  }else{
    int bfl = _btype_s[index_s[to_move]];
    for (int ii=to_move-1; ii>=0; ii--)
      if (_btype_s[index_s[ii]]==bfl) return _time_s[index_s[ii]];
    for (int ii=_time_s.size()-1; ii>=to_move; ii--)
      if (_btype_s[index_s[ii]]==bfl) return _time_s[index_s[ii]]-dbt;
  }
  return 0;
}

double nIntervals::NextTime(int type, int to_move, bool add_beta=false) const
{
  double dbt = add_beta ? common::beta : 0.0;// When it jumps around: should we add/subtract beta or not.
  if (type==c){
    int bfl = _btype_e[index_e[to_move]];
    for (int ii=to_move+1; ii<_time_e.size(); ii++)
      if (_btype_e[index_e[ii]]==bfl) return _time_e[index_e[ii]];
    for (int ii=0; ii<=to_move; ii++)
      if (_btype_e[index_e[ii]]==bfl) return _time_e[index_e[ii]]+dbt;
  }else{
    int bfl = _btype_s[index_s[to_move]];
    for (int ii=to_move+1; ii<_time_s.size(); ii++)
      if (_btype_s[index_s[ii]]==bfl) return _time_s[index_s[ii]];
    for (int ii=0; ii<=to_move; ii++)
      if (_btype_s[index_s[ii]]==bfl) return _time_s[index_s[ii]]+dbt;
  }
  return 0;
}

pair<double,double> nIntervals::Closest_Start_Times(int bfl, double t_s) const
{
  int i_previous=-1;
  int i_first=-1;
  int i_last=-1;

  for (int ii=0; ii<_time_s.size(); ii++){
    int is = index_s[ii];
    if (_btype_s[is]==bfl){
      if (i_last<0) i_first=ii;
      i_last=ii;
      if (_time_s[is]>t_s) break;
      i_previous=ii;
    }
  }

  if (i_first<0)                  // No such type in the list
    return pair<double,double>(-common::beta,-common::beta);
  
  if (i_previous<0 && i_last>=0){ // The new time is earlier than any other time in the list
    int ii=_time_s.size()-1;
    for (; ii>=0; --ii)
      if (_btype_s[index_s[ii]]==bfl) break;
    return pair<double,double>(_time_s[index_s[ii]]-common::beta, _time_s[index_s[i_last]]);
  }
  if (i_previous==i_last)        // The new time is later than any other time in the list
    return pair<double,double>(_time_s[index_s[i_previous]],_time_s[index_s[i_first]]+common::beta);
  
  return pair<double,double>(_time_s[index_s[i_previous]],_time_s[index_s[i_last]]);
}

pair<double,double> nIntervals::Closest_End_Times(int bfl, double t_e) const
{
  int i_previous=-1;
  int i_first=-1;
  int i_last=-1;

  for (int ii=0; ii<_time_e.size(); ii++){
    int ie = index_e[ii];
    if (_btype_e[ie]==bfl){
      if (i_last<0) i_first=ii;
      i_last=ii;
      if (_time_e[ie]>t_e) break;
      i_previous=ii;
    }
  }

  if (i_first<0)                  // No such type in the list
    return pair<double,double>(-common::beta,-common::beta);
  
  if (i_previous<0 && i_last>=0){ // The new time is earlier than any other time in the list
    int ii=_time_e.size()-1;
    for (; ii>=0; --ii)
      if (_btype_e[index_e[ii]]==bfl) break;
    return pair<double,double>(_time_e[index_e[ii]]-common::beta, _time_e[index_e[i_last]]);
  }
  if (i_previous==i_last)        // The new time is later than any other time in the list
    return pair<double,double>(_time_e[index_e[i_previous]],_time_e[index_e[i_first]]+common::beta);
  
  return pair<double,double>(_time_e[index_e[i_previous]],_time_e[index_e[i_last]]);
}

inline pair<double,double> nIntervals::Closest_Times(int type, int bfl, double t_current) const
{
  if (type==nIntervals::cd){
    return Closest_Start_Times(bfl, t_current);
  }else{
    return Closest_End_Times(bfl, t_current);
  }
}

void nIntervals::Closest_Times(pair<double,double>& times, pair<int,int>& index, int type, int bfl, double t_current)
{
  int i_previous=-1;
  int i_first=-1;
  int i_last=-1;
  int i_number=0;
  if (type==cd){
    for (int ii=0; ii<_time_s.size(); ii++){
      int is = index_s[ii];
      if (_btype_s[is]==bfl){
	if (i_last<0) i_first=ii;
	i_last=ii;
	if (_time_s[is]>t_current) break;
	i_previous=ii;
	i_number++;
      }
    }
    if (i_first<0){                  // No such type in the list
      times=make_pair(-common::beta,-common::beta);
      index=make_pair(-1,-1);
      return;
    }
    if (i_previous<0 && i_last>=0){ // The new time is earlier than any other time in the list
      int ii=_time_s.size()-1;
      for (; ii>=0; --ii)
	if (_btype_s[index_s[ii]]==bfl) break;
      times=make_pair(_time_s[index_s[ii]]-common::beta, _time_s[index_s[i_last]]);
      index=make_pair(ii,i_last);
      return;
    }
    if (i_previous==i_last){        // The new time is later than any other time in the list
      times=make_pair(_time_s[index_s[i_previous]],_time_s[index_s[i_first]]+common::beta);
      index=make_pair(i_previous,i_first);
      return;
    }
    times=make_pair(_time_s[index_s[i_previous]],_time_s[index_s[i_last]]);
    index=make_pair(i_previous,i_last);
    return;
  }else{
    for (int ii=0; ii<_time_e.size(); ii++){
      int ie = index_e[ii];
      if (_btype_e[ie]==bfl){
	if (i_last<0) i_first=ii;
	i_last=ii;
	if (_time_e[ie]>t_current) break;
	i_previous=ii;
	i_number++;
      }
    }
    if (i_first<0){                  // No such type in the list
      times=make_pair(-common::beta,-common::beta);
      index=make_pair(-1,-1);
      return;
    }
    if (i_previous<0 && i_last>=0){ // The new time is earlier than any other time in the list
      int ii=_time_e.size()-1;
      for (; ii>=0; --ii)
	if (_btype_e[index_e[ii]]==bfl) break;
      times=make_pair(_time_e[index_e[ii]]-common::beta, _time_e[index_e[i_last]]);
      index=make_pair(ii,i_last);
      return;
    }
    if (i_previous==i_last){        // The new time is later than any other time in the list
      times=make_pair(_time_e[index_e[i_previous]],_time_e[index_e[i_first]]+common::beta);
      index=make_pair(i_previous,i_first);
      return;
    }
    times=make_pair(_time_e[index_e[i_previous]],_time_e[index_e[i_last]]);
    index=make_pair(i_previous,i_last);
    return;
  }
}


bool nIntervals::CheckOnSegmentHere(int bfl, double t_current, double& next_t_start, double& next_t_end, double& previous_t_start, double& previous_t_end, long long istep)
{
  pair<double,double> ts = Closest_Start_Times(bfl, t_current); // ts.first is creation operator just before t_start
  if (ts.second<=-common::beta+1e-15){ // No kinks yet of this type
    next_t_start=t_current+common::beta;
    next_t_end=t_current+common::beta;
    return false;
  }
  next_t_start = ts.second;
  pair<double,double> te = Closest_End_Times(bfl, t_current);
  next_t_end=te.second;
  previous_t_start = ts.first;
  previous_t_end = te.first;
  if (te.first<ts.first && te.second<ts.second) return true;
  //if (te.second==t_current || ts.first==t_current) return true;// should also not be exactly equal to the neighbors
  else return false;
}

pair<int,int> nIntervals::NextTimeIndex_(int type, double t_current, int bfl)
{
  int i_previous=-1;
  int i_first=-1;
  int i_last=-1;
  int i_number=0;
  if (type==cd){
    for (int ii=0; ii<_time_s.size(); ii++){
      int is = index_s[ii];
      if (_btype_s[is]==bfl){
	if (i_last<0) i_first=ii;
	i_last=ii;
	if (_time_s[is]>t_current) break;
	i_previous=ii;
	i_number++;
      }
    }
  }else{
    for (int ii=0; ii<_time_e.size(); ii++){
      int ie = index_e[ii];
      if (_btype_e[ie]==bfl){
	if (i_last<0) i_first=ii;
	i_last=ii;
	if (_time_e[ie]>t_current) break;
	i_previous=ii;
	i_number++;
      }
    }
  }
  if (i_first<0)                  // No such type in the list
    return pair<int,int>(-1,-1);
  if (i_previous==i_last)        // The new time is later than any other time in the list
    return pair<int,int>(i_first,0);
  return pair<int,int>(i_last,i_number);
}


int nIntervals::FindSuccessive(int type, int ix, int bflx)
{
  int ii=0;
  if (type==c){
    for (int i=0; i<ix; i++) if (_btype_e[index_e[i]]==bflx) ii++;
  } else{
    for (int i=0; i<ix; i++) if (_btype_s[index_s[i]]==bflx) ii++;
  }
  return ii;
}


int nIntervals::FindWhichOneIs(int type, int bflx, int iop){
  int ii=0;
  if (type==c){
    for (int i=0; i<_time_e.size(); i++){
      if (_btype_e[index_e[i]]==bflx){
	if (ii==iop) return i;
	ii++;
      }
    }
  } else{
    for (int i=0; i<_time_s.size(); i++){
      if (_btype_s[index_s[i]]==bflx){
	if (ii==iop) return i;
	ii++;
      }
    }
  }
  cout<<"ERROR: Could not find this element in FindWhichOneIs! "<<ii<<endl;
  return -1;
}


void nIntervals::SetUp(int Nmax, const mesh1D& iom, const mesh1D& iomb, double beta_, int dim=1)
{
  beta = beta_;
  int nom = iom.size();
  int nomb = iomb.size();
  _time_e.resize(Nmax); _time_e.resize_virtual(0);
  _time_s.resize(Nmax); _time_s.resize_virtual(0);
  if (!common::Qsvd){
    _exp_e.resize(Nmax,nom);
    _exp_s.resize(Nmax,nom);
    b_exp_e.resize(Nmax,nomb);
    b_exp_s.resize(Nmax,nomb);
  }
  index_e.resize(Nmax);
  index_s.resize(Nmax);
  index_s_1.resize(Nmax);
  Empty_e.resize(Nmax);
  Empty_s.resize(Nmax);
  _btype_s.resize(Nmax);
  _btype_e.resize(Nmax);
  Nbtyp_s.resize(dim);
  Nbtyp_s=0;
  Nbtyp_e.resize(dim);
  Nbtyp_e=0;
  piom = &iom;
  biom = &iomb;
  for (int i=0; i<Nmax; i++) Empty_e[i]=i;
  for (int i=0; i<Nmax; i++) Empty_s[i]=i;
};

void nIntervals::Resize(int newNmax)
{
  int Nmax = index_e.size();
  int Nk = _time_e.size();
  _time_e.resize_preserve(newNmax);
  //_time_e.resize_virtual(Nk);
  Nk = _time_s.size();
  _time_s.resize_preserve(newNmax);
  //_time_s.resize_virtual(Nk);
  
  if (!common::Qsvd){
    int nom = piom->size();
    int nomb = biom->size();
    _exp_e.resize_preserve(newNmax,nom);   _exp_e.resize_virtual(newNmax,nom);
    _exp_s.resize_preserve(newNmax,nom);   _exp_s.resize_virtual(newNmax,nom);
    b_exp_e.resize_preserve(newNmax,nomb); b_exp_e.resize_virtual(newNmax,nomb);
    b_exp_s.resize_preserve(newNmax,nomb); b_exp_s.resize_virtual(newNmax,nomb);
  }
  index_e.resize_preserve(newNmax);   index_e.resize_virtual(newNmax);
  index_s.resize_preserve(newNmax);   index_s.resize_virtual(newNmax);
  index_s_1.resize_preserve(newNmax); index_s_1.resize_virtual(newNmax);
  _btype_s.resize_preserve(newNmax); _btype_s.resize_virtual(newNmax);
  _btype_e.resize_preserve(newNmax); _btype_e.resize_virtual(newNmax);
  
  for (int i=Nmax; i<newNmax; i++) Empty_e.push_back(i);
  for (int i=Nmax; i<newNmax; i++) Empty_s.push_back(i);
}


void nIntervals::ReleaseSomeTempMemory()
{
  if (!common::Qsvd){
    _exp_e.~function2D();
    _exp_s.~function2D();
    b_exp_e.~function2D();
    b_exp_s.~function2D();
  }
  _time_e.~function1D();
  _time_s.~function1D();
  index_e.~function1D();
  index_s.~function1D();
  _btype_e.~function1D();
  _btype_s.~function1D();
}



inline pair<int,int> nIntervals::InsertExponents(double t_start, int is, int btyp_s, double t_end, int ie, int btyp_e)
{
  int nsize = index_e.size()-Empty_e.size();
  
  int to_insert_c = Empty_e.front();
  Empty_e.pop_front();
  for (int i=nsize; i>ie; i--) index_e[i] = index_e[i-1];
  index_e[ie] = to_insert_c;
  
  int to_insert_cd = Empty_s.front();
  Empty_s.pop_front();
  for (int i=nsize; i>is; i--){
    index_s[i] = index_s[i-1];
    index_s_1[index_s[i]] = i;
  }
  index_s[is] = to_insert_cd;
  index_s_1[index_s[is]] = is;
  
  
  _time_e[to_insert_c] = t_end;
  _time_s[to_insert_cd] = t_start;

  _time_e.resize_virtual(_time_e.size()+1);
  _time_s.resize_virtual(_time_s.size()+1);
  
  _btype_s[to_insert_cd] = btyp_s;
  _btype_e[to_insert_c ] = btyp_e;
  Nbtyp_e[btyp_e]++;
  Nbtyp_s[btyp_s]++;


  if (! common::Qsvd ){
    for (int i=0; i<piom->size(); i++){
      double xe = t_end*(*piom)[i];
      _exp_e[to_insert_c][i] = dcomplex(cos(xe),sin(xe));
      double xs = t_start*(*piom)[i];
      _exp_s[to_insert_cd][i] = dcomplex(cos(xs),sin(xs));
    }
    if (common::SampleSusc){
      for (int i=0; i<biom->size(); i++){
	double xe = t_end*(*biom)[i];
	b_exp_e[to_insert_c][i] = dcomplex(cos(xe),sin(xe));
	double xs = t_start*(*biom)[i];
	b_exp_s[to_insert_cd][i] = dcomplex(cos(xs),sin(xs));
      }
    }
  }
  return make_pair(to_insert_cd, to_insert_c);
}

inline void nIntervals::RemoveExponents(int is, int ie)
{
  int nsize = index_e.size()-Empty_e.size();

  Nbtyp_e[_btype_e[index_e[ie]]]--;
  Nbtyp_s[_btype_s[index_s[is]]]--; 
  
  Empty_e.push_back(index_e[ie]);
  Empty_s.push_back(index_s[is]);
  for (int i=is; i<nsize-1; i++){
    index_s[i] = index_s[i+1];
    index_s_1[index_s[i]]=i;
  }
  for (int i=ie; i<nsize-1; i++) index_e[i] = index_e[i+1];
  
  _time_e.resize_virtual(_time_e.size()-1);
  _time_s.resize_virtual(_time_s.size()-1);
}

inline void nIntervals::MoveExponents(int type, double t_old, int i_old, double t_new, int i_new)
{
  int nsize = index_e.size()-Empty_e.size();
  
  if (type==c){
    int to_insert = index_e[i_old];
    _time_e[to_insert] = t_new;
    for (int i=i_old; i<nsize-1; i++) index_e[i] = index_e[i+1];
    for (int i=nsize-1; i>i_new; i--) index_e[i] = index_e[i-1];
    index_e[i_new] = to_insert;
    
    if (! common::Qsvd){
      for (int i=0; i<piom->size(); i++){
	double xe = t_new*(*piom)[i];
	_exp_e[to_insert][i] = dcomplex(cos(xe),sin(xe));
      }
      if (common::SampleSusc){
	for (int i=0; i<biom->size(); i++){
	  double xe = t_new*(*biom)[i];
	  b_exp_e[to_insert][i].Set(cos(xe),sin(xe));
	}
      }
    }
  }else{
    int to_insert = index_s[i_old];
    _time_s[to_insert] = t_new;
    for (int i=i_old; i<nsize-1; i++){
      index_s[i] = index_s[i+1];
      index_s_1[index_s[i]]=i;
    }
    for (int i=nsize-1; i>i_new; i--){
      index_s[i] = index_s[i-1];
      index_s_1[index_s[i]]=i;
    }
    index_s[i_new] = to_insert;
    index_s_1[index_s[i_new]]=i_new;

    if (!common::Qsvd){
      for (int i=0; i<piom->size(); i++){
	double xe = t_new*(*piom)[i];
	_exp_s[to_insert][i] = dcomplex(cos(xe),sin(xe));
      }
      if (common::SampleSusc){
	for (int i=0; i<biom->size(); i++){
	  double xe = t_new*(*biom)[i];
	  b_exp_s[to_insert][i].Set(cos(xe),sin(xe));
	}
      }
    }
  }
}


inline int nIntervals::FindIndex(double t_new, double t_old, int type)
{
  int is=0;
  if (type==cd){
    while (is<_time_s.size() && time_s(is)<t_new) is++;
  }else{
    while (is<_time_e.size() && time_e(is)<t_new) is++;
  }
  if (t_new>t_old) is--;
  return is;
}

//FindIndex_nearby(t_b, t_a, typea, ia)

inline int nIntervals::FindIndex_nearby(double t_new, int type, int i0)
{
  if (type==cd){
    if (time_s(i0) > t_new){
      int iis = i0-1;
      while (iis>=0 && time_s(iis) > t_new) --iis;
      return iis+1;
    }else{
      int iis = i0+1;
      while (iis<_time_s.size() && t_new > time_s(iis)) ++iis;
      return iis-1;
    }
  }else{
    if (time_e(i0) > t_new){
      int iis = i0-1;
      while (iis>=0 && time_e(iis) > t_new) --iis;
      return iis+1;
    }else{
      int iis = i0+1;
      while (iis<_time_e.size() && t_new > time_e(iis)) ++iis;
      return iis-1;
    }
  }
}

inline void nIntervals::Find_is_ie(double t_start, int& is, double t_end, int& ie)
{
  is=0;
  while (is<_time_s.size() && time_s(is)<t_start) is++;
  ie=0;
  while (ie<_time_e.size() && time_e(ie)<t_end) ie++;
}

  
inline void nIntervals::TestOrder()
{
  int nsize = index_e.size()-Empty_e.size();
  for (int i=0; i<nsize-1; i++)
    if (time_s(i)>time_s(i+1)) cout<<"Times are not ordered! "<<time_s(i)<<" "<<time_s(i+1)<<" "<<i<<" "<<i+1<<endl;
  for (int i=0; i<nsize-1; i++)
    if (time_e(i)>time_e(i+1)) cout<<"Times are not ordered! "<<time_e(i)<<" "<<time_e(i+1)<<" "<<i<<" "<<i+1<<endl;
}

inline void nIntervals::print(ostream& stream)
{
  stream<<setw(4)<<_time_s.size()<<":  s= ";
  for (int is=0; is<_time_s.size(); is++) stream<<setw(4)<<time_s(is)<<" ";
  cout<<" e= ";
  for (int ie=0; ie<_time_e.size(); ie++) stream<<setw(4)<<time_e(ie)<<" ";
  stream<<endl;
}

void nIntervals::copy_data(const nIntervals& s)
{
  //  beta= s.beta;
  _time_e.copy_full(s._time_e);
  _time_s.copy_full(s._time_s);

  if (! common::Qsvd){
    _exp_e = s._exp_e;
    _exp_s = s._exp_s;
    b_exp_e = s.b_exp_e;
    b_exp_s = s.b_exp_s;
  }
  
  index_e = s.index_e;
  index_s = s.index_s;
  index_s_1 = s.index_s_1;
  _btype_e.copy_full(s._btype_e);
  _btype_s.copy_full(s._btype_s);
  Empty_e = s.Empty_e;
  Empty_s = s.Empty_s;
  Nbtyp_s.copy_full(s.Nbtyp_s);
  Nbtyp_e.copy_full(s.Nbtyp_e);
}


void nIntervals::check_index_1()
{
  cout<<"size="<<size()/2<<endl;
  for (int i=0; i<size()/2; i++)
    if (index_s_1[index_s[i]]!=i){
      cout<<"index_1[i] not correct "<<i<<" "<<index_s_1[index_s[i]]<<" index_s_1["<<index_s[i]<<"]="<<index_s_1[index_s[i]]<<endl;
    }
}

void swap_data(functionb<dcomplex>& a_exp, functionb<dcomplex>& b_exp, int len)
{
  dcomplex* __restrict__ fa = a_exp.MemPt();
  dcomplex* __restrict__ fb = b_exp.MemPt();
  dcomplex z;
  for (int i=0; i<len; i++){
    z=fa[i];
    fa[i]=fb[i];
    fb[i]=z;
  }
}
void swapd(double& ta, double& tb)
{
  double tc = ta;
  ta = tb;
  tb = tc;
}

void IntervalsExchangeExponents(int typea, nIntervals& interval_a, int ia, int ia_new,
				int typeb, nIntervals& interval_b, int ib, int ib_new)
{
  int nsizea = interval_a._time_s.size();
  int nsizeb = interval_b._time_s.size();
  function1D<int>      *a_index, *b_index;
  function1D<double>   *a_time,  *b_time; 
  function2D<dcomplex> *a_exp,   *b_exp;  
  function2D<dcomplex> *a_b_exp, *b_b_exp;
  
  if (typea==nIntervals::c){
    a_index = &interval_a.index_e;	
    a_time  = &interval_a._time_e;
    a_exp   = &interval_a._exp_e;
    a_b_exp = &interval_a.b_exp_e;
  }else{
    a_index = &interval_a.index_s;	
    a_time  = &interval_a._time_s;
    a_exp   = &interval_a._exp_s;
    a_b_exp = &interval_a.b_exp_s;
  }
  if (typeb==nIntervals::c){
    b_index = &interval_b.index_e;	
    b_time  = &interval_b._time_e;
    b_exp   = &interval_b._exp_e;
    b_b_exp = &interval_b.b_exp_e;
  }else{
    b_index = &interval_b.index_s;	
    b_time  = &interval_b._time_s;
    b_exp   = &interval_b._exp_s;
    b_b_exp = &interval_b.b_exp_s;
  }

  int a_to_insert = (*a_index)[ia];
  int b_to_insert = (*b_index)[ib];
  
  swapd((*a_time)[a_to_insert],(*b_time)[b_to_insert]);

  if (ia!=ia_new){
    for (int i=ia; i<nsizea-1; i++) (*a_index)[i] = (*a_index)[i+1];
    for (int i=nsizea-1; i>ia_new; i--) (*a_index)[i] = (*a_index)[i-1];
    (*a_index)[ia_new] = a_to_insert;
  }
  if (ib!=ib_new){
    for (int i=ib; i<nsizeb-1; i++) (*b_index)[i] = (*b_index)[i+1];
    for (int i=nsizeb-1; i>ib_new; i--) (*b_index)[i] = (*b_index)[i-1];
    (*b_index)[ib_new] = b_to_insert;
  }

  if (! common::Qsvd){
    swap_data((*a_exp)[a_to_insert], (*b_exp)[b_to_insert], interval_a.piom->size());
    if (common::SampleSusc)
      swap_data((*a_b_exp)[a_to_insert], (*b_b_exp)[b_to_insert], interval_a.biom->size());
  }
  
  if (typea==nIntervals::cd && ia!=ia_new){
    for (int i=min(ia,ia_new); i<nsizea; i++) interval_a.index_s_1[interval_a.index_s[i]]=i;
  }
  if (typeb==nIntervals::cd && ib!=ib_new){
    for (int i=min(ib,ib_new); i<nsizeb; i++) interval_b.index_s_1[interval_b.index_s[i]]=i;
  }
}
