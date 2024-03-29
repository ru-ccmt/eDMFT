// @Copyright 2007 Kristjan Haule
// 
class ClusterData;

class NState{
public:
  int istate;
  function2D<double> M;
  double exponent;
public:
  //NState() : M(common::max_size,common::max_size), exponent(0.0) {};
  NState(int size1=1, int size2=1) : M(size1,size2), exponent(0.0) {};
  NState(const NState& C);
  void SetIstate(int istate_){istate=istate_;}
  void Evolve(const NState& C0, const ClusterData& cluster, const funProxy<double>& exp_);
  void apply(const functionb<function2D<double> >& FM, const functionb<int>& Fi, const NState& C);
  bool empty() const;
  void Print(ostream& out) const;
  Number TraceProject(const NState& s);
  Number Project_to_Pra(const NState& s, functionb<Number>& Proj);
  Number ScalarProduct(const NState& s);
  void SetPraState(int i, const ClusterData& cluster);
  int state_n() const {return istate;}
  void SetEqual(const NState& m);
  string TotalSize() const;
  string CurrentSize() const;
  static Timer t_evolve, t_apply;
  NState& operator=(const NState& C);
  friend std::ostream& operator<< (std::ostream& stream, const NState& st);
  void ReleaseMemory()
  { M.ReleaseMemory(); }
};

