// @Copyright 2007 Kristjan Haule
// 
#include <map>
#include <cmath>

class common{
public:
  static int N_flavors;
  static int N_ifl;
  static double beta;
  static double U;
  static double mu;
  static int max_size;
  static double minF; // If applying F makes the state smaller than that, we should assume it is zero.
  static int tsample;
  static int tsampleFast;
  static int warmup;
  static int CleanUpdate;
  static double PChangeOrder;
  static double PMove;
  static int Ncout;
  static long long Naver;
  static int my_rank;
  static int mpi_size;
  static double TwoKinks;
  static int GlobalFlip;
  static double treshold;
  static int SampleGtau;
  static int PreciseP;
  static double minDeltat;
  static bool SampleSusc;
  static bool cmp_vertex;
  static int SampleVertex;
  static double maxNoise;
  static int fastFilesystem;
  static bool LazyTrace;
  static bool QHB2;
  static double smallest_dt;
  static int Segment;
  static bool SampleTransitionP;
  static bool SaveStatus;
  static bool Qsvd;
  static double PlocalMoves;
  static void SetParameters(int rank, int size, double mu_, double U_, double beta_, int max_size_, int N_flavors_,
			    int N_ifl_, int tsample_, int tsampleFast_,
			    int warmup_, int CleanUpdate_, double minF_, double PChangeOrder_, double PMove_, double PlocalMoves_,
			    int Ncout_, long long Naver_, double TwoKinks_, int GlobalFlip_, double treshold_, int SampleGtau_,
			    int PreciseP_, double minDeltat_, bool SampleSusc_, int SampleVertex_, double maxNoise_,
			    bool LazyTrace_, int Segment_, int fastFilesystem_, bool SampleTransitionP_, bool Qsvd_)
  {
    my_rank = rank;
    mpi_size = size;
    mu = mu_; U = U_; beta=beta_;
    max_size = max_size_;
    N_flavors = N_flavors_;
    N_ifl = N_ifl_;
    tsample = tsample_;
    tsampleFast = tsampleFast_;
    warmup = warmup_;
    CleanUpdate = CleanUpdate_;
    minF = minF_;
    PChangeOrder = PChangeOrder_;
    PMove = PMove_;
    Ncout = Ncout_;
    Naver = Naver_;
    TwoKinks = TwoKinks_;
    GlobalFlip = GlobalFlip_;
    treshold = treshold_;
    SampleGtau = SampleGtau_;
    PreciseP = PreciseP_;
    minDeltat = minDeltat_;
    SampleSusc = SampleSusc_;
    cmp_vertex = (SampleVertex_>0);
    SampleVertex = SampleVertex_;
    maxNoise = maxNoise_;
    LazyTrace = LazyTrace_;
    fastFilesystem=fastFilesystem_;
    smallest_dt=3e-16*beta;
    Segment = Segment_;
    SampleTransitionP = SampleTransitionP_;
    PlocalMoves = PlocalMoves_;
    Qsvd = Qsvd_;
  }
};

inline double linear(double x, const std::pair<double,double>& x1, const std::pair<double,double>& x2)
{
  return x2.second + (x1.second-x2.second)*(x-x2.first)/(x1.first-x2.first);
}
inline double quadr(double x, const std::pair<double,double>& x1, const std::pair<double,double>& x2, const std::pair<double,double>& x3)
{
  return (x-x2.first)*(x-x3.first)*x1.second/((x1.first-x2.first)*(x1.first-x3.first))+
    (x-x1.first)*(x-x3.first)*x2.second/((x2.first-x1.first)*(x2.first-x3.first))+
    (x-x1.first)*(x-x2.first)*x3.second/((x3.first-x1.first)*(x3.first-x2.first));
}
inline double expin(double x, const std::pair<double,double>& x1, const std::pair<double,double>& x2)
{
  if (x2.second*x1.second<0) return x1.second + (x2.second-x1.second)*(x-x1.first)/(x2.first-x1.first);
  double beta = log(x2.second/x1.second)/(x1.first-x2.first);
  double alpha = x1.second*exp(beta*x1.first);
  return alpha*exp(-beta*x);
}
