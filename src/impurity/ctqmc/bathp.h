// @Copyright 2007 Kristjan Haule
// 

class BathProb{
public:
  int N_ifl;
  vector<double> BathP;
public:
  BathProb(const vector<double>& BathProbability, const ClusterData& cluster): N_ifl(cluster.N_ifl), BathP(N_ifl){
    for (int ifl=0; ifl<cluster.N_ifl; ifl++){
      for (int i1=0; i1<cluster.ifl_dim[ifl]; i1++){
	for (int i2=0; i2<cluster.ifl_dim[ifl]; i2++){
	  int ib = cluster.tfl_index[ifl][i1][i2];
	  int fl = cluster.bfl_index[ifl][ib];
	  BathP[ifl] = BathProbability[fl];
	}
      }
    }
    if (common::fastFilesystem) Print(cout);
    double dnorm=0;
    for (int ifl=0; ifl<cluster.N_ifl; ifl++) dnorm += BathP[ifl];
    for (int ifl=0; ifl<cluster.N_ifl; ifl++) BathP[ifl]/=dnorm;
    double dsum=0;
    for (int ifl=0; ifl<cluster.N_ifl; ifl++){
      dsum += BathP[ifl];
      BathP[ifl] = dsum;
    }
    if (common::fastFilesystem) Print(cout);
  }
  void Print(ostream& out){
    out<<"BathProbability=";
    for (int ifl=0; ifl<N_ifl; ifl++) out<<BathP[ifl]<<" ";
    out<<endl;
  }
  int ifl(double x){
    for (int i=0; i<BathP.size(); i++)
      if (x<BathP[i]) return i;
    return N_ifl;
  }
};

