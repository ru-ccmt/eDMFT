// @Copyright 2007 Kristjan Haule
// 
#ifdef _MPI
#include <mpi.h>
using namespace std;

void Reduce(int my_rank, int Master, int mpi_size, function1D<double>& histogram, function2D<dcomplex>& Gd, function2D<dcomplex>& Sd,
	    function2D<double>& AverageProbability, double& asign, double& asign_fast, function1D<double>& nlc, function1D<double>& kaver, function2D<double>& susc,
	    function2D<dcomplex>& asuscg, int DOsize, function2D<double>& Gtau, function5D<dcomplex>& VertexH, function5D<dcomplex>& VertexF, function1D<int>& Gd_deg,
	    function2D<double>& AP_transition, bool cmp_vertex, bool QHB2, bool SampleSusc, bool SampleTransitionP)
{
  function2D<double> cAverageProbability;
  function1D<double> chistogram;
  function1D<double> cnlc;
  function1D<double> ckaver;
  function2D<dcomplex> cGd;
  function2D<dcomplex> cSd;
  function2D<double> cSusc;
  function2D<dcomplex> casuscg;
  function2D<double> cGtau;
  function1D<int> cGd_deg;
  function2D<double> cAP_transition;
  function1D<double> casign(2), asign_(2);
  if (my_rank==Master){
    cnlc.resize(nlc.size());
    cAverageProbability.resize(AverageProbability.fullsize_N(),AverageProbability.fullsize_Nd());
    cGd.resize(Gd.fullsize_N(),Gd.fullsize_Nd());
    if (QHB2) cSd.resize(Sd.fullsize_N(),Sd.fullsize_Nd());
    cGtau.resize(Gtau.fullsize_N(),Gtau.fullsize_Nd());
    ckaver.resize(kaver.size());
    if (SampleSusc) cSusc.resize(susc.fullsize_N(), susc.fullsize_Nd());
    if (SampleSusc && DOsize>0) casuscg.resize(asuscg.fullsize_N(),asuscg.fullsize_Nd());
    cGd_deg.resize(Gd_deg.size());
    if (SampleTransitionP) cAP_transition.resize(AP_transition.fullsize_N(),AP_transition.fullsize_Nd());
  }
  asign_[0]=asign; asign_[1]=asign_fast;
  
  // It turns out that different processors work with different size histograms. We need to sum them up,
  // but first we need to make all histogram of equal size
  int global_Nmax;
  int histogram_Nmax = histogram.size();
  
  //MPI::COMM_WORLD.Allreduce(&histogram_Nmax, &global_Nmax, 1, MPI_INT, MPI_MAX);
  MPI_Allreduce(&histogram_Nmax, &global_Nmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    
  function1D<double> histogram_copy(global_Nmax);
  histogram_copy=0.0;
  for (int i=0; i<histogram.size(); i++) histogram_copy[i]=histogram[i];
  if (my_rank==Master) chistogram.resize(global_Nmax);
  
  //MPI::COMM_WORLD.Reduce(histogram_copy.MemPt(), chistogram.MemPt(), global_Nmax, MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(histogram_copy.MemPt(), chistogram.MemPt(), global_Nmax, MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);

  //MPI::COMM_WORLD.Reduce(AverageProbability.MemPt(), cAverageProbability.MemPt(), AverageProbability.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(AverageProbability.MemPt(), cAverageProbability.MemPt(), AverageProbability.fullsize2(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);

  //MPI::COMM_WORLD.Reduce(asign_.MemPt(), casign.MemPt(), 2, MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(asign_.MemPt(), casign.MemPt(), 2, MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  
  //MPI::COMM_WORLD.Reduce(nlc.MemPt(), cnlc.MemPt(), nlc.size(), MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(nlc.MemPt(), cnlc.MemPt(), nlc.size(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  
  //MPI::COMM_WORLD.Reduce(kaver.MemPt(), ckaver.MemPt(), kaver.size(), MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(kaver.MemPt(), ckaver.MemPt(), kaver.size(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
	     
  //MPI::COMM_WORLD.Reduce(Gtau.MemPt(), cGtau.MemPt(), Gtau.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(Gtau.MemPt(), cGtau.MemPt(), Gtau.fullsize2(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
	     
  //MPI::COMM_WORLD.Reduce(Gd.MemPt(), cGd.MemPt(), Gd.fullsize2()*2, MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(Gd.MemPt(), cGd.MemPt(), Gd.fullsize2()*2, MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);

  if (QHB2)
    //MPI::COMM_WORLD.Reduce(Sd.MemPt(), cSd.MemPt(), Sd.fullsize2()*2, MPI_DOUBLE, MPI_SUM, Master);
    MPI_Reduce(Sd.MemPt(), cSd.MemPt(), Sd.fullsize2()*2, MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  
  if (SampleSusc)
    //MPI::COMM_WORLD.Reduce(susc.MemPt(), cSusc.MemPt(), susc.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);
    MPI_Reduce(susc.MemPt(), cSusc.MemPt(), susc.fullsize2(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  if (SampleSusc && DOsize>0)
    MPI_Reduce(asuscg.MemPt(), casuscg.MemPt(), asuscg.fullsize2()*2, MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  
  //MPI::COMM_WORLD.Reduce(Gd_deg.MemPt(), cGd_deg.MemPt(), Gd_deg.size(), MPI_INT, MPI_SUM, Master);
  MPI_Reduce(Gd_deg.MemPt(), cGd_deg.MemPt(), Gd_deg.size(), MPI_INT, MPI_SUM, Master, MPI_COMM_WORLD);
  
  if (SampleTransitionP)
    //MPI::COMM_WORLD.Reduce(AP_transition.MemPt(), cAP_transition.MemPt(), AP_transition.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);
    MPI_Reduce(AP_transition.MemPt(), cAP_transition.MemPt(), AP_transition.fullsize2(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
    
  if (cmp_vertex){
    function2D<dcomplex> cVertex(VertexH.N3, VertexH.N4);
    int psize = VertexH.N3*VertexH.N4;
    for (int i0=0; i0<VertexH.N0; i0++){
      for (int i1=0; i1<VertexH.N1; i1++){
	for (int i2=0; i2<VertexH.N2; i2++){
	  
	  cVertex=0.0;
	  dcomplex* f = &VertexH(i0,i1,i2,0,0);
	  //MPI::COMM_WORLD.Reduce(f, cVertex.MemPt(), psize*2, MPI_DOUBLE, MPI_SUM, Master);
	  MPI_Reduce(f, cVertex.MemPt(), psize*2, MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
	  
	  //	  for (int i3=0; i3<VertexH.N1; i3++)   !!! WAS THIS A BUG FOR MANY MANY YEARS????
	  //	    for (int i4=0; i4<VertexH.N2; i4++)
	  for (int i3=0; i3<VertexH.N3; i3++) 
	    for (int i4=0; i4<VertexH.N4; i4++)
	      VertexH(i0,i1,i2,i3,i4) = cVertex(i3,i4)*(1./mpi_size);
	  
	  cVertex=0.0;
	  f = &VertexF(i0,i1,i2,0,0);
	  //MPI::COMM_WORLD.Reduce(f, cVertex.MemPt(), psize*2, MPI_DOUBLE, MPI_SUM, Master);
	  MPI_Reduce(f, cVertex.MemPt(), psize*2, MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
	  
	  for (int i3=0; i3<VertexH.N3; i3++)
	    for (int i4=0; i4<VertexH.N4; i4++)
	      VertexF(i0,i1,i2,i3,i4) = cVertex(i3,i4)*(1./mpi_size);
	}
      }
    }
  }

  if (my_rank==Master){
    histogram.resize(global_Nmax);
    histogram = chistogram;
    histogram *= (1./mpi_size);
    AverageProbability = cAverageProbability;
    asign = casign[0];
    asign_fast = casign[1];
    AverageProbability *= (1./mpi_size);
    nlc = cnlc;
    nlc *= (1./mpi_size);
    kaver = ckaver;
    kaver *= (1./mpi_size);
    Gd = cGd;
    if (QHB2) Sd = cSd;
    Gtau = cGtau;
    Gtau *= (1./mpi_size);
    if (SampleSusc){
      susc = cSusc;
      susc *= (1./mpi_size);
    }
    if (SampleSusc && DOsize>0){
      asuscg = casuscg;
      asuscg *= (1./mpi_size);
    }
    asign *= (1./mpi_size);
    asign_fast *= (1./mpi_size);
    Gd_deg = cGd_deg;
    if (SampleTransitionP){
      AP_transition = cAP_transition;
      AP_transition *=  (1./mpi_size);
    }
  }
}

void ReduceS(int my_rank, int Master, int mpi_size, function1D<double>& histogram, function2D<double>& Gd, function2D<double>& Ft,
	     function2D<double>& AverageProbability, double& asign, double& asign_fast, function1D<double>& nlc, function1D<double>& kaver,
	     function2D<double>& Gtau, function5D<double>& VH, function1D<int>& Gd_deg,
	     function2D<double>& AP_transition, vector<function2D<double> >& Gsvd,
	     bool cmp_vertex, bool QHB2, bool SampleTransitionP, ostream& clog)
{
  function1D<double> chistogram;
  function2D<double> cGd;
  function2D<double> cFt;
  function2D<double> cAverageProbability;
  function1D<double> cnlc;
  function1D<double> ckaver;
  function2D<double> cGtau;
  function1D<int> cGd_deg;
  function2D<double> cAP_transition;
  vector<function2D<double> > cGsvd(Gsvd.size());
  function1D<double> casign(2), asign_(2);
  if (my_rank==Master){
    cGd.resize(Gd.fullsize_N(),Gd.fullsize_Nd());
    if (QHB2) cFt.resize(Ft.fullsize_N(),Ft.fullsize_Nd());
    cAverageProbability.resize(AverageProbability.fullsize_N(),AverageProbability.fullsize_Nd());
    cnlc.resize(nlc.size());
    ckaver.resize(kaver.size());
    if (Gtau.size_Nd()) cGtau.resize(Gtau.fullsize_N(),Gtau.fullsize_Nd());
    cGd_deg.resize(Gd_deg.size());
    if (SampleTransitionP) cAP_transition.resize(AP_transition.fullsize_N(),AP_transition.fullsize_Nd());
    for (int ifl=0; ifl<Gsvd.size(); ifl++) cGsvd[ifl].resize( Gsvd[ifl].fullsize_N(), Gsvd[ifl].fullsize_Nd() );
  }else{
    for (int ifl=0; ifl<Gsvd.size(); ifl++) cGsvd[ifl].resize( 1, 1 );
  }
  asign_[0]=asign; asign_[1]=asign_fast;

  // It turns out that different processors work with different size histograms. We need to sum them up,
  // but first we need to make all histogram of equal size
  int global_Nmax;
  int histogram_Nmax = histogram.size();
  //MPI::COMM_WORLD.Allreduce(&histogram_Nmax, &global_Nmax, 1, MPI_INT, MPI_MAX);
  MPI_Allreduce(&histogram_Nmax, &global_Nmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  function1D<double> histogram_copy(global_Nmax);
  histogram_copy=0.0;
  for (int i=0; i<histogram.size(); i++) histogram_copy[i]=histogram[i];
  
  if (my_rank==Master) chistogram.resize(global_Nmax);
  
  //MPI::COMM_WORLD.Reduce(histogram_copy.MemPt(), chistogram.MemPt(), global_Nmax, MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(histogram_copy.MemPt(), chistogram.MemPt(), global_Nmax, MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  
  //MPI::COMM_WORLD.Reduce(AverageProbability.MemPt(), cAverageProbability.MemPt(), AverageProbability.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(AverageProbability.MemPt(), cAverageProbability.MemPt(), AverageProbability.fullsize2(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);

  //MPI::COMM_WORLD.Reduce(asign_.MemPt(), casign.MemPt(), 2, MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(asign_.MemPt(), casign.MemPt(), 2, MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  
  //MPI::COMM_WORLD.Reduce(nlc.MemPt(), cnlc.MemPt(), nlc.size(), MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(nlc.MemPt(), cnlc.MemPt(), nlc.size(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  
  //MPI::COMM_WORLD.Reduce(kaver.MemPt(), ckaver.MemPt(), kaver.size(), MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(kaver.MemPt(), ckaver.MemPt(), kaver.size(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  
  if (Gtau.size_Nd())
    //MPI::COMM_WORLD.Reduce(Gtau.MemPt(), cGtau.MemPt(), Gtau.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);
    MPI_Reduce(Gtau.MemPt(), cGtau.MemPt(), Gtau.fullsize2(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);

  //MPI::COMM_WORLD.Reduce(Gd.MemPt(), cGd.MemPt(), Gd.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);
  MPI_Reduce(Gd.MemPt(), cGd.MemPt(), Gd.fullsize2(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);

  if (QHB2)
    //MPI::COMM_WORLD.Reduce(Ft.MemPt(), cFt.MemPt(), Ft.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);
    MPI_Reduce(Ft.MemPt(), cFt.MemPt(), Ft.fullsize2(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  
  //MPI::COMM_WORLD.Reduce(Gd_deg.MemPt(), cGd_deg.MemPt(), Gd_deg.size(), MPI_INT, MPI_SUM, Master);
  MPI_Reduce(Gd_deg.MemPt(), cGd_deg.MemPt(), Gd_deg.size(), MPI_INT, MPI_SUM, Master, MPI_COMM_WORLD);
  
  if (SampleTransitionP)
    //MPI::COMM_WORLD.Reduce(AP_transition.MemPt(), cAP_transition.MemPt(), AP_transition.fullsize2(), MPI_DOUBLE, MPI_SUM, Master);
    MPI_Reduce(AP_transition.MemPt(), cAP_transition.MemPt(), AP_transition.fullsize2(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
  
  if (cmp_vertex){
    if (my_rank==Master) cout<<"Reducing svd.VH"<<endl;

    function3D<double> cVertex(VH.N2, VH.N3, VH.N4);
    int psize = VH.N2*VH.N3*VH.N4;
    
    for (int i0=0; i0<VH.N0; i0++){
      for (int i1=0; i1<VH.N1; i1++){
	
	double* f = &VH(i0,i1,0,0,0);
	double* cf = &cVertex(0,0,0);
	//MPI::COMM_WORLD.Reduce(f, cf, psize, MPI_DOUBLE, MPI_SUM, Master);
	MPI_Reduce(f, cf, psize, MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);

	if (my_rank==Master){
	  for (int i2=0; i2<VH.N2; i2++){
	    for (int i3=0; i3<VH.N3; i3++)
	      for (int i4=0; i4<VH.N4; i4++)
		VH(i0,i1,i2,i3,i4) = cVertex(i2,i3,i4)*(1./mpi_size);
	  }
	}
      }
    }
  }

  for (int ifl=0; ifl<Gsvd.size(); ifl++){
    MPI_Reduce(Gsvd[ifl].MemPt(), cGsvd[ifl].MemPt(), Gsvd[ifl].fullsize2(), MPI_DOUBLE, MPI_SUM, Master, MPI_COMM_WORLD);
    
    if (my_rank==Master){
      for (int i1=0; i1<Gsvd[ifl].size_N(); i1++)
	for (int i2=0; i2<Gsvd[ifl].size_Nd(); i2++)
	  Gsvd[ifl](i1,i2) = cGsvd[ifl](i1,i2)*(1./mpi_size);
    }
  }

    
  if (my_rank==Master){
    histogram.resize(global_Nmax);
    histogram = chistogram;
    histogram *= (1./mpi_size);
    AverageProbability = cAverageProbability;
    AverageProbability *= (1./mpi_size);
    asign = casign[0];
    asign_fast = casign[1];
    nlc = cnlc;
    nlc *= (1./mpi_size);
    kaver = ckaver;
    kaver *= (1./mpi_size);
    if (Gtau.size_Nd()){
      Gtau = cGtau;
      Gtau *= (1./mpi_size);
    }
    asign *= (1./mpi_size);
    asign_fast *= (1./mpi_size);
    Gd = cGd;
    if (QHB2) Ft = cFt;
    Gd_deg = cGd_deg;
    if (SampleTransitionP){
      AP_transition = cAP_transition;
      AP_transition *=  (1./mpi_size);
    }
  }
}

void MPI_Init(int argc, char* argv[], int& my_rank, int& mpi_size, int& Master)
{
  /* // Unfortunately C++ binding is deprecated...
  MPI::Init(argc, argv);
  my_rank = MPI::COMM_WORLD.Get_rank();
  mpi_size = MPI::COMM_WORLD.Get_size();
  */
  ::MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  Master = 0;
}

void MPI_finalize()
{
  //MPI::Finalize();
  MPI_Finalize();
}

#else

using namespace std;

void Reduce(int my_rank, int Master, int mpi_size, function1D<double>& histogram, function2D<dcomplex>& Gd, function2D<dcomplex>& Sd,
	    function2D<double>& AverageProbability, double& asign, double& asign_fast, function1D<double>& nlc, function1D<double>& kaver, function2D<double>& susc,
	    function2D<dcomplex>& asuscg, int DOsize, function2D<double>& Gtau, function5D<dcomplex>& VertexH, function5D<dcomplex>& VertexF, function1D<int>& Gd_deg,
	    function2D<double>& AP_transition, bool cmp_vertex, bool QHB2, bool SampleSusc, bool SampleTransitionP){}


void ReduceS(int my_rank, int Master, int mpi_size, function1D<double>& histogram, function2D<double>& Gd, function2D<double>& Ft,
	     function2D<double>& AverageProbability, double& asign, double& asign_fast, function1D<double>& nlc, function1D<double>& kaver,
	     function2D<double>& Gtau, function5D<double>& VH, function1D<int>& Gd_deg,
	     function2D<double>& AP_transition, vector<function2D<double> >& Gsvd, bool cmp_vertex, bool QHB2, bool SampleTransitionP, ostream& clog){};

void MPI_Init(int argc, char* argv[], int& my_rank, int& mpi_size, int& Master)
{
  my_rank = 0;
  mpi_size = 1;
  Master = 0;
}

void MPI_finalize(){}
#endif
