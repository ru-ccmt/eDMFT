#include <blitz/array.h>
using namespace std;
using namespace blitz;

#ifdef _MPI
#include <mpi.h>
#include <iostream>

// This is a stream which does not output anything
class NullStream : public std::ostream
{
    // streambuffer doing nothing
    class NullBuffer : public std::streambuf
    {
    public:
        int overflow(int c) noexcept override { return c; }
    } nullBuffer;

public:
#pragma warning(suppress: 26455)
    NullStream() : std::ostream(&nullBuffer) {}
    NullStream(const NullStream&) = delete;
    NullStream& operator=(const NullStream&) = delete;
};

class MPIcls{
public:
  int rank, size, master;
  ostream* out;
  Array<int,1> sendcounts;
  Array<int,1> displacements;
  //
  MPIcls(int argc, char* argv[])
  {
    ::MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    master = 0;
    if (rank==master){
      out = new ofstream("ctqmc.log");
    }else{
      out = new NullStream;
    }
    (*out)<<"Parallel start of bold_ctqmc on process " << rank << " of " << size << std::endl;
  }
  ~MPIcls()
  {
    delete out; out=NULL;
    MPI_Finalize();
  }
  void Abort(){MPI_Abort(MPI_COMM_WORLD, 1);}
  void Print(const string& message){
    if (rank==master) (*out)<<message<<endl;
  }
  bool I_am_root(){ return rank==master;}
  void Pbin_Pnorm_sum_allreduce(Array<double,2>& Pbin_Pnorm_sum)
  {
    if (size>1){
      MPI_Allreduce(MPI_IN_PLACE,Pbin_Pnorm_sum.data(),Pbin_Pnorm_sum.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      Pbin_Pnorm_sum *= 1.0/size;
    }
  }
  void bcast_dSml(Array<double,3>& dSml, Array<double,1>& V0norm, Array<double,1>& Pnorm, Array<double,1>& Phi_accumulate,
		  Array<double,1>& Fake_real_overlap, Array<double,2>& Pbin_Pnorm_sum, 
		  Array<double,1>& N_measure, Array<double,1>& aver_sign)
  {
    int Norders = dSml.extent(0);
    int Npseudo = dSml.extent(1);
    int lmax    = dSml.extent(2);
    MPI_Bcast(dSml.data(),         dSml.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(V0norm.data(),           Norders, MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(Pnorm.data(),            Norders, MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(Phi_accumulate.data(),   Norders, MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(Fake_real_overlap.data(),Norders, MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(Pbin_Pnorm_sum.data(), 2*Norders, MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(N_measure.data(),        Norders, MPI_DOUBLE, master, MPI_COMM_WORLD); 
    MPI_Bcast(aver_sign.data(),        Norders, MPI_DOUBLE, master, MPI_COMM_WORLD); 
  }
  void bcast_dSm_again(Array<double,1>& Phi_tot, Array<double,3>& Stval, bool cmpPhysG, Array<complex<double>,3>& dG00, Array<complex<double>,3>& dF00R){
    int Norders = Stval.extent(0);
    int Npseudo = Stval.extent(1);
    int Nbin    = Stval.extent(2);
    MPI_Bcast(Stval.data(),  Stval.size(),   MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(Phi_tot.data(),Phi_tot.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
    if (cmpPhysG){
      MPI_Bcast(dG00.data(), dG00.size(),  MPI_C_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
      MPI_Bcast(dF00R.data(),dF00R.size(), MPI_C_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
      //MPI_Bcast(dF00L.data(),dF00L.size(), MPI_C_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
      //MPI_Bcast(dI00.data(),dI00.size(), MPI_C_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
    }
  }
  void bcast_dSm_again(Array<double,1>& Phi_tot, Array<double,3>& Stval, bool cmpPhysG, Array<double,3>& dG00, Array<double,3>& dF00R){
    int Norders = Stval.extent(0);
    int Npseudo = Stval.extent(1);
    int Nbin    = Stval.extent(2);
    MPI_Bcast(Stval.data(),   Stval.size(),  MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(Phi_tot.data(),Phi_tot.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
    if (cmpPhysG){
      MPI_Bcast(dG00.data(),dG00.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
      MPI_Bcast(dF00R.data(),dF00R.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
      //MPI_Bcast(dF00L.data(),dF00L.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
      //MPI_Bcast(dI00.data(),dI00.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
    }
  }
  void reduce_dSml(Array<double,3>& dSml, Array<double,1>& Pnorm, Array<double,1>& Pnorm2, Array<double,1>& Phi_accumulate, Array<double,1>& Fake_real_overlap,
		   Array<double,2>& Pbin_Pnorm_sum, //Array<double,1>& Pnorm_sum, Array<double,1>& Pbin_sum,
		   Array<double,1>& N_measure, Array<double,1>& aver_sign)
  {
    int Norders = dSml.extent(0);
    int Npseudo = dSml.extent(1);
    int lmax    = dSml.extent(2);
    if (rank==master){
      MPI_Reduce(MPI_IN_PLACE, dSml.data(), Norders*Npseudo*lmax, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      //
      MPI_Reduce(MPI_IN_PLACE, Pnorm.data(),             Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, Pnorm2.data(),            Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, Phi_accumulate.data(),    Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, Fake_real_overlap.data(), Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, Pbin_Pnorm_sum.data(),  2*Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, N_measure.data(),         Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, aver_sign.data(),         Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    }else{
      MPI_Reduce(dSml.data(), dSml.data(), Norders*Npseudo*lmax, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      //
      MPI_Reduce(Pnorm.data(),             Pnorm.data(),             Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(Pnorm2.data(),            Pnorm2.data(),            Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(Phi_accumulate.data(),    Phi_accumulate.data(),    Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(Fake_real_overlap.data(), Fake_real_overlap.data(), Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(Pbin_Pnorm_sum.data(),    Pbin_Pnorm_sum.data(),  2*Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(N_measure.data(),         N_measure.data(),         Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(aver_sign.data(),         aver_sign.data(),         Norders, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    }
    if (rank==master){
      double norm = 1./size;
      dSml   *= norm;
      Pnorm  *= norm;
      Pnorm2 *= norm;
      Phi_accumulate *= norm;
      Fake_real_overlap *= norm;
      Pbin_Pnorm_sum *= norm;
      N_measure *= norm;
      aver_sign *= norm;
    }
  }
  void reduce_dG00(Array<complex<double>,3>& dG00l,Array<complex<double>,3>& dF00l){
    if (rank==master){
      MPI_Reduce(MPI_IN_PLACE, dG00l.data(), dG00l.size(), MPI_C_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, dF00l.data(), dF00l.size(), MPI_C_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD);
      double norm = 1./size;
      dG00l *= 1./size;
      dF00l *= 1./size;
    }else{
      MPI_Reduce(dG00l.data(), dG00l.data(), dG00l.size(), MPI_C_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(dF00l.data(), dF00l.data(), dF00l.size(), MPI_C_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD);
    }
  }
  void reduce_dG00(Array<double,3>& dG00l,Array<double,3>& dF00l){
    if (rank==master){
      MPI_Reduce(MPI_IN_PLACE, dG00l.data(), dG00l.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, dF00l.data(), dF00l.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      double norm = 1./size;
      dG00l *= 1./size;
      dF00l *= 1./size;
    }else{
      MPI_Reduce(dG00l.data(), dG00l.data(), dG00l.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(dF00l.data(), dF00l.data(), dF00l.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    }
  }
  void bcast_one_boolean(bool& read_succeeded)
  {
    MPI_Bcast(&read_succeeded, 1, MPI_C_BOOL, master, MPI_COMM_WORLD);
  }
  void bcast_two_boolean(bool& bool1, bool& bool2)
  {
    bool buff[2]={bool1,bool2};
    MPI_Bcast(buff, 2, MPI_C_BOOL, master, MPI_COMM_WORLD);
    bool1 = buff[0];
    bool2 = buff[1];
  }
  void bcast_ints(int* array, int n)
  {
    MPI_Bcast(array, n, MPI_INT, master, MPI_COMM_WORLD);
  }
  void bcast_double(double* Sm_data, int size){
    MPI_Bcast(Sm_data, size, MPI_DOUBLE, master, MPI_COMM_WORLD);
  }
  void bcast_complex(complex<double>* Sm_data, int size){
    MPI_Bcast(Sm_data, size, MPI_C_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
  }
  void BcastDeltaw(Array<double,1>& iom, Array<complex<double>,2>& Deltaw, int Nbath, const string& fDelta)
  {
    int size[2]={static_cast<int>(iom.size()),static_cast<int>(Deltaw.size())};
    MPI_Bcast(&size, 2, MPI_INT, master, MPI_COMM_WORLD);
    if (size[0]==0){
      cerr<<"file not found: "<<fDelta<<endl;
      (*out)<<"file not found: "<<fDelta<<endl;
      this->~MPIcls();
      exit(1);
    }
    if (rank!=master){
      iom.resize(size[0]);
      Deltaw.resize(Nbath,size[0]);
    }
    MPI_Bcast(iom.data(), size[0], MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(Deltaw.data(), size[1], MPI_C_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
  }
  void bcast_Sm(Array<double,1>& xn, Array<complex<double>,2>& Sm, int Npseudo)
  {
    int w_size = xn.size();
    MPI_Bcast(&w_size, 1, MPI_INT, master, MPI_COMM_WORLD);
    //cout<<"Rank="<< rank <<" Bcast Sm data of size="<<w_size<<endl;
    if (rank!=master){
      xn.resize(w_size);
      Sm.resize(Npseudo,w_size);
    }
    bcast_complex(Sm.data(), Sm.size());
    bcast_double(xn.data(), xn.size());
    //cout<<"Rank="<<rank<<" Bcast Sm finished"<<endl;
  }
  void bcast_Sm(Array<double,1>& xn, Array<double,2>& Sm, int Npseudo)
  {
    int w_size = xn.size();
    MPI_Bcast(&w_size, 1, MPI_INT, master, MPI_COMM_WORLD);
    //cout<<"Rank="<< rank <<" Bcast Sm data of size="<<w_size<<endl;
    if (rank!=master){
      xn.resize(w_size);
      Sm.resize(Npseudo,w_size);
    }
    bcast_double(Sm.data(), Sm.size());
    bcast_double(xn.data(), xn.size());
    //cout<<"Rank="<<rank<<" Bcast Sm finished"<<endl;
  }
  template<int Nrank>
  void reduce_mweight(Array<double,Nrank>& gt, Array<double,Nrank>& t_hist){
    if (rank==master){
      MPI_Reduce(MPI_IN_PLACE,      gt.data(),     gt.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE,  t_hist.data(), t_hist.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    }else{
      MPI_Reduce(gt.data(),         gt.data(),     gt.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
      MPI_Reduce(t_hist.data(), t_hist.data(), t_hist.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    }
    gt     *= (1./size);
    t_hist *= (1./size);
  };
  void bcast_Vol_diag(Array<double,1>& V0norm, Array<double,1>& Vol_diagrams0, Array<double,1>& Vol_diagrams1)
  {
    MPI_Bcast(V0norm.data(), V0norm.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(Vol_diagrams0.data(), Vol_diagrams0.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(Vol_diagrams1.data(), Vol_diagrams1.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
  }
  void bcast_diags(Array<int,3>& diags,Array<int,2>& pseudos, int stored_size, int norder)
  {
    if (rank!=master){
      diags.resize(stored_size,norder,2);
      pseudos.resize(stored_size,2*norder);
    }
    MPI_Bcast(diags.data(),   stored_size*norder*2, MPI_INT, master, MPI_COMM_WORLD);
    MPI_Bcast(pseudos.data(), stored_size*norder*2, MPI_INT, master, MPI_COMM_WORLD);
  }
  void bcast_mweight2(Array<double,2>& gt, Array<double,2>& t_hist, int Nsym, int Nbin, bool& self_consistent, bool& read_histog_succ) const{
    MPI_Bcast(&read_histog_succ, 1, MPI_C_BOOL, master, MPI_COMM_WORLD);
    if (read_histog_succ){
      self_consistent=true;
      if (rank!=master) gt.resize(Nsym,Nbin);
      MPI_Bcast(gt.data(), Nsym*Nbin, MPI_DOUBLE, master, MPI_COMM_WORLD);
      MPI_Bcast(t_hist.data(), Nsym*Nbin, MPI_DOUBLE, master, MPI_COMM_WORLD);
    }
  }
  void bcast_mweight(Array<double,1>& gt, Array<double,1>& t_hist, int Nbin, bool& self_consistent, bool& read_histog_succ) const
  {
    MPI_Bcast(&read_histog_succ, 1, MPI_C_BOOL, master, MPI_COMM_WORLD);
    if (read_histog_succ){
      self_consistent=true;
      if (rank!=master) gt.resize(Nbin);
      MPI_Bcast(gt.data(), Nbin, MPI_DOUBLE, master, MPI_COMM_WORLD);
      MPI_Bcast(t_hist.data(), Nbin, MPI_DOUBLE, master, MPI_COMM_WORLD);
    }
  }
  auto SplitArray(int irank, int size, int leng)
  {
    struct retVals{
      int iqs,iqe;
    };
    int pr_proc = (leng % size==0) ? int(leng/size) : int(leng/size+1);
    retVals r;
    if (size<=leng){
      r.iqs = min(irank*pr_proc,leng);
      r.iqe = min((irank+1)*pr_proc,leng);
    }else{
      int rstep=(size+1)/leng;
      if (irank%rstep==0 and irank/rstep<leng){
	r.iqs = irank/rstep;
	r.iqe = r.iqs+1;
      }else{
	if (irank/rstep<leng){
	  r.iqs = irank/rstep;
	  r.iqe = irank/rstep;
	 }else{
	  r.iqs = leng-1;
	  r.iqe = leng-1;
	}
      }
    }
    return r;
  }
  auto mpiSplitArray(int leng)
  {
    struct retVals{
      int iqs,iqe;
    };
    sendcounts.resize(size);
    displacements.resize(size);
    for (int irank=0; irank<size; ++irank){
      auto [iqs,iqe] = SplitArray(irank,size,leng);
      sendcounts(irank) = iqe-iqs;
      displacements(irank) = iqs;
    }
    //cout<<"rank="<<rank<<" sendcounts="<<sendcounts<<" displacement="<<displacements<<endl;
    return SplitArray(rank,size,leng);
  }
  void AllGatherOCA(int Npseudo_unique, int Nbath_unique, int iqs, int iqe, Array<double,2>& TSmt, Array<complex<double>,2>& mGlc, Array<complex<double>,2>& mFlc){
    sendcounts    *= Npseudo_unique;
    displacements *= Npseudo_unique;
    MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Npseudo_unique,MPI_DOUBLE,TSmt.data(),sendcounts.data(),displacements.data(),MPI_DOUBLE,MPI_COMM_WORLD);
  
    sendcounts    /= Npseudo_unique;
    displacements /= Npseudo_unique;
    sendcounts    *= Nbath_unique;
    displacements *= Nbath_unique;
    MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,mGlc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,mFlc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
  }
  void AllGather1OCA(int Npseudo_unique, int iqs, int iqe, Array<double,2>& TSmt){
    sendcounts    *= Npseudo_unique;
    displacements *= Npseudo_unique;
    MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Npseudo_unique,MPI_DOUBLE,TSmt.data(),sendcounts.data(),displacements.data(),MPI_DOUBLE,MPI_COMM_WORLD);
  }
  void AllGather2OCA(int Nbath_unique, int iqs, int iqe, Array<complex<double>,2>& Glc, Array<complex<double>,2>& mGlc, Array<complex<double>,2>& Flc, Array<complex<double>,2>& mFlc){
    sendcounts    *= Nbath_unique;
    displacements *= Nbath_unique;
    MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,Glc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,Flc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,mGlc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,mFlc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
  }
  void AllGather3OCA(int Nbath_unique, int iqs, int iqe, Array<double,2>& Glc, Array<double,2>& Flc){
    sendcounts    *= Nbath_unique;
    displacements *= Nbath_unique;
    MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Nbath_unique,MPI_DOUBLE,Glc.data(),sendcounts.data(),displacements.data(),MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Nbath_unique,MPI_DOUBLE,Flc.data(),sendcounts.data(),displacements.data(),MPI_DOUBLE,MPI_COMM_WORLD);
  }
  void AllGatherNCA1(int Npseudo_unique, int iqs, int iqe, Array<double,2>& TSmt){
    if (size>1){
      sendcounts *= Npseudo_unique;
      displacements *= Npseudo_unique;
      MPI_Allgatherv(MPI_IN_PLACE,(iqe-iqs)*Npseudo_unique,MPI_DOUBLE,TSmt.data(),sendcounts.data(),displacements.data(),MPI_DOUBLE,MPI_COMM_WORLD);
    }
  }
  void AllGatherNCA2(int Nbath_unique, int iws, int iwe, Array<complex<double>,2>& TGc, Array<complex<double>,2>& TFc, Array<complex<double>,2>& TGcm, Array<complex<double>,2>& TFcm)
  {
    if (size>1){
      sendcounts *= Nbath_unique;
      displacements *= Nbath_unique;
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TGc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TGcm.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFcm.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    }
  }
  void AllGatherNCA3(int Nbath_unique, int iws, int iwe, Array<complex<double>,2>& TGc, Array<complex<double>,2>& TGcm,
		     Array<complex<double>,2>& TFc, Array<complex<double>,2>& TFcm, Array<complex<double>,2>& TFc2, Array<complex<double>,2>& TFcm2)
  {
    if (size>1){
      sendcounts *= Nbath_unique;
      displacements *= Nbath_unique;
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TGc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TGcm.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFcm.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFc2.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFcm2.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    }
  }
  void AllGatherNCA4(int Nbath_unique, int iws, int iwe, Array<complex<double>,2>& TGc, Array<complex<double>,2>& TGcm,
		     Array<complex<double>,2>& TFc, Array<complex<double>,2>& TFcm, Array<complex<double>,2>& TFc2, Array<complex<double>,2>& TFcm2,
		     Array<complex<double>,2>& TFc3, Array<complex<double>,2>& TFcm3, Array<complex<double>,2>& TFc4, Array<complex<double>,2>& TFcm4)
  {
    if (size>1){
      sendcounts *= Nbath_unique;
      displacements *= Nbath_unique;
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TGc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TGcm.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFc.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFcm.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFc2.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFcm2.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFc3.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFcm3.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFc4.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE,(iwe-iws)*Nbath_unique,MPI_C_DOUBLE_COMPLEX,TFcm4.data(),sendcounts.data(),displacements.data(),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    }
  }  
  void Bcast_New_mesh(Array<double,1>& om, Array<double,1>& few, Array<std::complex<double>,2>& Sm, Array<std::complex<double>,2>& Gm)
  {
    int dims[] = {0,0,0};
    if (rank == master) {
      dims[0] = static_cast<int>(om.size());
      dims[1] = static_cast<int>(Sm.extent(0));
      dims[2] = static_cast<int>(Gm.extent(0));
    }
    MPI_Bcast(dims,   3, MPI_INT, master, MPI_COMM_WORLD);
    // 2) Resize on receivers
    if (rank != master) {
        om.resize(dims[0]);
        few.resize(dims[0]);
        Sm.resize(dims[1], dims[0]);
        Gm.resize(dims[2], dims[0]);
    }
    /*
    if (rank == master) {
      std::fprintf(stderr, "[root] om=%d Sm=(%d,%d) Gm=(%d,%d)\n",
		   dims[0], dims[1], dims[0], dims[2], dims[0]);
    } else {
      std::fprintf(stderr, "[rank %d pre] om=%d Sm=(%d,%d) Gm=(%d,%d)\n",
		   rank, (int)om.size(), (int)Sm.extent(0), (int)Sm.extent(1),
		   (int)Gm.extent(0), (int)Gm.extent(1));
    }
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);
    */
    // Optional: assert contiguity if your array type can be non-contiguous
    assert(om.isStorageContiguous() && few.isStorageContiguous() &&  Sm.isStorageContiguous() && Gm.isStorageContiguous());
    // 3) Broadcast raw data
    // Counts are element counts, not bytes.
    const int cnt_om  = static_cast<int>(om.size());
    const int cnt_Sm  = static_cast<int>(Sm.size()); // total elements = rows*cols
    const int cnt_Gm  = static_cast<int>(Gm.size());

    MPI_Bcast(om.data(),  cnt_om, MPI_DOUBLE,            master, MPI_COMM_WORLD);
    MPI_Bcast(few.data(), cnt_om, MPI_DOUBLE,            master, MPI_COMM_WORLD);
    MPI_Bcast(reinterpret_cast<void*>(Sm.data()), cnt_Sm, MPI_CXX_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
    MPI_Bcast(reinterpret_cast<void*>(Gm.data()), cnt_Gm, MPI_CXX_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
  }

  void Bcast_New_mesh(Array<double,1>& om, Array<std::complex<double>,2>& Sm, Array<std::complex<double>,2>& Gm, Array<std::complex<double>,2>& Dm)
  {
    int dims[] = {0,0,0,0,0,0,0};
    if (rank == master) {
      dims[0] = static_cast<int>(om.size());
      dims[1] = static_cast<int>(Sm.extent(0));
      dims[2] = static_cast<int>(Sm.extent(1));
      dims[3] = static_cast<int>(Gm.extent(0));
      dims[4] = static_cast<int>(Gm.extent(1));
      dims[5] = static_cast<int>(Dm.extent(0));
      dims[6] = static_cast<int>(Dm.extent(1));
    }
    MPI_Bcast(dims,   7, MPI_INT, master, MPI_COMM_WORLD);
    // 2) Resize on receivers
    if (rank != master) {
      om.resize(dims[0]);
      Sm.resize(dims[1], dims[2]);
      Gm.resize(dims[3], dims[4]);
      Dm.resize(dims[5], dims[6]);
    }
    // Optional: assert contiguity if your array type can be non-contiguous
    assert(om.isStorageContiguous() && Sm.isStorageContiguous() &&  Gm.isStorageContiguous() && Dm.isStorageContiguous());
    // 3) Broadcast raw data
    MPI_Bcast(om.data(), om.size(), MPI_DOUBLE,             master, MPI_COMM_WORLD);
    MPI_Bcast(Sm.data(), Sm.size(), MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
    MPI_Bcast(Gm.data(), Gm.size(), MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
    MPI_Bcast(Dm.data(), Dm.size(), MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
  }
  
  void AllReduce_data(double* data, int size)
  {
    MPI_Allreduce(MPI_IN_PLACE,data,size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }
  void allgather_nca(double& Phi_tot, double* Sm_data, int Nt_Npseudo, double* G00_data, double* F00_data, int Nt_Nbath)
  {
    MPI_Allreduce(MPI_IN_PLACE,&Phi_tot,1,         MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,Sm_data, Nt_Npseudo,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,G00_data,Nt_Nbath,  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,F00_data,Nt_Nbath,  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }

  template<int N>
  void Reduce_sum_array(Array<double,N>& arr)
  {
    if (rank==master){
      MPI_Reduce(MPI_IN_PLACE, arr.data(), arr.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    }else{
      MPI_Reduce(arr.data(), arr.data(), arr.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    }
  }
  template<int N>
  void Reduce_sum_array(Array<complex<double>,N>& arr)
  {
    if (rank==master){
      MPI_Reduce(MPI_IN_PLACE, arr.data(), arr.size(), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD);
    }else{
      MPI_Reduce(arr.data(), arr.data(), arr.size(), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD);
    }
  }
  
  void Reduce_sum_number(double& value)
  {
    if (rank==master){
      MPI_Reduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    }else{
      MPI_Reduce(&value, &value, 1, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    }
  }
  void Reduce_sum_number(complex<double>& value)
  {
    if (rank==master){
      MPI_Reduce(MPI_IN_PLACE, &value, 1, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD);
    }else{
      MPI_Reduce(&value, &value, 1, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD);
    }
  }
  
  template<int N>
  void Reduce_sum_array(Array<double,N>& res, const Array<double,N>& arr)
  {
    if (rank==master) res.resize(arr.shape());
    MPI_Reduce(arr.data(), res.data(), arr.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
  }
  template<int N>
  void Bcast_simple_array(Array<double,N>& arr)
  {
    MPI_Bcast(arr.data(),  arr.size(),  MPI_DOUBLE, master, MPI_COMM_WORLD);
  }
  template<int N>
  void Bcast_simple_array(Array<complex<double>,N>& arr)
  {
    MPI_Bcast(arr.data(),  arr.size(),  MPI_CXX_DOUBLE_COMPLEX, master, MPI_COMM_WORLD);
  }
  void bcast_string(string& message)
  {
    int length = message.size();
    MPI_Bcast(&length, 1, MPI_INT, master, MPI_COMM_WORLD);
    if (rank != master)  message.resize(length);
    MPI_Bcast(&message[0], length, MPI_CHAR, master, MPI_COMM_WORLD);
  }
  template<int N>
  void AllReduceSum_simple_array(Array<double,N>& arr)
  {
    MPI_Allreduce(MPI_IN_PLACE, arr.data(), arr.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  template<int N>
  void AllReduceSum_simple_array(Array<complex<double>,N>& arr)
  {
    MPI_Allreduce(MPI_IN_PLACE, arr.data(), arr.size(), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  }
  void AllReduceSum_rawdata(double* arr, int N)
  {
    MPI_Allreduce(MPI_IN_PLACE, arr, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  void Assign_jobs(std::vector<std::vector<BlockJob>>& assign, std::vector<BlockJob>& jobs)
  {
    std::vector<double> load(size, 0.0);
    // sorting with most expensive first
    std::sort(jobs.begin(), jobs.end(), [](auto& a, auto& b){ return a.cost > b.cost; });
    for (auto& job : jobs) {
      int target = std::min_element(load.begin(), load.end()) - load.begin();
      assign[target].push_back(job);
      load[target] += job.cost;
    }
    for (int icore=0; icore<size; ++icore)
      std::sort(assign[icore].begin(),assign[icore].end(), [](auto& a, auto& b){ return a.id < b.id; });
  }
  std::pair<double,int> FindSmallest_lambda0(std::pair<double,int> local)
  {
    // local.first  = λ₀ on this rank
    // local.second = followPeak (any int payload you want to keep)
    struct { double val; int idx; } in { local.first, local.second };
    struct { double val; int idx; } out;
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    return { out.val, out.idx };
  }

  void WhichIndicesComputed(Array<int,1>& which_ii,
			    const std::vector<std::vector<BlockJob>>& assign,
			    const vector<Array<int,1>>& cix_Block,
			    const vector<Array<int,2>>& cix_Sgm_symmetry,
			    const vector<deque<TinyVector<int,3>>>& cix_Sgm_indx, int cix_Npseudo_unique)
  {
    deque<int> which_ii_;
    for (int ij=0; ij<assign[rank].size(); ij++){
      int ib = assign[rank][ij].id;
      int N = cix_Block[ib].size(); // block size
      for (int i=0; i<N; ++i){
	for (int j=0; j<N; ++j){
	  int ii = abs(cix_Sgm_symmetry[ib](i,j))-1;
	  if (ii>=0 and  bl::all(cix_Sgm_indx[ii][0] == TinyVector<int,3>(ib,i,j))){
	    // here we ignore the case that cix_Sgm_indx[ii][j](0) could be -ib. But the first entry in cix_Sgm_indx should always be positive, hence we can ignore that case.
	    which_ii_.push_back(ii);
	  }
	}
      }
    }
    // sorting is to make sure that when we have single core, no extra reordering is needed of the data.
    sort(which_ii_.begin(),which_ii_.end());
    which_ii.resize(which_ii_.size());
    for (int i=0; i<which_ii_.size(); ++i) which_ii(i) = which_ii_[i];
  }

  void Gather_Gm(Array<double,2>& Gm,
		 const Array<int,1>& which_ii, int cix_Npseudo_unique, int tau_size)
  {
    if (size>1){
      int  count_rows = which_ii.size();
      Array<int,1> recv_counts_rows(1), recv_counts_dbl(1);
      if (rank==master) {
	recv_counts_rows.resize(size);
	recv_counts_dbl.resize(size);
      }
      MPI_Gather(&count_rows, 1, MPI_INT, recv_counts_rows.data(), 1, MPI_INT, master, MPI_COMM_WORLD); // root receives counts from everyone
      if (rank==master){
	recv_counts_dbl = recv_counts_rows*tau_size;
	//clog<<"recv_counts_rows="<<recv_counts_rows<<" recv_counts_dbl="<<recv_counts_dbl<<endl;
      }
      Array<int,1> disp_rows(size+1), disp_dbl(size+1), which_all(1);
      int Ntotal;
      if (rank==master){
	disp_rows(0)=0;
	for (int i=1; i<=size; ++i) disp_rows(i) = disp_rows(i-1) + recv_counts_rows(i-1);
	disp_dbl(Range::all()) = disp_rows(Range::all())*tau_size;
	Ntotal = disp_rows(disp_rows.extent(0)-1);
	which_all.resize(Ntotal);
	//clog<<"Ntotal="<<Ntotal<<" Npseudo="<<cix_Npseudo_unique<<" disp_rows="<<disp_rows<<" disp_dbl="<<disp_dbl<<endl;
      }
      MPI_Gatherv(which_ii.data(), count_rows, MPI_INT, which_all.data(), recv_counts_rows.data(), disp_rows.data(), MPI_INT, master, MPI_COMM_WORLD);
      Array<double,2> Gm_combined;
      if (rank==master){
	//clog<<"which_all="<<which_all<<endl;
	Gm_combined.resize(Ntotal,tau_size);    // total doubles
      }
      MPI_Gatherv(Gm.data(),  count_rows*tau_size, MPI_DOUBLE, Gm_combined.data(), recv_counts_dbl.data(), disp_dbl.data(), MPI_DOUBLE, master, MPI_COMM_WORLD);
      Gm.resize(cix_Npseudo_unique, tau_size);
      if (rank == master){
	for (int i=0; i<Ntotal; ++i)
	  Gm(which_all(i),Range::all()) = Gm_combined(i,Range::all());
      }
      MPI_Bcast(Gm.data(), Gm.size(), MPI_DOUBLE, master, MPI_COMM_WORLD);
    }
  }
  
  void Gather_Prob(Array<double,1>& Prob_eig)
  {
    if (rank == master)
      MPI_Reduce(MPI_IN_PLACE, Prob_eig.data(), Prob_eig.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    else
      MPI_Reduce(Prob_eig.data(), nullptr, Prob_eig.size(), MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
  }
  
};

#else
class MPIcls{
public:
  int rank, master, size;
  ostream* out;
  MPIcls(int argc, char* argv[])
  {
    rank = 0;
    size = 1;
    master = 0;
    out = &cout;
    std::cout<<"Not parallel!"<<std::endl;
  }
  ~MPIcls(){}

  bool I_am_root(){ return rank==master;}
  
  void Bcast(double& U, double& mu, double& beta, int Ms[], int& max_order, double& max_diff, double& max_diff_over_error, int& Ntau, int& Nitt, int& tsample, double& nom_beta, double& Q0, double& StartLambda, double& dLambda, double& mix, int& Optimize_tau_grid, int& Nt_small){}

  void allgather(double* g_at_beta, int Nf){}
  
  void allgather2(complex<double>* selfs, int Ns, int Nf){}
  
  void allgather3(double* Gm, int Nf, int m){}

  void allgather4(double* Gc, double* Fc, double& Phi_tot, int Nt, int m){}
  void allgather_nca(double& Phi_tot, double* Sm_data, int Nt_Npseudo, double* G00_data, double* F00_data, int Nt_Nbath){}

  //void reduceall(function2D<double>& dSm, function2D<double>& dG00, function2D<double>& dF00, double& aPhi, double& errPhi, double& succ, double& asign, double& error){}

  void bcast2(bool& should_break){}
  void BcastDeltaw(Array<double,1>& iom, Array<complex<double>,2>& Deltaw, int Nbath, const string& fDelta){}
  void bcast_ints(int* array, int n){}
  void bcast_double(double* Sm_data, int size){}
  void bcast_complex(complex<double>* Sm_data, int size){}
  void bcast_Sm(Array<double,1>& xn, Array<complex<double>,2>& Sm, int Npseudo){};
  void bcast_Sm(Array<double,1>& xn, Array<double,2>& Sm, int Npseudo){};
  //void bcast_ret(int& return_value){}
  //void bcast_stored_size(int& stored_size){};
  void bcast_diags(Array<int,3>& diags,Array<int,2>& pseudos, int stored_size, int norder){};
  void bcast_one_boolean(bool& read_succeeded){}
  void bcast_two_boolean(bool& bool1, bool& bool2){}
  void bcast_dSml(Array<double,3>& dSml, Array<double,1>& V0norm, Array<double,1>& Pnorm, Array<double,1>& Phi_accumulate,
		  Array<double,1>& Fake_real_overlap, Array<double,2>& Pbin_Pnorm_sum, 
		  Array<double,1>& N_measure, Array<double,1>& aver_sign){}
  void Print(const string& message){
    (*out)<<message<<endl;
  }
  void Abort(){};
  void bcast_mweight2(Array<double,2>& gt, Array<double,2>& t_hist, int Nsym, int Nbin, bool& self_consistent, bool& read_histog_succ) const{};
  void bcast_mweight(Array<double,1>& gt, Array<double,1>& t_hist, int Nbin, bool& self_consistent, bool& read_histog_succ) const{};
  void bcast_Vol_diag(Array<double,1>& V0norm, Array<double,1>& Vol_diagrams0, Array<double,1>& Vol_diagrams1){};
  template<class T>void bcast_dSm_again(Array<double,1>& Phi_tot, Array<double,3>& dSm, bool cmpPhysG, Array<T,3>& dG00, Array<T,3>& dF00R/*, Array<T,3>& dF00L, Array<T,3>& dI00*/){};
  template<int Nrank> void reduce_mweight(Array<double,Nrank>& gt, Array<double,Nrank>& t_hist){};
  void Pbin_Pnorm_sum_allreduce(Array<double,2>& Pbin_Pnorm_sum){};
  //void reduce_dG00R(Array<complex<double>,3>& dG00l,Array<complex<double>,3>& dF00l){};
  void reduce_dSml(Array<double,3>& dSml, Array<double,1>& Pnorm, Array<double,1>& Pnorm2, Array<double,1>& Phi_accumulate, Array<double,1>& Fake_real_overlap,
		   Array<double,2>& Pbin_Pnorm_sum, Array<double,1>& N_measure, Array<double,1>& aver_sign){};
  void Bcast_New_mesh(Array<double,1>& om, Array<double,1>& few, Array<complex<double>,2>& Sm, Array<complex<double>,2>& Gm){};
  void Bcast_New_mesh(Array<double,1>& om, Array<std::complex<double>,2>& Sm, Array<std::complex<double>,2>& Gm, Array<std::complex<double>,2>& Dm){};
  void AllGatherNCA1(int Npseudo_unique, int iqs, int iqe, Array<double,2>& TSmt){};
  void AllGatherNCA2(int Nbath_unique, int iws, int iwe, Array<complex<double>,2>& TGc, Array<complex<double>,2>& TFc, Array<complex<double>,2>& TGcm, Array<complex<double>,2>& TFcm){};
  void AllGatherNCA3(int Nbath_unique, int iws, int iwe, Array<complex<double>,2>& TGc, Array<complex<double>,2>& TGcm,
		     Array<complex<double>,2>& TFc, Array<complex<double>,2>& TFcm, Array<complex<double>,2>& TFc2, Array<complex<double>,2>& TFcm2){};
  void AllGatherOCA(int Npseudo_unique, int Nbath_unique, int iqs, int iqe, Array<double,2>& TSmt, Array<complex<double>,2>& mGlc, Array<complex<double>,2>& mFlc){};
  void AllGather1OCA(int Npseudo_unique, int iqs, int iqe, Array<double,2>& TSmt){};
  void AllGather2OCA(int Nbath_unique, int iqs, int iqe, Array<complex<double>,2>& Glc, Array<complex<double>,2>& Flc, Array<complex<double>,2>& mGlc, Array<complex<double>,2>& mFlc){};
  void AllGather3OCA(int Nbath_unique, int iqs, int iqe, Array<double,2>& Glc, Array<double,2>& Flc){};
  void AllReduce_data(double* data, int size){};
  template<class T>
  void reduce_dG00(Array<T,3>& dG00l,Array<T,3>& dF00l){};

  auto mpiSplitArray(int leng){
    struct retVals{
      int iqs,iqe;
      retVals(int _iqs_, int _iqe_):iqs(_iqs_),iqe(_iqe_){}
    };
    return retVals(0,leng);
  }

  template<int N>
  void Reduce_sum_array(Array<double,N>& arr){};
  template<int N>
  void Reduce_sum_array(Array<complex<double>,N>& arr){};  
  template<int N>
  void Reduce_sum_array(Array<double,N>& res, const Array<double,N>& arr){};
  template<int N>
  void Bcast_simple_array(Array<double,N>& arr){};
  template<int N>
  void Bcast_simple_array(Array<complex<double>,N>& arr){};
  void Reduce_sum_number(double& value){};
  void Reduce_sum_number(complex<double>& value){};
  void bcast_string(string& message){};
  template<int N>
  void AllReduceSum_simple_array(Array<double,N>& arr){};
  template<int N>
  void AllReduceSum_simple_array(Array<complex<double>,N>& arr){};
  void AllReduceSum_rawdata(double* arr, int N){};

  void Assign_jobs(std::vector<std::vector<BlockJob>>& assign, std::vector<BlockJob>& jobs) { assign[0]=jobs; }
  std::pair<double,int> FindSmallest_lambda0(std::pair<double,int> local) { return local; }

  void WhichIndicesComputed(Array<int,1>& which_ii,
			    const std::vector<std::vector<BlockJob>>& assign, const vector<Array<int,1>>& cix_Block, const vector<Array<int,2>>& cix_Sgm_symmetry, const vector<deque<TinyVector<int,3>>>& cix_Sgm_indx, int cix_Npseudo_unique)
  {
    which_ii.resize(cix_Npseudo_unique);
    for (int i=0; i<cix_Npseudo_unique; ++i) which_ii(i)=i;
  }
  
  void Gather_Gm(Array<double,2>& Gm, const Array<int,1>& which_ii, int cix_Npseudo_unique, int tau_size){}
  void Gather_Prob(Array<double,1>& Prob_eig){};

};

#endif
