// @Copyright 2007 Kristjan Haule
// 
#ifdef _MPI
#include <mpi.h>

int ClusterData::BcastClusterData(int my_rank, int Master)
{// To avoid reading CIX file by all processors in massively parallel run
 // we broadcast data using MPI calls. Firts we create two large arrays
 // and combine all integer variables into one array, and all double variables
 // into the second array. We then broadcast the two arrays to all processor 
 // in a few MPI calls only.
  int isize, dsize, ssize;
  if (my_rank==Master){
    // CommonStructure:
    isize = 8;
    // ifl_dim
    isize += ifl_dim.size();
    // fl_from_ifl
    isize += fl_from_ifl.size_N()*fl_from_ifl.size_Nd();
    // tfl_index
    for (int i=0; i<tfl_index.size(); i++) isize += tfl_index[i].size_N()*tfl_index[i].size_Nd();
    // vfl_index
    for (int i=0; i<vfl_index.size(); i++) isize += vfl_index[i].size_N()*vfl_index[i].size_Nd();
    // v2fl_index
    for (int i=0; i<v2fl_index.size(); i++) isize += v2fl_index[i].size();
    // sign
    for (int i=0; i<sign.size(); i++) isize += sign[i].size();
    // conjg
    for (int i=0; i<conjg.size(); i++) isize += conjg[i].size();
    // bfl_index
    for (int i=0; i<bfl_index.size(); i++) isize += bfl_index[i].size();
    // gflip_index
    isize += gflip_index.size();
    // Ns
    isize += Ns.size();
    // Ks
    isize += Ks.size();
    // msize
    isize += msize_.size();
    // F_i
    isize += F_i.size_N()*F_i.size_Nd();
    // F_M(,).sizes()
    isize += F_M.size_N()*F_M.size_Nd()*2;
    // Id.size
    isize += 1;
    // Id.first
    isize += Id.size();
    // ssize0
    isize += 1; 
    // Enes.size()
    isize += 1;
    // DOsize
    isize += 1;
    if (DOsize>0){
      // DF_i
      isize += DF_i.size_N()*DF_i.size_Nd();
      // DF_inv
      isize += DF_inv.size_N()*DF_inv.size_Nd();
      // DF_M(,).sizes()
      isize += DF_M.size_N()*DF_M.size_Nd()*2;
    }
    if (common::QHB2){
      // UC
      isize += 4*UC.size() + 1;
    }
    
    // double variables
    // Sz
    dsize  = Sz.size();
    // epsk
    dsize += epsk.size();
    // Ene
    //dsize += Ene.size_N()*Ene.size_Nd();
    // Enes
    dsize += Enes.size();
    // Ene0
    dsize += Ene0.size_N()*Ene0.size_Nd();
    // Spin
    dsize += Spin.size_N()*Spin.size_Nd();
    // F_M

    //common::gout<<"1)dsize before F_M="<<dsize<<endl;
    
    for (int i=0; i<F_M.size_N(); i++)
      for (int j=0; j<F_M.size_Nd(); j++)
	dsize += F_M(i,j).size_N() * F_M(i,j).size_Nd();


    //common::gout<<"1)dsize before HF_M="<<dsize<<endl;
    
    for (int i=0; i<HF_M.size(); i++)
      for (int j=0; j<HF_M[i].size_N(); j++)
	for (int k=0; k<HF_M[i].size_Nd(); k++)
	  dsize += HF_M[i](j,k).size_N()*HF_M[i](j,k).size_Nd();

    //common::gout<<"1)dsize after HF_M="<<dsize<<endl;
    
    // Id.second
    dsize += Id.size();

    // DF_M
    for (int i=0; i<DF_M.size_N(); i++)
      for (int j=0; j<DF_M.size_Nd(); j++)
	dsize += DF_M(i,j).size_N() * DF_M(i,j).size_Nd();
    
    // Uc
    if (common::QHB2){
      dsize += UC.size();
    }
    
    ssize = RealSigma.size()+1+ImagSigma.size()+1;
  }
  
  //MPI::COMM_WORLD.Bcast(&isize, 1, MPI::INT, Master);
  MPI_Bcast(&isize, 1, MPI_INT, Master, MPI_COMM_WORLD);
  //MPI::COMM_WORLD.Bcast(&dsize, 1, MPI::INT, Master);
  MPI_Bcast(&dsize, 1, MPI_INT, Master, MPI_COMM_WORLD);
  //MPI::COMM_WORLD.Bcast(&ssize, 1, MPI::INT, Master);  
  MPI_Bcast(&ssize, 1, MPI_INT, Master, MPI_COMM_WORLD);  
  //MPI::COMM_WORLD.Bcast(&QHB1, 1, MPI::BOOL, Master);
  int HB1 =  static_cast<int>(QHB1); // need to convert to int, as C does not have bool type
  MPI_Bcast(&HB1, 1, MPI_INT, Master, MPI_COMM_WORLD);
  QHB1 = static_cast<bool>(HB1);
  //MPI::COMM_WORLD.Bcast(&common::QHB2, 1, MPI::BOOL, Master);
  int HB2 = static_cast<int>(common::QHB2); // need to convert to int, as C does not have bool type
  MPI_Bcast(&HB2, 1, MPI_INT, Master, MPI_COMM_WORLD);
  common::QHB2 = static_cast<bool>(HB2);
  
  int* ibuff = new int[isize];
  double* dbuff = new double[dsize];
  char* sbuff = new char[ssize];
  
  if (my_rank==Master){
    
    int max_dim = fl_from_ifl.size_Nd();
    
    ibuff[0] = Nm;
    ibuff[1] = nsize;
    ibuff[2] = N_ifl;
    ibuff[3] = N_flavors;
    ibuff[4] = Nvfl;
    ibuff[5] = max_size;
    ibuff[6] = Osize;
    ibuff[7] = max_dim;
    int ip=8-1;
    // ifl_dim
    for (int i=0; i<ifl_dim.size(); i++) ibuff[++ip] = ifl_dim[i];
    // fl_from_ifl
    for (int i=0; i<fl_from_ifl.size_N(); i++){                   
      for (int j=0; j<fl_from_ifl.size_Nd(); j++){
	ibuff[++ip] = fl_from_ifl(i,j);
      }
    }
    // tfl_index
    for (int i=0; i<tfl_index.size(); i++)
      for (int j=0; j<tfl_index[i].size_N(); j++)
	for (int k=0; k<tfl_index[i].size_Nd(); k++)
	  ibuff[++ip] = tfl_index[i](j,k);
    // vfl_index
    for (int i=0; i<vfl_index.size(); i++)
      for (int j=0; j<vfl_index[i].size_N(); j++)
	for (int k=0; k<vfl_index[i].size_Nd(); k++)
	  ibuff[++ip] = vfl_index[i](j,k);
    // v2fl_index
    for (int i=0; i<v2fl_index.size(); i++)
      for (int j=0; j<v2fl_index[i].size(); j++)
	ibuff[++ip] = v2fl_index[i][j];
    // sign
    for (int i=0; i<sign.size(); i++)
      for (int j=0; j<sign[i].size(); j++)
	ibuff[++ip] = sign[i][j];
    // conjg
    for (int i=0; i<conjg.size(); i++)
      for (int j=0; j<conjg[i].size(); j++)
	ibuff[++ip] = conjg[i][j];
    // bfl_index
    for (int i=0; i<bfl_index.size(); i++)
      for (int j=0; j<bfl_index[i].size(); j++)
	ibuff[++ip] = bfl_index[i][j];
    // gflip_index
    for (int i=0; i<gflip_index.size(); i++) ibuff[++ip] = gflip_index[i];
    // Ns
    for (int i=0; i<Ns.size(); i++) ibuff[++ip] = Ns[i];
    // Ks
    for (int i=0; i<Ks.size(); i++) ibuff[++ip] = Ks[i];
    // msize
    for (int i=0; i<msize_.size(); i++) ibuff[++ip] = msize_[i];
    // F_i
    for (int i=0; i<F_i.size_N(); i++)
      for (int j=0; j<F_i.size_Nd(); j++)
	ibuff[++ip] = F_i(i,j);
    // F_M(,).sizes()
    for (int i=0; i<F_M.size_N(); i++){
      for (int j=0; j<F_M.size_Nd(); j++){
	//common::gout<<my_rank<<") F_M("<<i<<","<<j<<").size="<<F_M(i,j).size_N()<<" "<<F_M(i,j).size_Nd()<<endl;
	ibuff[++ip] = F_M(i,j).size_N();
	ibuff[++ip] = F_M(i,j).size_Nd();
      }
    }
    // Id.size
    ibuff[++ip] = Id.size();
    // Id.first
    for (map<int,double>::const_iterator i=Id.begin(); i!=Id.end(); i++) ibuff[++ip] = i->first;
    // RealSigma
    ibuff[++ip] = RealSigma.size();
    // Enes
    ibuff[++ip] = Enes.size();
    // DOsize
    ibuff[++ip] = DOsize;
    
    if (DOsize>0){
      // DF_i
      for (int i=0; i<DF_i.size_N(); i++)
	for (int j=0; j<DF_i.size_Nd(); j++)
	  ibuff[++ip] = DF_i(i,j);
      // DF_inv
      for (int i=0; i<DF_inv.size_N(); i++)
	for (int j=0; j<DF_inv.size_Nd(); j++)
	  ibuff[++ip] = DF_inv(i,j);
      // DF_M(,).sizes()
      for (int i=0; i<DF_M.size_N(); i++){
	for (int j=0; j<DF_M.size_Nd(); j++){
	  ibuff[++ip] = DF_M(i,j).size_N();
	  ibuff[++ip] = DF_M(i,j).size_Nd();
	}
      }
    }
    
    if (common::QHB2){
      ibuff[++ip] = UC.size();
      for (int k=0; k<UC.size(); k++){
	ibuff[++ip] = UC[k].a;
	ibuff[++ip] = UC[k].b;
	ibuff[++ip] = UC[k].c;
	ibuff[++ip] = UC[k].d;
      }
    }
    
    //common::gout<<"Id.first=";
    //for (map<int,double>::const_iterator i=Id.begin(); i!=Id.end(); i++) common::gout<<i->first<<" ";
    //common::gout<<endl;
    //common::gout<<"ssize0="<<RealSigma.size();
    
    // double variables
    int dp=-1;
    // Sz
    for (int i=0; i<Sz.size(); i++) dbuff[++dp] = Sz[i];
    // epsk
    for (int i=0; i<epsk.size(); i++) dbuff[++dp] = epsk[i];
    // Ene
    /*
    for (int i=0; i<Ene.size_N(); i++)
      for (int j=0; j<Ene.size_Nd(); j++)
	dbuff[++dp] = Ene(i,j);
    */
    for (int i=0; i<Enes.size(); i++)
      dbuff[++dp] = Enes[i];
    // Ene0
    for (int i=0; i<Ene0.size_N(); i++)
      for (int j=0; j<Ene0.size_Nd(); j++)
	dbuff[++dp] = Ene0(i,j);
    // Spin
    for (int i=0; i<Spin.size_N(); i++)
      for (int j=0; j<Spin.size_Nd(); j++)
	dbuff[++dp] = Spin(i,j);

    //common::gout<<"2)dsize before F_M="<<dp+1<<endl;
    
    // F_M
    for (int i=0; i<F_M.size_N(); i++){
      for (int j=0; j<F_M.size_Nd(); j++){
	for (int k=0; k<F_M(i,j).size_N(); k++){
	  for (int l=0; l<F_M(i,j).size_Nd(); l++){
	    dbuff[++dp] = F_M(i,j)(k,l);
	  }
	}
      }
    }
    //common::gout<<"2)dsize before HF_M="<<dp+1<<endl;
    
    // HF_M
    for (int i=0; i<HF_M.size(); i++)
      for (int j=0; j<HF_M[i].size_N(); j++)
	for (int k=0; k<HF_M[i].size_Nd(); k++)
	  for (int l=0; l<HF_M[i](j,k).size_N(); l++)
	    for (int m=0; m<HF_M[i](j,k).size_Nd(); m++)
	      dbuff[++dp] = HF_M[i](j,k)(l,m);

    //common::gout<<"2)dsize after HF_M="<<dp+1<<endl;

    // Id.second
    for (map<int,double>::const_iterator i=Id.begin(); i!=Id.end(); i++) dbuff[++dp] = i->second;

    //common::gout<<"Id.second=";
    //for (map<int,double>::const_iterator i=Id.begin(); i!=Id.end(); i++) common::gout<<i->second<<" ";
    //common::gout<<endl;

    // DF_M
    if (DOsize>0){
      for (int i=0; i<DF_M.size_N(); i++){
	for (int j=0; j<DF_M.size_Nd(); j++){
	  for (int k=0; k<DF_M(i,j).size_N(); k++){
	    for (int l=0; l<DF_M(i,j).size_Nd(); l++){
	      dbuff[++dp] = DF_M(i,j)(k,l);
	    }
	  }
	}
      }
    }
    // Uc
    if (common::QHB2){
      for (int k=0; k<UC.size(); k++)  dbuff[++dp] = UC[k].u;
    }
    
    // strings
    strcpy(sbuff, RealSigma.c_str());
    strcpy(sbuff+(RealSigma.size()+1), ImagSigma.c_str());

    //common::gout<<my_rank<<":i-sizes="<<isize<<" "<<ip<<endl;
    //common::gout<<my_rank<<":d-sizes="<<dsize<<" "<<dp<<endl;

  }

  //common::gout<<my_rank<<":isize="<<isize<<endl;
  //common::gout<<my_rank<<":dsize="<<dsize<<endl;

  
  //common::gout<<my_rank<<": Before second bcast"<<endl;
  //MPI::COMM_WORLD.Bcast(ibuff, isize, MPI::INT, Master);
  MPI_Bcast(ibuff, isize, MPI_INT, Master, MPI_COMM_WORLD);
  //MPI::COMM_WORLD.Bcast(dbuff, dsize, MPI::DOUBLE, Master);
  MPI_Bcast(dbuff, dsize, MPI_DOUBLE, Master, MPI_COMM_WORLD);
  //MPI::COMM_WORLD.Bcast(sbuff, ssize, MPI::CHAR, Master);
  MPI_Bcast(sbuff, ssize, MPI_CHAR, Master, MPI_COMM_WORLD);
  //common::gout<<my_rank<<": After second bcast"<<endl;
  
  if (my_rank!=Master){
    
    Nm = ibuff[0];
    nsize = ibuff[1];
    N_ifl = ibuff[2];
    N_flavors = ibuff[3];
    Nvfl = ibuff[4];
    max_size = ibuff[5];
    Osize = ibuff[6];
    int max_dim = ibuff[7];
    
    int ip=8-1;
    // ifl_dim
    ifl_dim.resize(N_ifl);
    for (int i=0; i<ifl_dim.size(); i++) ifl_dim[i] = ibuff[++ip];
    // fl_from_ifl
    fl_from_ifl.resize(N_ifl,max_dim);
    for (int i=0; i<fl_from_ifl.size_N(); i++){                   
      for (int j=0; j<fl_from_ifl.size_Nd(); j++){
	fl_from_ifl(i,j) = ibuff[++ip];
      }
    }
    // tfl_index
    tfl_index.resize(N_ifl);
    for (int i=0; i<tfl_index.size(); i++){
      tfl_index[i].resize(ifl_dim[i],ifl_dim[i]);
      for (int j=0; j<tfl_index[i].size_N(); j++)
	for (int k=0; k<tfl_index[i].size_Nd(); k++)
	  tfl_index[i](j,k) = ibuff[++ip];
    }
    // vfl_index
    vfl_index.resize(N_ifl);
    for (int i=0; i<vfl_index.size(); i++){
      vfl_index[i].resize(ifl_dim[i],ifl_dim[i]);
      for (int j=0; j<vfl_index[i].size_N(); j++)
	for (int k=0; k<vfl_index[i].size_Nd(); k++)
	  vfl_index[i](j,k) = ibuff[++ip];
    }
    // v2fl_index
    v2fl_index.resize(N_ifl);
    for (int i=0; i<v2fl_index.size(); i++){
      int n = ifl_dim[i]*ifl_dim[i];
      for (int j=0; j<n; j++) v2fl_index[i].push_back( ibuff[++ip] );
    }
    // sign
    sign.resize(N_ifl);
    for (int i=0; i<sign.size(); i++){
      int n = ifl_dim[i]*ifl_dim[i];
      for (int j=0; j<n; j++) sign[i].push_back(ibuff[++ip]);
    }
    // conjg
    conjg.resize(N_ifl);
    for (int i=0; i<conjg.size(); i++){
      int n = ifl_dim[i]*ifl_dim[i];
      for (int j=0; j<n; j++) conjg[i].push_back(ibuff[++ip]);
    }
    // bfl_index
    bfl_index.resize(N_ifl);
    for (int i=0; i<bfl_index.size(); i++){
      int n = ifl_dim[i]*ifl_dim[i];
      for (int j=0; j<n; j++) bfl_index[i].push_back(ibuff[++ip]);
    }
    // gflip_index
    gflip_index.resize(N_ifl);
    for (int i=0; i<gflip_index.size(); i++) gflip_index[i] = ibuff[++ip];
    // Ns
    Ns.resize(nsize+1);
    for (int i=0; i<Ns.size(); i++) Ns[i] = ibuff[++ip];
    // Ks
    Ks.resize(nsize+1);
    for (int i=0; i<Ks.size(); i++) Ks[i] = ibuff[++ip];
    // msize
    msize_.resize(nsize+1);
    for (int i=0; i<msize_.size(); i++) msize_[i] = ibuff[++ip];
    // F_i
    F_i.resize(2*N_flavors,nsize+1);
    for (int i=0; i<F_i.size_N(); i++)
      for (int j=0; j<F_i.size_Nd(); j++)
	F_i(i,j) = ibuff[++ip];

    // F_M(,).sizes()
    function2D<pair<int,int> > F_M_sizes(2*N_flavors,nsize+1);
    for (int i=0; i<F_M_sizes.size_N(); i++){
      for (int j=0; j<F_M_sizes.size_Nd(); j++){
	int size_N = ibuff[++ip];
	int size_Nd = ibuff[++ip];
	
	//common::gout<<my_rank<<") F_M("<<i<<","<<j<<").size="<<size_N<<" "<<size_Nd<<endl;
	
	F_M_sizes(i,j) = make_pair( size_N, size_Nd );
      }
    }

    // Id.size
    int Idsize = ibuff[++ip];
    vector<int> Id_first(Idsize);
    // Id.first
    for (int i=0; i<Idsize; i++) Id_first[i] = ibuff[++ip];
    
    int ssize0 = ibuff[++ip];
    int Enes_size = ibuff[++ip];
    
    // DOsize
    DOsize = ibuff[++ip];

    if (DOsize>0){
      // DF_i(DOsize,nsize+1)
      DF_i.resize(DOsize,nsize+1);
      for (int i=0; i<DOsize; i++)
	for (int j=0; j<nsize+1; j++)
	  DF_i(i,j) = ibuff[++ip];
      // DF_inv
      DF_inv.resize(DOsize,nsize+1);
      for (int i=0; i<DOsize; i++)
	for (int j=0; j<nsize+1; j++)
	  DF_inv(i,j) = ibuff[++ip];
      // DF_M(,).sizes()
      DF_M.resize(DOsize,nsize+1);
      for (int i=0; i<DOsize; i++){
	for (int j=0; j<nsize+1; j++){
	  int _N_ = ibuff[++ip];
	  int _Nd_ = ibuff[++ip];
	  DF_M(i,j).resize(_N_,_Nd_);
	}
      }
    }
    
    if (common::QHB2){
      int UC_size = ibuff[++ip];
      UC.resize(UC_size);
      for (int k=0; k<UC.size(); k++){
	UC[k].a = ibuff[++ip];
	UC[k].b = ibuff[++ip];
	UC[k].c = ibuff[++ip];
	UC[k].d = ibuff[++ip];
      }
    }
    
    //debug
    //common::gout<<"Id.first=";
    //for (int i=0; i<Idsize; i++) common::gout<<Id_first[i]<<" ";
    //common::gout<<endl;
    //common::gout<<"ssize0="<<ssize0<<endl;
    
    if (ip+1!=isize) {cerr<<"ERROR: During BcastClusterData ip!=isize "<<ip<<" "<<isize<<endl; return 1;}

    //common::gout<<"ip+1="<<ip+1<<" isize="<<isize<<endl;
    
    // Done: N_unique_fl
    N_unique_fl = 0;
    for (int ifl=0; ifl<N_ifl; ifl++)
      for (size_t b=0; b<bfl_index[ifl].size(); b++)
	if (bfl_index[ifl][b]>N_unique_fl) N_unique_fl = bfl_index[ifl][b];
    N_unique_fl++;
    
    
    
    int dp=-1;
    // Sz
    Sz.resize(nsize+1);
    for (int i=0; i<Sz.size(); i++) Sz[i] = dbuff[++dp];
    // epsk
    epsk.resize(N_unique_fl);
    for (int i=0; i<epsk.size(); i++) epsk[i] = dbuff[++dp];
    // Ene
    /*
    Ene.resize(nsize+1,max_size);
    for (int i=0; i<Ene.size_N(); i++)
      for (int j=0; j<Ene.size_Nd(); j++)
	Ene(i,j) = dbuff[++dp];
    */
    Enes.resize(Enes_size);
    for (int i=0; i<Enes.size(); i++) Enes[i] = dbuff[++dp];
    // Ene0
    Ene0.resize(nsize+1,max_size);
    for (int i=0; i<Ene0.size_N(); i++)
      for (int j=0; j<Ene0.size_Nd(); j++)
	Ene0(i,j) = dbuff[++dp];
    // Spin
    Spin.resize(nsize+1,max_size);
    for (int i=0; i<Spin.size_N(); i++)
      for (int j=0; j<Spin.size_Nd(); j++)
	Spin(i,j) = dbuff[++dp];
    
    //common::gout<<"3)dsize before F_M="<<dp+1<<endl;
    // F_M
    F_M.resize(2*N_flavors,nsize+1);
    for (int i=0; i<F_M.size_N(); i++){
      for (int j=0; j<F_M.size_Nd(); j++){
	F_M(i,j).resize(F_M_sizes(i,j).first,F_M_sizes(i,j).second);
	for (int k=0; k<F_M(i,j).size_N(); k++){
	  for (int l=0; l<F_M(i,j).size_Nd(); l++){
	    F_M(i,j)(k,l) = dbuff[++dp];
	  }
	}
      }
    }
    //common::gout<<"3)dsize before HF_M="<<dp+1<<endl;
    
    // HF_M
    HF_M.resize(Osize);
    for (int i=0; i<HF_M.size(); i++){
      HF_M[i].resize(N_unique_fl,nsize+1);
      for (int j=0; j<HF_M[i].size_N(); j++){
	for (int k=0; k<HF_M[i].size_Nd(); k++){
	  HF_M[i](j,k).resize(msize_[k],msize_[k]);
	  for (int l=0; l<HF_M[i](j,k).size_N(); l++){
	    for (int m=0; m<HF_M[i](j,k).size_Nd(); m++){
	      HF_M[i](j,k)(l,m) = dbuff[++dp];
	    }
	  }
	}
      }
    }
    //common::gout<<"3)dsize after HF_M="<<dp+1<<endl;

    // Id.second
    for (int i=0; i<Idsize; i++) Id[Id_first[i]] = dbuff[++dp];

    if (DOsize>0){
      // DF_M
      for (int i=0; i<DF_M.size_N(); i++){
	for (int j=0; j<DF_M.size_Nd(); j++){
	  for (int k=0; k<DF_M(i,j).size_N(); k++){
	    for (int l=0; l<DF_M(i,j).size_Nd(); l++){
	      DF_M(i,j)(k,l) = dbuff[++dp];
	    }
	  }
	}
      }
    }
    
    // Uc
    if (common::QHB2){
      for (int k=0; k<UC.size(); k++)  UC[k].u = dbuff[++dp];
    }
    
    if (dp+1!=dsize) {cerr<<"ERROR: During BcastClusterData dp!=dsize "<<dp<<" "<<dsize<<endl; return 1;}

    RealSigma = sbuff;
    ImagSigma = &(sbuff[ssize0+1]);

    // Need to set up Eind index for Enes!!
    Eind.resize(nsize+1,max_size);
    int iind=0;
    for (int i=1; i<=nsize; i++){
      for (int j=0; j<msize_[i]; j++){
	Eind(i,j) = iind;
	iind++;
      }
    }
  }
  
  delete[] ibuff;
  delete[] dbuff;
  delete[] sbuff;
  
  if (my_rank!=Master){
    // Done: vfli_index
    vfli_index.resize(Nvfl);
    for (int i=0; i<N_ifl; i++){
      int ii=0;
      for (int i1=0; i1<ifl_dim[i]; i1++){
	for (int i2=0; i2<ifl_dim[i]; i2++){
	  int jj = vfl_index[i][i1][i2];
	  vfli_index[jj] = TBath(i,ii);
	  ii++;
	}
      }
    }

    // Done: ifl_from_fl, bfl_from_fl
    ifl_from_fl.resize(N_flavors);
    bfl_from_fl.resize(N_flavors);
    for (int ifl=0; ifl<N_ifl; ifl++){
      for (int jb=0; jb<ifl_dim[ifl]; jb++){
	int tfl = fl_from_ifl(ifl,jb);
	ifl_from_fl[tfl] = ifl;
	bfl_from_fl[tfl] = jb;
      }
    }

    // Done: fl_deg
    fl_deg.resize(N_unique_fl);
    fl_deg=0;
    for (int ifl=0; ifl<N_ifl; ifl++)
      for (size_t b=0; b<bfl_index[ifl].size(); b++) fl_deg[bfl_index[ifl][b]]++;


    // gfl_tmp
    map<int,deque<int> > gfl_tmp;
    for (int ifl=0; ifl<N_ifl; ifl++){
      int ii = gflip_index[ifl];
      gfl_tmp[ii].push_back(ifl);
    }
    int sz=0;
    for (map<int,deque<int> >::iterator ia=gfl_tmp.begin(); ia!=gfl_tmp.end(); ia++){
      deque<int>& d = ia->second;
      sz += d.size()*(d.size()-1)/2;
    }
  
    // gflip_fl, gflip_ifl
    gflip_fl.resize(sz+1,N_flavors);
    gflip_ifl.resize(sz);
    for (int iu=0; iu<sz+1; iu++){
      for (int j=0; j<gflip_fl.size_Nd(); j++) gflip_fl[iu][j] = j;
    }
    int iu=0;
    for (map<int,deque<int> >::iterator ia=gfl_tmp.begin(); ia!=gfl_tmp.end(); ia++){
      deque<int>& d = ia->second;
      for (size_t j1=0; j1<d.size(); j1++){
	for (size_t j2=j1+1; j2<d.size(); j2++){
	  int ifa = d[j1];
	  int ifb = d[j2];
	  if (ifl_dim[ifa]!=ifl_dim[ifb]){cerr<<"Dimensions of similar baths have to be the same!"<<endl; return 1;}
	  gflip_ifl[iu] = make_pair(ifa,ifb);
	  for (int ib=0; ib<ifl_dim[ifa]; ib++){
	    int fa = fl_from_ifl[ifa][ib];
	    int fb = fl_from_ifl[ifb][ib];
	    gflip_fl[iu][fa] = fb;
	    gflip_fl[iu][fb] = fa;
	    gflip_fl[sz][fa] = fb;
	    gflip_fl[sz][fb] = fa;
	  }
	  iu++;
	}
      }
    }
    // ifl_equivalent
    ifl_equivalent.resize(N_ifl,N_ifl);
    for (int ifl=0; ifl<N_ifl; ifl++){
      for (int jfl=0; jfl<N_ifl; jfl++){
	ifl_equivalent[ifl][jfl] = bfl_index[ifl]==bfl_index[jfl];
      }
    }
    
  }
  //DEBUGGING

  //common::gout<<"Nm="<<Nm<<endl;
  //common::gout<<"nsize="<<nsize<<endl;
  //common::gout<<"N_ifl="<<N_ifl<<endl;
  //common::gout<<"N_flavors="<<N_flavors<<endl;
  //common::gout<<"Nvfl="<<Nvfl<<endl;
  //common::gout<<"max_size="<<max_size<<endl;
  //common::gout<<"Osize="<<Osize<<endl;
  
  //common::gout<<"Sz=";
  //for (int i=0; i<Sz.size(); i++) common::gout<<Sz[i]<<" ";
  //common::gout<<endl;
  //common::gout<<"epsk=";
  //for (int i=0; i<epsk.size(); i++) common::gout<<epsk[i]<<" ";
  //common::gout<<endl;
  //common::gout<<"Ene=";
  //for (int i=0; i<Ene.size_N(); i++)
  //  for (int j=0; j<Ene.size_Nd(); j++)
  //    common::gout<<Ene(i,j)<<" ";
  //common::gout<<endl;
  
  //common::gout<<"RealSigma="<<RealSigma<<endl;
  //common::gout<<"ImagSigma="<<ImagSigma<<endl;
  return 0;
}

void BCastDelta(int my_rank, int Master, mesh1D& iom_large, function2D<dcomplex>& Deltaw, vector<pair<double,double> >& ah)
{
  int isize[4];
  double* ah_buf;
  if (my_rank==Master){
    isize[0] = iom_large.size();
    isize[1] = Deltaw.size_N();
    isize[2] = Deltaw.size_Nd();
    isize[3] = ah.size();
  }

  // only sizes are brodcast
  //MPI::COMM_WORLD.Bcast(isize, 4, MPI::INT, Master);
  MPI_Bcast(isize, 4, MPI_INT, Master, MPI_COMM_WORLD);
  // Vectors might not be implemented in a simple way. Changing vector to plane array.
  ah_buf = new double[isize[3]*2];
  
  if (my_rank==Master){
    for (int i=0; i<ah.size(); i++){
      ah_buf[2*i] = ah[i].first;
      ah_buf[2*i+1]=ah[i].second;
    }
  }else{
    iom_large.resize(isize[0]);
    Deltaw.resize(isize[1],isize[2]);
  }

  //MPI::COMM_WORLD.Bcast(iom_large.MemPt(), isize[0], MPI::DOUBLE, Master);
  MPI_Bcast(iom_large.MemPt(), isize[0], MPI_DOUBLE, Master, MPI_COMM_WORLD);
  int Deltaw_size = isize[1]*isize[2];
  //MPI::COMM_WORLD.Bcast(Deltaw.MemPt(), Deltaw_size, MPI::DOUBLE_COMPLEX, Master);
  MPI_Bcast(Deltaw.MemPt(), Deltaw_size, MPI_DOUBLE_COMPLEX, Master, MPI_COMM_WORLD);
  //MPI::COMM_WORLD.Bcast(ah_buf, isize[3]*2, MPI::DOUBLE, Master);
  MPI_Bcast(ah_buf, isize[3]*2, MPI_DOUBLE, Master, MPI_COMM_WORLD);
  
  if (my_rank!=Master){
    iom_large.SetUp(0);
    ah.resize(isize[3]);
    for (int i=0; i<ah.size(); i++){
      ah[i].first = ah_buf[2*i];
      ah[i].second = ah_buf[2*i+1];
    }
  }
  delete[] ah_buf;
};


#else
int ClusterData::BcastClusterData(int my_rank, int Master)
{return 0;}

void BCastDelta(int my_rank, int Master, mesh1D& iom_large, function2D<dcomplex>& Deltaw, vector<pair<double,double> >& ah){};

#endif
