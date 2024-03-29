// @Copyright 2007 Kristjan Haule
// 
#ifdef _MPI
void MPI_Init(int argc, char *argv[], int& my_rank, int& mpi_size)
{
  /* // Unfortunately C++ binding is deprecated...
  MPI::Init(argc, argv);
  my_rank = MPI::COMM_WORLD.Get_rank();
  mpi_size = MPI::COMM_WORLD.Get_size();
  */
  ::MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
}
void MPI_Finish()
{
  //MPI::Finalize();
  MPI_Finalize();
}
void MPI_Allgather_Self_energy(int my_rank, function1D<int>& send_size, function1D<int>& start, double* au_Sign)
{
  //MPI::COMM_WORLD.Allgatherv(MPI_IN_PLACE, send_size[my_rank], MPI_DOUBLE, au_Sign, send_size.MemPt(), start.MemPt(), MPI_DOUBLE);
  MPI_Allgatherv(MPI_IN_PLACE, send_size[my_rank], MPI_DOUBLE, au_Sign, send_size.MemPt(), start.MemPt(), MPI_DOUBLE, MPI_COMM_WORLD);
}
void MPI_Gather_A00(int my_rank, int Master, double* ph_P_A00, function1D<int>& send_size, function1D<int>& start)
{
  //MPI::COMM_WORLD.Gatherv((my_rank==Master) ? MPI_IN_PLACE : ph_P_A00+start[my_rank], send_size[my_rank], MPI_DOUBLE, ph_P_A00, send_size.MemPt(), start.MemPt(), MPI_DOUBLE, Master);
  MPI_Gatherv((my_rank==Master) ? MPI_IN_PLACE : ph_P_A00+start[my_rank], send_size[my_rank], MPI_DOUBLE, ph_P_A00, send_size.MemPt(), start.MemPt(), MPI_DOUBLE, Master, MPI_COMM_WORLD);
}
#else
void MPI_Init(int argc, char *argv[], int& my_rank, int& mpi_size)
{
  my_rank = 0;
  mpi_size = 1;
}
void MPI_Finish(){}
void MPI_Allgather_Self_energy(int my_rank, function1D<int>& send_size, function1D<int>& start, double* au_Sign){}
void MPI_Gather_A00(int my_rank, int Master, double* ph_P_A00, function1D<int>& send_size, function1D<int>& start){}
#endif
