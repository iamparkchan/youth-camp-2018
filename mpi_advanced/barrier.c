#include <stdio.h>
#include <mpi.h>
int main (int argc, char *argv[])
{
  int nRank, nProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  if (nRank == 0) {
     printf("############\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  printf("nRank(%d) : Hello World.\n", nRank);
  MPI_Finalize();
  return 0;
}

