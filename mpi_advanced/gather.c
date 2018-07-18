#include <stdio.h>
#include <mpi.h>
int main(int argc, char *argv[])
{
  int i, nprocs, nrank;
  int isend, irecv[4];
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &nrank);
  isend = nrank + 1;
  printf("rank (%d) : isend = %d ", nrank, isend);
  MPI_Gather(&isend, 1, MPI_INT, irecv, 1, MPI_INT, 0,
    MPI_COMM_WORLD);
  if (nrank == 0) {
    printf("\n");
    for (i=0; i<3; i++)
    printf("rank (%d) : irecv[%d] = %d\n",
    nrank, i, irecv[i]);
  }
  printf("\n");
  MPI_Finalize();
  return 0;
}

