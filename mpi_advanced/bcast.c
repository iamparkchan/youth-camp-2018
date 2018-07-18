#include <stdio.h>
#include <mpi.h>
int main(int argc, char *argv[])
{
  int i, nrank, nprocs, ROOT = 0;
  int buf[4] = { 0, 0, 0, 0 };
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &nrank);
  if (nrank == ROOT) {
     buf[0] = 1; buf[1] = 2; buf[2] = 3; buf[3] = 4;
  }
  printf("rank (%d) : Before : ", nrank);
  for (i=0; i<4; i++) printf(" %d", buf[i]);
  printf("\n");
  MPI_Bcast(buf, 4, MPI_INT, ROOT, MPI_COMM_WORLD);
  printf("rank (%d) : After : ", nrank);
  for (i=0; i<4; i++) printf(" %d", buf[i]);
  printf("\n");
  MPI_Finalize();
  return 0;
}

