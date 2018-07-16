#include <stdio.h>
#include <mpi.h>
#define N 100
int main(int argc, char *argv[])
{
	char ping[N] = "ping", pong[N] = "pong", received[N];
	int my_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	if (0 == my_rank) {
		MPI_Send(ping, N, MPI_CHAR, 1, 9, MPI_COMM_WORLD);
		MPI_Recv(received, N, MPI_CHAR, 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	} else {
		MPI_Recv(received, N, MPI_CHAR, 0, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(pong, N, MPI_CHAR, 0, 99, MPI_COMM_WORLD);
	}
	printf("%d: I received %s.\n", my_rank, received);
	MPI_Finalize();
	return 0;
}
