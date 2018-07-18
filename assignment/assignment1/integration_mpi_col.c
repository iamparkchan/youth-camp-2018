#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define N 1000000000
#define A 0.0
#define B 1.0

#define f(x) (x*x)

/// example
// $ mpirun -np 4 ./integration 
// result=3.333333e-01
// error=5.000362e-10
// elapsed time=1.561990

double solution(double a, double b)
{
	return b*b*b/3.0 - a*a*a/3.0;
}

int main(int argc, char *argv[])
{
	int i;
	double my_sum = 0.0, sum = 0., tmp=0., step = (B - A) / N;
	double start, end;
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	start = MPI_Wtime();
	for (i = rank; i < N; i += size) {
		my_sum += f(i*step);	
	}
	MPI_Reduce(&my_sum, &tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	sum = tmp*step;
/*	if (0 != rank) {
		MPI_Send(&my_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	} else {
		sum += my_sum;
		for (i = 1; i < size; i++) { 
			MPI_Recv(&my_sum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			sum += my_sum;
		}
		sum *= step;
	}*/
	end = MPI_Wtime();

	if (0 == rank) {
		printf("result=%e\n", sum);
		printf("error=%e\n", fabs(sum - solution(A, B)));
		printf("elapsed time=%f\n", end - start);
	}

	MPI_Finalize();

	return 0;
}
