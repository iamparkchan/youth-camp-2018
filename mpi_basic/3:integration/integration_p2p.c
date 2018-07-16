#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define N 100000000
#define A 0.0
#define B 1.0

double f(double x)
{
    return x*x;
}

double solution(double a, double b)
{
    return b*b*b/3.0 - a*a*a/3.0;
}

int main(int argc, char *argv[])
{
    int i;
    double my_sum = 0.0, sum = 0.0, h = (B - A) / N;
    double start, end;
    int rank, size;
    int unit, i_s, i_e;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    start = MPI_Wtime();
    unit = (N + size - 1) / size;
    i_s = rank*unit;
    if (rank == size - 1) {
        i_e = N;
    } else {
        i_e = i_s + unit;
    }
    for (i = i_s; i < i_e; i++) {
        my_sum += f(A + i*h);
    }
    if (0 != rank) {
        MPI_Send(&my_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else {
        sum += my_sum;
        for (i = 1; i < size; i++) {
            MPI_Recv(&my_sum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += my_sum;
        }
        sum *= h;
    }
    end = MPI_Wtime();

    if (0 == rank) {
        printf("result=%e\n", sum);
        printf("error=%e\n", fabs(sum - solution(A, B)));
        printf("elapsed time=%f\n", end - start);
    }

    MPI_Finalize();

    return 0;
}
