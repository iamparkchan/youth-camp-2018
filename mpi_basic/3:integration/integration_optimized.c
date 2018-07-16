#include <stdio.h>
#include <math.h>
#include <time.h>

#define N 100000000
#define A 0.0
#define B 1.0

#define f(x) (x*x)

double solution(double a, double b)
{
	return b*b*b/3.0 - a*a*a/3.0;
}

int main()
{
	int i;
	double sum = 0.0, step = (B - A) / N;
	clock_t start, end;

	start = clock();
	for (i = 0; i < N; i++) {
		sum += f(i*step);	
	}
	sum *= step;
	end = clock();

	printf("result=%e\n", sum);
	printf("error=%e\n", fabs(sum - solution(A, B)));
	printf("elapsed time=%f\n", (double)(end - start)/CLOCKS_PER_SEC);

	return 0;
}
