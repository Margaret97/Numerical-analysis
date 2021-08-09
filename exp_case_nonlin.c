#include "pch.h"
#include <iostream>
#include <stdio.h> 
#include <stdlib.h> 
#define _USE_MATH_DEFINES 
#include <math.h> 
#include<float.h>
int sign(double v); // the sign function

int sign(double v) 
{
	if (v > 1e-30)
	{
		return 1;
	}
	if (v < -1e-30)
	{
		return -1;
	}
	else return 0;
}
int main(void)
{
	FILE*nonlin_case;
	FILE*norm;
	double h, tau, u; // tau, h are the integration steps in time and space, u is the analytical solution
	int n, m; // the numbers of integration steps in time and space
	double *v; // the numerical solution
	double prev, current;
	int i, j;
	i = 0; j = 1;
	double rate_Ch, rate_Lh, v_Ch, v_Lh; // the C,L1-norms of the difference between v and u, the C,L1-norms of v  (description in the file "partial_differential_equation.pdf")
	rate_Ch = 0;
	rate_Lh = 0;
	v_Lh = 0; v_Ch = 0;
	fopen_s(&nonlin_case, "C://Users//Margarita//data//task_2.txt", "w");
	fopen_s(&norm, "C://Users//Margarita//data//norm.txt", "w");
	if (nonlin_case == NULL)
	{
		printf("Error opening file nonlin_case.\n");
		return -1;
	}
	if (norm == NULL)
	{
		printf("Error opening file norm.\n");
		return -1;
	}
	printf("Enter h:\n");
	scanf_s("%lf", &h);
	printf("Enter tau:\n");
	scanf_s("%lf", &tau);
	m = int(4 / h);
	printf("m=%d\n", m);
	n = int(1 / tau);
	printf("n=%d\n", n);
	v = (double*)malloc((m + 2*n + 1) * sizeof(double));
	for (i = 0; i <= n + m / 2; i++) // we give the initial condition
	{
		v[i] = 1;
	}
	for (i = m / 2 + n + 1; i <= m + 2*n; i++) // we give the initial condition
	{
		v[i] = 0;
	}
	for (j = 1; j <= n; j++) // we integrate the equation
	{
		prev = v[j - 1];
		for (i = j; i <= m + 2*n - j; i++) 
		{
			switch (sign(v[i]))
			{
				case 1:
					prev = v[i];
					v[i] = v[i] + tau * 0.5*(v[i + 1] * v[i + 1] - v[i] * v[i]) / h;
					if (fabs(v[i])>1e300)
					{
						printf("The solution goes to infinity\n");
						printf("The layer %d step %d\n", j, i);
						return -1;
					}
					break;
				case -1:
					current = v[i];
					v[i] = v[i] + tau * 0.5*(v[i] * v[i] - prev * prev) / h;
					if (fabs(v[i]) > 1e300)
					{
						printf("The solution goes to infinity\n");
						printf("The layer %d step %d\n", j, i);
						return -1;
					}
					prev = current;
					break;
				case 0:
					current = v[i];
					v[i] = v[i] + tau * 0.25*(v[i + 1] * v[i + 1] - prev * prev) / h;
					if (fabs(v[i]) > 1e300)
					{
						printf("The solution goes to infinity\n");
						printf("The layer %d step %d\n", j, i);
						return -1;
					}
					prev = current;
					break;
				default: printf("Error for calculating sign\n");
					return -1;
			}
			
		}

	}
	for (i = 0; i <= m; i++) // we calculate the norms taking into account the solution u
	{
		if ((-2 + i * h) <= -1)
		{
			u = 1;

		}
		if ((-2 + i * h > -1) && (-2 + i * h <= 0))
		{
			u = -i * h + 2;
		}

		if (-2 + i * h > 0)
		{
			u = 0;
		}

		if (fabs(v[i+n] - u) > rate_Ch)
		{
			rate_Ch = fabs(v[i+n] - u);
		}
		if (fabs(v[i+n]) > v_Ch)
		{
			v_Ch = fabs(v[i+n]);
		}
		rate_Lh = h * fabs(v[i+n] - u) + rate_Lh;
		v_Lh = h * fabs(v[i+n]) + v_Lh;
	}
	printf("rate_Ch=%lf\n", rate_Ch);
	printf("v_Ch=%lf\n", v_Ch);
	fprintf_s(norm, "&$%e$\t&$%e$\t&$%e$\t&$%e$\n", rate_Ch, rate_Lh, (rate_Ch / v_Ch), (rate_Lh / v_Lh));
	for (i = 0; i <= m; i++)
	{
		fprintf_s(nonlin_case, "%lf\t%lf\n", -2 + i * h, v[i+n]);
	}
	

	printf("Everything is ok!\n");
	free(v);
	fclose(nonlin_case);
	fclose(norm);
	return 0;

}



