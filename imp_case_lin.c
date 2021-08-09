#include "pch.h"
#include <iostream>
#include <stdio.h> 
#include <stdlib.h> 
#define _USE_MATH_DEFINES 
#include <math.h> 
#include <float.h>
#include <malloc.h>


int main(void)
{
	FILE*data;
	FILE*norm;
	double h, tau; // tau, h are the integration steps in time and space
	int n, m, omega; // n, m are the numbers of integration steps in time and space, omega is the number in the implicit scheme
	double *v = NULL; // the numerical solution v
	double *lambda = NULL; // the elements below the main diagonal of the lower triangular matrix L
	double *u_i = NULL; // the elements above the main diagonal of the upper triangular matrix U, all elements of the main diagonal are units
	double *l = NULL; // the main diagonal of the lower triangular matrix L
	int i, j, u; // u is the analytical solution
	i = 0; j = 1; 
	double prev_1, prev_2, current;
	omega = 0;
	double rate_Ch, rate_Lh, v_Ch, v_Lh; // the C,L1-norms of the difference between v and u, the C,L1-norms of v  (description in the file "partial_differential_equation.pdf")
	rate_Ch = 0;
	rate_Lh = 0;
	v_Lh = 0; v_Ch = 0;
	
	fopen_s(&data, "C://Users//Margarita//data//task_2.txt", "w");
	fopen_s(&norm, "C://Users//Margarita//data//norm.txt", "w");
	if (data == NULL)
	{
		printf("Error opening file data.\n");
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
	v = (double*)malloc((m + 1) * sizeof(double));
	l = (double*)malloc((m + 1) * sizeof(double));
	u_i = (double*)malloc(m * sizeof(double));
	lambda = (double*)malloc(m * sizeof(double));
	for (i = 0; i <= m; i++) // LU decomposition
	{
		if (i == 0)
		{
			l[i] = 1;
		}
		else l[i] = 1 - lambda[i - 1] * u_i[i - 1];
		if (fabs(l[i]) < 1e-30)
		{
			printf("There is the problem: l[%d] is small\n", i);
			return -1;
		}
		if (i < m)
		{
			u_i[i] = -0.25*tau/(h*l[i]);
			lambda[i] = 0.25*tau / h;

		}

	}
	for (i = 0; i <= m / 2; i++) // we give the initial condition
	{
		v[i] = 1;
	}
	for (i = m / 2 + 1; i <= m; i++) // we give the initial condition
	{
		v[i] = 0;
	}

	
	for (j = 1; j <= n; j++) // we integrate the equation
	{
			prev_1 = 1;
			prev_2 = 1;
			for (i = 0; i <= m; i++)
			{
				if (i == 0)
				{
					current= v[i];
					v[i] = (v[i] - 0.25*tau/h + 0.125*omega*(v[i + 2] - 4 * v[i + 1] + 6 * v[i] - 4*prev_1 + prev_2)) / l[i];
					if (fabs(v[i]) > 1e300)
					{
						printf("The solution goes to infinity\n");
						printf("The layer %d step %d\n", j, i);
						return -1;
					}
					prev_1 = current;
					prev_2 = 1;
				}
				else
				{
						if (i == m - 1)
						{
							current = v[i];
							v[i] = (v[i] + 0.125*omega*(-4 * v[i + 1] + 6 * v[i] - 4 * prev_1 + prev_2) - lambda[i - 1] * v[i - 1]) / l[i];
							if (fabs(v[i]) > 1e300)
							{
								printf("The solution goes to infinity\n");
								printf("The layer %d step %d\n", j, i);
								return -1;
							}
							prev_2 = prev_1;
							prev_1 = current;

						}
						else
						{
							if (i == m)
							{	
								v[i] = (v[i]  + 0.125*omega*(6 * v[i] - 4 * prev_1 + prev_2) - lambda[i - 1] * v[i - 1]) / l[i];
								if (fabs(v[i]) > 1e300)
								{
									printf("The solution goes to infinity\n");
									printf("The layer %d step %d\n", j, i);
									return -1;
								}

							}
							else 
								current = v[i];
								v[i] = (v[i] + 0.125*omega*(v[i + 2] - 4 * v[i + 1] + 6 * v[i] - 4 * prev_1 + prev_2) - lambda[i - 1] * v[i - 1]) / l[i];
								if (fabs(v[i]) > 1e300)
								{
									printf("The solution goes to infinity\n");
									printf("The layer %d step %d\n", j, i);
									return -1;
								}
								prev_2 = prev_1;
								prev_1 = current; 
								
						}

				}

				

			}
			
			for (i = m - 1; i >= 0; i--)
			{
				v[i] = v[i] - u_i[i] * v[i + 1];
			}


			

		}

		
	for(i = 0; i <= m; i++) // we calculate the norms taking into account the solution u
	{
		if ((-2 + i * h) <= -0.5)
		{
			u = 1;

		}
		else u = 0;
		if (fabs(v[i] - u) > rate_Ch)
		{
			rate_Ch = fabs(v[i] - u);
		}
		if (fabs(v[i]) > v_Ch)
		{
			v_Ch = fabs(v[i]);
		}
		rate_Lh = h * fabs(v[i] - u) + rate_Lh;
		v_Lh = h * fabs(v[i]) + v_Lh;
	}

	printf("rate_Ch=%lf\n", rate_Ch);

	for (i = 0; i <= m; i++)
	{
		fprintf_s(data, "%lf\t%e\n", -2 + i * h, v[i]);
	}
	fprintf_s(norm, "$%e$\t$%e$\t$%e$\t$%e$\n", rate_Ch, rate_Lh, (rate_Ch / v_Ch), (rate_Lh / v_Lh));
	
	
	free(lambda);
	free(l);
	free(u_i);
	free(v);
	fclose(data);
	fclose(norm);
	return 0;

}



