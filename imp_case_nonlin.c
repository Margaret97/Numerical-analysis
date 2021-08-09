#include "pch.h"
#include <iostream>
#include <stdio.h> 
#include <stdlib.h> 
#define _USE_MATH_DEFINES 
#include <math.h> 
#include <float.h>
#include <malloc.h>
#define  eps 1e-9
double F(double v);
double F(double v)
{
	return -0.5*v*v;
}
double dF(double v);
double dF(double v)
{
	return -v;
}



int main(void)
{
	FILE*data;
	FILE*norm;
	double h, tau, omega;
	int n, m, k;
	double *v = NULL;
	double *r = NULL;
	double *lambda = NULL;
	double *u_i = NULL;
	double *l = NULL;
	double*step = NULL;
	int i, j;
	double u;
	i = 0; j = 1; k = 1;
	double max; 
	double prev;
	omega = 1;
	double rate_Ch, rate_Lh, v_Ch, v_Lh;
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
	r = (double*)malloc((m + 1) * sizeof(double));
	l = (double*)malloc((m + 1) * sizeof(double));
	u_i = (double*)malloc(m * sizeof(double));
	lambda = (double*)malloc(m * sizeof(double));
	step= (double*)malloc((m + 1) * sizeof(double));
	for (i = 0; i <= m / 2 ; i++)
	{
		v[i] = 1;
	}
	for (i = m / 2 + 1; i <= m; i++)
	{
		v[i] = 0;
	}

	for (i = 0; i <= m; i++)
	{
		step[i] = v[i];
	}

	for (j = 1; j <= n; j++)
	{
		max = 1;
		k = 1;
		while(max>eps)
		{
			for (i = 0; i <= m; i++)
			{
				if (i == 0)
				{
					prev = 1;
					l[i] = 1;
					r[i] = (-step[i] - 0.5*tau*F(step[i + 1]) / h + v[i] + 0.5*tau*F(prev)/h + 0.125*omega*(v[i+2]-4*v[i+1]+6*v[i]-4+1))/ l[i];
					if (fabs(r[i])>1e300)
					{
						printf("The solution goes to infinity\n");
						printf("Iteration %d the layer %d point %d", k, j, i);
						return -1;
					}
				}
				else
				{
					l[i] = 1 - lambda[i - 1] * u_i[i - 1];
					if (fabs(l[i]) < 1e-30)
					{
						printf("There is the problem: l_i is small\n");
						printf("Layer %d\t point %d\t step %d\n", j, i, k);
						return -1;
					}

					if (i == 1)
					{
						r[i] = (-step[i] - 0.5*tau*F(step[i + 1]) / h + v[i] + 0.5*tau*F(prev)/ h + 0.125*omega*(v[i + 2] - 4 * v[i + 1] + 6 * v[i] - 4*v[i-1] + 1) - lambda[i - 1] * r[i - 1]) / l[i];
						if (fabs(r[i]) > 1e300)
						{
							printf("The solution goes to infinity\n");
							printf("Iteration %d the layer %d point %d", k, j, i);
							return -1;
						}
						//printf("r[%d]=%e\n", i, r[i]);
					}
					else
					{
						if (i == m - 1)
						{
							r[i] = (-step[i] - 0.5*tau*F(step[i+1]) / h + v[i] + 0.5*tau*F(prev) / h + 0.125*omega*(-4 * v[i + 1] + 6 * v[i] - 4 * v[i - 1] + v[i - 2]) - lambda[i - 1] * r[i - 1]) / l[i];
							if (fabs(r[i]) > 1e300)
							{
								printf("The solution goes to infinity\n");
								printf("Iteration %d the layer %d point %d", k, j, i);
								return -1;
							}
						}
						else 
						{
							if (i == m)
							{
								r[i] = (-step[i] - 0.5*tau*F(0) / h + v[i] + 0.5*tau*F(prev) / h + 0.125*omega*(6 * v[i] - 4 * v[i - 1] + v[i - 2]) - lambda[i - 1] * r[i - 1]) / l[i];
								if (fabs(r[i]) > 1e300)
								{
									printf("The solution goes to infinity\n");
									printf("Iteration %d the layer %d point %d", k, j, i);
									return -1;
								}

							}
							else r[i] = (-step[i] - 0.5*tau*F(step[i + 1]) / h + v[i] + 0.5*tau*F(prev) / h + 0.125*omega*(v[i + 2] - 4 * v[i + 1] + 6 * v[i] - 4 * v[i - 1] + v[i - 2]) - lambda[i - 1] * r[i - 1]) / l[i];
							if (fabs(r[i]) > 1e300)
							{
								printf("The solution goes to infinity\n");
								printf("Iteration %d the layer %d point %d", k, j, i);
								return -1;
							}
						}
						
					}
					
				}
				
				prev = step[i];
				step[i] = step[i] + r[i];
				if (i < m)
				{
					u_i[i] = 0.5*tau*dF(step[i+1]) / (h*l[i]);
					lambda[i] = -0.5*tau*dF(prev)/ h;
					
				}
			}
			//printf("\n");
			max = fabs(r[m]);
			//printf("k=%d\n", k);
			//printf("r[%d]=%lf\n", m, r[m]);
			for (i = m - 1; i >= 0; i--)
			{
				r[i] = r[i] - u_i[i] * r[i + 1];
				//printf("r[%d]=%e\n", i, r[i]);
				if (fabs(r[i]) > max)
				{
					max = fabs(r[i]);
				}
			}
			

			k++;
			if (k > 100)
			{
				printf("Method diverges\n");
				return -1;
			}

		}

		for (i = 0; i <= m; i++)
		{
			v[i] = step[i];
		}
		/*fprintf_s(data,"Layer %d\t amount of iterations %d\n", j, k);
		for (i = 0; i <= m; i++)
		{
			fprintf_s(data,"%e\n", v[i]);
		}*/
		
	}
	

	for (i = 0; i <= m; i++)
	{
		if ((-2 + i * h) <= -1)
		{
			u = 1;

		}
		else
		{
			if ((-2 + i * h) > -1 && (-2 + i * h) < 0)
				u = 2 - i * h;
			else u = 0;
		}
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

	for(i=0;i<=m;i++)
	{
		fprintf_s(data, "%lf\t%e\n", -2 + i * h, v[i]);
	}
	fprintf_s(norm, "$%e$\t&$%e$\t&$%e$\t&$%e$\n", rate_Ch, rate_Lh, (rate_Ch / v_Ch), (rate_Lh / v_Lh));
	/*for(i=0;i<=m;i++)
	{
		printf("%lf\n", l[i]);
	}*/
		free(r);
		free(step);
		free(lambda);
		free(l);
		free(u_i);
		free(v);
		fclose(data);
		fclose(norm);
		return 0;

}



