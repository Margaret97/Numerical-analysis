#include "pch.h"
#include <iostream>
#include <stdio.h> 
#include <stdlib.h> 
#define _USE_MATH_DEFINES 
#include <math.h> 
#include<float.h>

int main(void)
{	
	FILE*lin_case;
	FILE*norm;
	double h, tau; // tau, h are the integration steps in time and space
	int n, m; // the numbers of integration steps in time and space
	double *v; // the numerical solution
	int i, j, u; // u is the analytical solution
	i = 0; j = 1;
	double rate_Ch, rate_Lh, v_Ch, v_Lh; // the C,L1-norms of the difference between v and u, the C,L1-norms of v  (description in the file "partial_differential_equation.pdf")
	rate_Ch = 0;
	rate_Lh = 0;
	v_Lh = 0; v_Ch = 0;
	fopen_s(&lin_case,"C://Users//Margarita//data//task_2.txt", "w");
	fopen_s(&norm, "C://Users//Margarita//data//norm.txt", "w");
	if (lin_case == NULL)
	{
		printf("Error opening file lin_case.\n");
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
	v =(double*)malloc((m+n+1)*sizeof(double));
	for (i = 0; i <= m / 2 ; i++) // we give the initial condition
	{
		v[i] = 1;
	}
	for (i = m / 2 + 1; i <= m + n; i++) // we give the initial condition
	{
		v[i] = 0;
	}
	for (j = 1; j <= n; j++) // we integrate the equation
	{
		for(i=0;i<=m+n-j;i++)
		{
			v[i] = v[i] + tau * 0.5 * (v[i + 1] - v[i])/ h;
		}
		
	}
	for (i = 0; i <= m; i++) // we calculate the norms taking into account the solution u
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
		if(fabs(v[i]) > v_Ch)
		{
			v_Ch = fabs(v[i]);
		}
		rate_Lh = h * fabs(v[i] - u) + rate_Lh;
		v_Lh = h * fabs(v[i]) + v_Lh;
	}
	printf("rate_Ch=%e\n", rate_Ch);
	fprintf_s(norm, "$%e$\t&$%e$\t&$%e$\t&$%e$\n", rate_Ch, rate_Lh, (rate_Ch / v_Ch), (rate_Lh / v_Lh));
	for (i = 0; i <= m; i++)
	{
		fprintf_s(lin_case, "%lf\t%e\n", -2 + i * h, v[i]);
	}

	free(v);
	fclose(lin_case);
	fclose(norm);
	return 0;
    
}



