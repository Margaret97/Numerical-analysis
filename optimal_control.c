#include "pch.h"
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include<float.h> 
#define    eps1 1e-8
#define    eps2 1e-10
#define    eps3 1e-12
#define  accuracy 1e-12
double f_1(double y); // the right side of the first equation
double f_2(double t, double lambda, double alpha); // the right side of the second equation
void Runge_Kutta_8(double t, double *x, double *y, double lambda, double step, double alpha); // the Runge-Kutta method of the 8th-order
double local_error(double t, double actual_x, double actual_y, double lambda, double step, double alpha); // the local error for the Runge-Kutta method
double max(double a, double b); // the maximum function
double min(double a, double b); // the minimum function
void solution(double lambda, double alpha, double epsilon, double* difference); // we find the zero of the function x(lambda, 1) + 11/24
double f_1(double y)
{

	return  y;

}
double f_2(double t, double lambda, double alpha)
{
	if(lambda*(1 - t)*(1 + alpha * pow(t, 4)) > 1)
	{
		return 1;
	}
	if (lambda*(1 - t)*(1 + alpha * pow(t, 4)) < 1)
	{
		return -1;
	}
	else return 0;

}
void Runge_Kutta_8(double t, double *x, double *y, double lambda, double step, double alpha)
{
	double k[13], q[13];
	k[0] = step * f_1(*y);
	q[0] = step * f_2(t, lambda, alpha);
	k[1] = step * f_1(*y + (1 / 18.)*q[0]);
	q[1] = step * f_2(t + (1 / 18.)*step, lambda, alpha);
	k[2] = step * f_1(*y + (1 / 48.)*q[0] + (1 / 16.)*q[1]);
	q[2] = step * f_2(t + (1 / 12.)*step, lambda, alpha);
	k[3] = step * f_1(*y + (1 / 32.)*q[0] + (3 / 32.)*q[2]);
	q[3] = step * f_2(t + (1 / 8.)*step, lambda, alpha);
	k[4] = step * f_1(*y + (5 / 16.)*q[0] - (75 / 64.)*q[2] + (75 / 64.)*q[3]);
	q[4] = step * f_2(t + (5 / 16.)*step, lambda, alpha);
	k[5] = step * f_1(*y + (3 / 80.)*q[0] + (3 / 16.)*q[3] + (3 / 20.)*q[4]);
	q[5] = step * f_2(t + (3 / 8.)*step, lambda, alpha);
	k[6] = step * f_1(*y + (29443841 / 614563906.)*q[0] + (77736538 / 692538347.)*q[3] - (28693883 / 1125000000.)*q[4] + (23124283 / 1800000000.)*q[5]);
	q[6] = step * f_2(t + (59/400.)*step, lambda, alpha);
	k[7] = step * f_1(*y + (16016141/946692911.)*q[0] + (61564180 / 158732637.)*q[3] + (22789713 / 633445777.)*q[4] + (545815736 / 2771057229.)*q[5] - (180193667 / 1043307555.)*q[6]);
	q[7] = step * f_2(t + (93/200.)*step, lambda, alpha);
	k[8] = step * f_1(*y + (39632708 / 573591083.)*q[0] - (433636366 / 683701615.)*q[3] - (421739975 / 2616292301.)*q[4] + (100302831 / 723423059.)*q[5] + (790204164 / 839813087.)*q[6] + (800635310 / 3783071287.)*q[7]);
	q[8] = step * f_2(t + (5490023248/9719169821.)*step, lambda, alpha);
	k[9] = step * f_1(*y + (246121993/1340847787.)*q[0] - (37695042795 / 15268766246.)*q[3] - (309121744 / 1061227803.)*q[4] - (12992083 / 490766935.)*q[5] + (6005943493 / 2108947869.)*q[6] + (393006217 / 1396673457.)*q[7] + (123872331 / 1001029789.)*q[8]);
	q[9] = step * f_2(t + (13/20.)*step, lambda, alpha);
	k[10] = step * f_1(*y - (1028468189 / 846180014.)*q[0] + (8478235783 / 508512852.)*q[3] + (1311729495 / 1432422823.)*q[4] - (10304129995 / 1701304382.)*q[5] - (48777925059 / 3047939560.)*q[6] + (15336726248 / 1032824649.)*q[7] - (45442868181 / 3398467696.)*q[8] + (3065993473 / 597172653.)*q[9]);
	q[10] = step * f_2(t + (1201146811/1299019798.)*step, lambda, alpha);
	k[11] = step * f_1(*y + (185892177 / 718116043.)*q[0] - (3185094517 / 667107341.)*q[3] - (477755414 / 1098053517.)*q[4] - (703635378 / 230739211.)*q[5] + (5731566787 / 1027545527.)*q[6] + (5232866602 / 850066563.)*q[7] - (4093664535 / 808688257.)*q[8] + (3962137247 / 1805957418.)*q[9] + (65686358 / 487910083.)*q[10]);
	q[11] = step * f_2(t + step, lambda,alpha);
	k[12] = step * f_1(*y + (403863854 / 491063109.)*q[0] - (5068492393 / 434740067.)*q[3] - (411421997 / 543043805.)*q[4] + (652783627 / 914296604.)*q[5] + (11173962825 / 925320556.)*q[6] - (13158990841 / 6184727034.)*q[7] + (3936647629 / 1978049680.)*q[8] - (160528059 / 685178525.)*q[9] + (248638103 / 1413531060.)*q[10]);
	q[12] = step * f_2(t + step, lambda, alpha);
	*x = *x + (14005451/335480064.)*k[0] + (-59238493/1068277825.)*k[5] + (181606767/758867731.)*k[6] + (561292985/797845732.)*k[7] - (1041891430/1371343529.)*k[8] + (760417239/1151165299.)*k[9] + (118820643/751138087.)*k[10] - (528747749/2220607170.)*k[11] + (1/4.)*k[12];
	*y = *y + (14005451/335480064.)*q[0] + (-59238493/1068277825.)*q[5] + (181606767/758867731.)*q[6] + (561292985/797845732.)*q[7] - (1041891430/1371343529.)*q[8] + (760417239/1151165299.)*q[9] + (118820643/751138087.)*q[10] - (528747749/2220607170.)*q[11] + (1/4.)*q[12];
}

double local_error(double t, double actual_x, double actual_y, double lambda, double step, double alpha)
{
	double value_x[2], value_y[2];
	double error_x, error_y, error;
	value_x[0] = actual_x;
	value_y[0] = actual_y;
	Runge_Kutta_8(t, &value_x[0], &value_y[0], lambda, step*(1 / 2.), alpha); 
	Runge_Kutta_8(t, &value_x[0], &value_y[0], lambda, step*(1 / 2.), alpha);
	value_x[1] = actual_x;
	value_y[1] = actual_y;
	Runge_Kutta_8(t,&value_x[1], &value_y[1],lambda, step, alpha); 
	error_x = (fabs(value_x[1] - value_x[0]))*(1 / 255.);
	error_y = (fabs(value_y[1] - value_y[0]))*(1 / 255.);
	error = max(error_x, error_y);
	return error;
}

double max(double a, double b)
{
	if (a > b)
		return a;
	else return b;
}

double min(double a, double b)
{
	if (a < b)
		return a;
	else return b;
}

void solution(double lambda, double alpha, double epsilon, double* difference)
{
	FILE*t_x_y;
	double t, x, y, prev_t, current_t;
	double h, interval, err;
	int i, j;
	i = 1; j = 1;
	x = 0;
	y = 0;
	t = 0; prev_t = 0;
	h = 0.01; // the integration step
	fopen_s(&t_x_y, "C://Users//Margarita//data//t_x_y.txt", "w");
	if (t_x_y == NULL)
	{
		printf("Error opening file t_x_y.\n");
		return;
	}
	fprintf_s(t_x_y, "%lf\t%lf\t%lf\n", t, x, y);
	while (t < 1) 
	{	
			if ((t + h )<= t && (t + h) >= t)
			{
				h = 0.01;
			}
		
			if (t + h > 1)
			{
				h = 1 - t;
			}
		err = local_error(t, x, y, lambda, h, alpha);
		if (err <= epsilon)
		{
			if ((lambda*(1 - t)*(1 + alpha * pow(t, 4)) - 1)*(lambda*(1 - (t+h))*(1 + alpha * pow(t+h, 4)) - 1) < 0)
			{
				prev_t = t;
				interval = h;
				do
				{
					current_t = prev_t + interval / 1.618;

					if (prev_t==current_t)
					{
						printf("For lambda=%.16lf very small accuracy\n", lambda);
							return;
					}
					if((lambda*(1 - current_t)*(1 + alpha * pow(current_t, 4)) - 1)*(lambda*(1 -prev_t-h)*(1 + alpha * pow(prev_t+h, 4)) - 1) > 0)
					{
						interval = current_t-prev_t;

					}
					else
					{
						interval= prev_t + interval - current_t;
						prev_t = current_t;
					}
				} 
				while ((lambda*(1 - current_t)*(1 + alpha * pow(current_t, 4)) - 1) > accuracy);
				printf("For lambda=%.16lf the root is found\n", lambda);
				printf("t*=%.16e\n", current_t);
			}
			Runge_Kutta_8(t, &x, &y, lambda, 0.5*h, alpha);
			Runge_Kutta_8(t, &x, &y, lambda, 0.5*h, alpha);
			t = t + h;
			fprintf_s(t_x_y, "%.16lf\t%.16e\t%.16e\n", t, x, y);
		
			if (err < 1e-30)
			{
				h = 2 * h;
			}

			else
			{
				h = h * min(2, max((1 / 3.), 0.8*pow((epsilon / err), (1 / 9.))));
			}
		}
		else
		{
			i = 1;
			while (err > epsilon)
			{
				h = h * min(2, max((1 / 3.), 0.8*pow((epsilon / err), (1 / 9.))));
				err = local_error(t, x, y, lambda, h, alpha);
				i++;
			}
			
			if ((lambda*(1 - t)*(1 + alpha * pow(t, 4)) - 1)*(lambda*(1 - t-h)*(1 + alpha * pow(t+h, 4)) - 1) < 0)
			{
				prev_t = t;
				interval = h;
				do
				{
					current_t = prev_t + interval / 1.618;

					if (current_t == prev_t)
					{
						printf("For lambda=%.16lf very small accuracy\n", lambda);
							return;
					}
					if ((lambda*(1 - current_t)*(1 + alpha * pow(current_t, 4)) - 1)*(lambda*(1 - prev_t-h)*(1 + alpha * pow(prev_t+h, 4)) - 1) > 0)
					{
						interval = current_t-prev_t;

					}
					else
					{
						interval = prev_t + interval - current_t;
						prev_t = current_t;
					}
				} while ((lambda*(1 - current_t)*(1 + alpha * pow(current_t, 4)) - 1) > accuracy);
				printf("For lambda=%.16lf the root is found\n", lambda);
				printf("t*=%.16e\n", current_t);
			}
			Runge_Kutta_8(t,&x, &y, lambda, 0.5*h, alpha);
			Runge_Kutta_8(t, &x, &y, lambda, 0.5*h, alpha);
			t = t + h;
			fprintf_s(t_x_y, "%16.lf\t%.16e\t%.16e\n", t, x, y);
			h = h * min(2, max((1 / 3.), 0.8*pow((epsilon / err), (1 / 9.))));
		}
		j++;
	}
	*difference = x + (11 / 24.);
	printf("For lambda=%.16lf difference=%.16e\n", lambda, *difference);
	fclose(t_x_y);
	return;
}

int main(void)
{
	int i;
	double epsilon; 
	double alpha, lambda, difference, lambda_start, lambda_end, difference_end; 
	lambda_start = 0.9; // the beginning of the segment to find the zero of the function x(lambda, 1) + 11/24
	lambda_end = 1.1; // the end of the segment to find the zero of the function x(lambda, 1) + 11/24
	i = 1;
	printf("Enter epsilon (local error):\n");
	scanf_s("%lf", &epsilon);
	printf("Enter alpha:\n");
	scanf_s("%lf", &alpha);
			
	solution(lambda_end, alpha, epsilon, &difference_end);
		do
		{
			lambda = lambda_start + (lambda_end-lambda_start)/ 1.618;
			printf("Step %d\n", i);
			printf("lambda_start=%.16e\n", lambda_start);
			printf("lambda=%.16e\n", lambda);
			printf("lambda_end=%.16e\n", lambda_end);
			if (((lambda <= lambda_start)&&(lambda >= lambda_start)) || ((lambda <= lambda_end)&&(lambda >= lambda_end)))
			{
				printf("Very small accuracy for lambda\n");
				return -1;
			}
			solution(lambda, alpha, epsilon, &difference);
			if (difference*difference_end > 0)
			{
				lambda_end=lambda;
				solution(lambda_end, alpha, epsilon, &difference_end);
			}
			else
			{
				lambda_start = lambda;
		
			}
			i++;
		} while (fabs(difference) > accuracy);
		printf("difference=%.16e\n", difference);
		printf("The lambda is found\n");
		printf("lambda=%.16e\n", lambda);
		return 0;
}
	