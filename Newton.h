#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include "function.h"

double D(int i, int j, double* y, double tau, double* temp)
{   // i - номер переменной по которой дефференцируем
	// i - номер уравнения системы
	double eps = 1e-7;
	double* y_right = new double[dim];
	for (size_t k = 0; k < dim; ++k)
	{
		y_right[k] = y[k];

		if (k == j)
			y_right[k] = y[k] + eps;
	}

	return ((Function_system(i, tau, y_right, temp) - Function_system(i, tau, y, temp)) / eps);
}

double** Jacobi(double* y, double tau, double* temp)
{
	double** J = new double* [dim];
	for (size_t i = 0; i < dim; ++i)
		J[i] = new double[dim];

	for (size_t i = 0; i < dim; ++i) 
	{
		for (size_t j = 0; j < dim; ++j)
			J[i][j] = D(i, j, y, tau, temp);
	}

	return J;
}

void Obrat_hod(double** A, double* B, double* x)
{
	double Sum = 0;
	for (int i = dim - 1; i >= 0; i--)
	{
		for (int j = i + 1; j < dim; j++)
			Sum = Sum + A[i][j] * x[j];

		x[i] = (B[i] - Sum) / (A[i][i]);
		Sum = 0;
	}
}

void Gauss(double** A, double* B, double* X)
{
	double max_matrix_element;
	int k = 0;
	double g;
	int imax = 0; //индекс максимального элемента матрицы
	while (k < dim)
	{
		max_matrix_element = fabs(A[k][k]);
		imax = k;

		for (size_t i = k + 1; i < dim; ++i)
			if (fabs(A[i][k]) > max_matrix_element)
			{
				max_matrix_element = fabs(A[i][k]);
				imax = i;
			}

		if (max_matrix_element < eps)
		{
			dim = 0;
			std::cout << "Degenerate matrix" << std::endl;
			return;
		}

		for (size_t j = 0; j < dim; ++j)
			std::swap(A[k][j], A[imax][j]);

		std::swap(B[k], B[imax]);

		for (size_t j = k + 1; j < dim; ++j)
		{
			g = A[j][k] / A[k][k];

			for (size_t i = k + 1; i < dim; i++)
				A[j][i] = A[j][i] - g * A[k][i];
			B[j] = B[j] - g * B[k];
			A[j][k] = 0;
		}

		k++;
	}
	Obrat_hod(A, B, X);
}


double* Newton_System (double* temp, double tau, double t0, double T, int method)
{
	double** J;
	double* Y = new double[dim];
	double* F = new double[dim]; // F' Y = F
	
	double* x = new double[dim];
	double* temp_x = new double[dim];//!!!!!!!!!!!!!!!!!!

	x[0] = t0;
	x[1] = T;

			do {

				J = Jacobi(x, tau, temp);

				if (J[0][0] * J[1][1] - J[1][0] * J[0][1] < eps)//проверка вырожденности для матрицы 2 на 2
					break;
				/*std::cout << "Jucobi matrix iter = " << iter << std::endl;
				Print_Matrix(J);*/

				if(method == 0) //for Not_Evident_Euler
				for (size_t i = 0; i < dim; ++i)
					F[i] = -Function_system(i, tau, x, temp);


				if (method == 1)
					for (size_t i = 0; i < dim; ++i)
						F[i] = -Function_system_for_symmetrical(i, tau, x, temp);

				/*std::cout << "F " << std::endl;
				Print_Vector(F);*/

				Gauss(J, F, Y);

				if ((Y[0]!= 0) || (Y[1] != 0))
				{
					temp_x[0] = x[0];
					temp_x[1] = x[1];

					x[0] = temp_x[0] + Y[0];
					x[1] = temp_x[1] + Y[1];
				}
				else break; 
			} while (norma(x, temp_x) > eps);

			delete[] temp_x,F,Y;
			Delete_Matrix(J);
			return x;
}
