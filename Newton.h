#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include "function.h"

const int dim = 2;
double eps = 1e-7;


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

double* Gauss(double** A, double* b, int n) // Метод Гаусса
{
	double** a = new double* [n];
	for (int i = 0; i < n; i++)
		a[i] = new double[n + 1];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			a[i][j] = A[i][j];

	for (int i = 0; i < n; i++)
		a[i][n] = b[i];

	//ShowAr(a, n, n+1);

	double max;
	double m;
	int jmax;
	bool flag, flag2;
	flag2 = true;
	for (int i = 0; i < n; i++)
	{
		if (flag2)
		{
			max = fabs(a[i][i]);

			for (int j = i + 1; j < n; j++) //Выбор главного элемента
			{
				flag = false;
				if (fabs(a[j][i]) > fabs(max))
				{
					max = fabs(a[j][i]);
					jmax = j;
					flag = true;
				}

				if (flag == true)
				{
					double* a1 = new double[n + 1];
					for (int k = 0; k < n + 1; k++)//перемена строк местами
					{
						a1[k] = a[i][k];
						a[i][k] = a[jmax][k];
						a[jmax][k] = a1[k];
					}
					delete a1;
				}
			}


			for (int j = i + 1; j < n; j++)
			{

				m = a[j][i] / a[i][i];
				for (int k = 0; k < n + 1; k++)
				{
					a[j][k] = a[j][k] - m * a[i][k];
				}
			}

		}

		if (fabs(a[i][i]) < eps)
		{
			flag2 = false;
		}
	}

	double* x = new double[n];
	double sum;
	if (flag2)
	{
		for (int i = n - 1; i >= 0; i--)
		{

			{
				sum = 0;
				for (int j = n - 1; j > i; j--)
					sum += x[j] * a[i][j];
				x[i] = (a[i][n] - sum) / a[i][i];
			}
		}

		/*cout << "Solution vector" << endl;
		for (int i = 0; i<n; i++)
		cout << x[i] << endl;*/
	}
	else
	{
		for (int i = 0; i < n; i++)
			x[i] = 0;
		//	cout << "Not consistent" << endl;
	}

	return x;

}

double* Newton_System (double* temp, double tau, double t0, double T)
{
	double** J;
	double* Y = new double[dim];
	double* F = new double[dim]; // F' Y = F
	
	double* x = new double[dim];
	double* tempx = new double[dim];//!!!!!!!!!!!!!!!!!!
	int iter = 0;
	x[0] = t0; x[1] = T;

			do {
				++iter;
				
				J = Jacobi(x, tau, temp);
				std::cout << "Jucobi matrix iter = " << iter << std::endl;
				Matrix_Print(J);

				for (size_t i = 0; i < dim; ++i)
					F[i] = -Function_system(i, tau, x, temp);

				std::cout << "F " << std::endl;
				Vector_Print(F);

				Y = Gauss(J, F, 2);
				std::cout << "Rule = " << Y[0] * Y[0] + Y[0] * Y[0] << std::endl;

				if (Y[0] * Y[0] + Y[0] * Y[0] != 0) {

					tempx[0] = x[0];
					tempx[1] = x[1];

					x[0] = tempx[0] + Y[0];
					x[1] = tempx[1] + Y[1];

				}
				else break; 
			} while (norma(x, tempx) > eps);

			return x;
}

//double* gauss(double** pA, double* b, int n)
//{
//	double* x, max;
//	int k, index;
//	const double eps = 0.00001; // точность 
//	x = new double[n];
//	k = 0;
//	while (k < n)
//	{
//		// Поиск строки с максимальным a[i][k] 
//		max = abs(pA[k][k]);
//		index = k;
//		for (int i = k + 1; i < n; i++)
//		{
//			if (abs(pA[i][k]) > max)
//			{
//				max = abs(pA[i][k]);
//				index = i;
//			}
//		}
//		// Перестановка строк 
//		if (max < eps)
//		{
//			// нет ненулевых диагональных элементов 
//			/*std::cout « "resh net ";
//			std::cout « index « " matrix A";*/
//			return 0;
//		}
//		for (int j = 0; j < n; j++)
//		{
//			double temp = pA[k][j];
//			pA[k][j] = pA[index][j];
//			pA[index][j] = temp;
//		}
//		double temp = b[k];
//		b[k] = b[index];
//		b[index] = temp;
//		// Нормализация уравнений 
//		for (int i = k; i < n; i++)
//		{
//			double temp = pA[i][k];
//			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить 
//			for (int j = 0; j < n; j++)
//				pA[i][j] = pA[i][j] / temp;
//			b[i] = b[i] / temp;
//			if (i == k) continue; // уравнение не вычитать само из себя 
//			for (int j = 0; j < n; j++)
//				pA[i][j] = pA[i][j] - pA[k][j];
//			b[i] = b[i] - b[k];
//		}
//		k++;
//	}
//	for (k = n - 1; k >= 0; k--)
//	{
//		x[k] = b[k];
//		for (int i = 0; i < k; i++)
//			b[i] = b[i] - pA[i][k] * x[k];
//	}
//	return x;
//}