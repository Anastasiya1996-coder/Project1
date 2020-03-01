#pragma once
#include <iostream>
#include <fstream>
#include<cmath>
#include "function.h"

const int dim = 2;
double eps = 1e-7;

double** Wnet(double a0, double b0, int N) // Создание сетки
{
	double** W = new double* [2];
	for (int i = 0; i < 2; i++)
		W[i] = new double[N];

	for (int i = 0; i < N; i++)
	{
		W[0][i] = -a0 + i * 2 * a0 / N;  // X1
		W[1][i] = -b0 + i * 2 * b0 / N;  // X2
		//std::cout << i  << " " << W[0][i] << " " << W[1][i] << std::endl;
		//std::cout << i << " " << W[0][i] << " " << W[1][i] << std::endl;
	}
	return W;
}

void Matrix_Print(double** A)
{
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
			std::cout << A[i][j] << " ";
		std::cout << std::endl;
	}
}

void Vector_Print(double* A)
{
	for (int i = 0; i < dim; i++)
			std::cout << A[i] << " ";
		std::cout << std::endl;
}

double norma(double* y, double* temp)
{
	double norma = fabs(y[0] - temp[0]);

	for (size_t i = 1; i < dim; ++i)
	{
		if (fabs(y[i] - temp[i]) > norma)
			norma = fabs(y[i] - temp[i]);
	}
	return norma;
}

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

double* NewtonSys(double a0, double b0, double tau, double* ty, int method)
{
	double* Resh = new double[dim];
	int N = 2; // КОЛИЧЕСТВО ЯЧЕЕК
	double** W, ** J;
	double* x = new double[dim];
	double* tempx = new double[dim];

	double tempiter, ii;
	int iter = 0;
	bool stop;
	double norm1;
	double* Y = new double[2];

	W = Wnet(a0, b0, N); // Создали сетку

	for (int i = 0; i < N + 1; i++)
		for (int j = 0; j < N + 1; j++)
		{
			x[0] = W[0][i];
			x[1] = W[1][j];
			std::cout << "i = " << i << "  j = " << j << std::endl;
			tempiter = iter; // Для вычисления итераций на шаге
			stop = true;
			while (stop)
			{
				J = Jacobi(x, tau, ty);	// F'(xk)

				std::cout << "Jucobi matrix iter = " << iter << std::endl;
				Matrix_Print(J);

				Y[0] = -Function_system(0, tau, x, ty); Y[1] = -Function_system(1, tau, x, ty); // F(xk)

				std::cout << "Y_begin " << std::endl;
				Vector_Print(Y);

				Y = Gauss(J, Y, 2);		
				std::cout << "Rule = " << Y[0] * Y[0] + Y[0] * Y[0] << std::endl;
							
				if (Y[0] * Y[0] + Y[0] * Y[0] != 0)		// Проверка, чтобы если вырожденная F'(xk), то ничего не делал
				{
					iter++;
					//	cout << "Y1 = " << Y[0] << "   Y2 = " << Y[1] << endl;
					tempx[0] = x[0];
					tempx[1] = x[1];
					x[0] = x[0] + Y[0]; // САМ МЕТОД НЬЮТОНА
					x[1] = x[1] + Y[1];
					ii = iter - tempiter;
					if (fabs(x[0] - tempx[0]) <= fabs(x[1] - tempx[1]))
					{
						norm1 = fabs(x[1] - tempx[1]);
					}
					else { norm1 = fabs(x[0] - tempx[0]); }
					if (norm1 < eps)
					{
						stop = false;
						//std::cout << "x[0] = " << x[0] << "   x[1] = " << x[1] << std::endl;
					}
					if ((iter - tempiter) >= 30) {
						stop = false; //cout << "Iter 30+" << endl; ii = 100;
					};
				}
				else { stop = false, ii = 100; };
			}
			//Запись функции для рисования в wolfram

		//	fout << W[0][i] << " " << W[1][j] << " " << ii << std:: endl;

		//	cout << "[ " << W[0][i] << " ; " << W[1][j] << " ]" << "   ITER = " << ii << "   X,Y = " << x[0] << "  " << x[1] << endl;
			if (ii < 30)
			{
				Resh[0] = x[0], Resh[1] = x[1];

				std::cout << "[ " << W[0][i] << " ; " << W[1][j] << " ]" << "   ITER = " << ii << "   X,Y = " << x[0] << "  " << x[1] << std::endl;
				j = N + 2;
				i = N + 2;
			}
		}
	//fout.close();

	std::cout << "iter " << iter << std::endl;
	return Resh;
}


double* gauss(double** pA, double* b, int n)
{
	double* x, max;
	int k, index;
	const double eps = 0.00001; // точность 
	x = new double[n];
	k = 0;
	while (k < n)
	{
		// Поиск строки с максимальным a[i][k] 
		max = abs(pA[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(pA[i][k]) > max)
			{
				max = abs(pA[i][k]);
				index = i;
			}
		}
		// Перестановка строк 
		if (max < eps)
		{
			// нет ненулевых диагональных элементов 
			/*std::cout « "resh net ";
			std::cout « index « " matrix A";*/
			return 0;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = pA[k][j];
			pA[k][j] = pA[index][j];
			pA[index][j] = temp;
		}
		double temp = b[k];
		b[k] = b[index];
		b[index] = temp;
		// Нормализация уравнений 
		for (int i = k; i < n; i++)
		{
			double temp = pA[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить 
			for (int j = 0; j < n; j++)
				pA[i][j] = pA[i][j] / temp;
			b[i] = b[i] / temp;
			if (i == k) continue; // уравнение не вычитать само из себя 
			for (int j = 0; j < n; j++)
				pA[i][j] = pA[i][j] - pA[k][j];
			b[i] = b[i] - b[k];
		}
		k++;
	}
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = b[k];
		for (int i = 0; i < k; i++)
			b[i] = b[i] - pA[i][k] * x[k];
	}
	return x;
}