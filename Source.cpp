#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>
#include "Newton.h"

//#include "Header1.h"//test

double tt = 1;

double EPSILON = 1e-7;

double* Runge_Kutta_4 (double t0, double T, double tau, bool record);
double* Evident_Euler (double t0, double T, double tau);
double* Runge_Kutta_2 (double t0, double T, double tau);
double* Not_Evident_Euler (double t0, double T, double tau);
double* Symmetrical(double t0, double T, double tau);
double* Evident_Adams(double t0, double T, double tau);
double* Predictor_Corrector(double t0, double T, double tau);
double* Runge_Rule(double t0, double T, double tau);

int main(int argc, char* argv[])
{
	double t0, T, tau, y_0;
	t0 = 0.; T = 1.; tau = 0.02;

	Evident_Euler(t0, T, tau);
	//std::cout << "Precision of solution E-2 = " << fabs(Evident_Euler(t0, T, 0.0006)[0] - Reference_solution(T)) << std::endl;// tau/15
	//std::cout << "Precision of solution E-4 = " << fabs(Evident_Euler(t0, T, tau/1000)[0] - Reference_solution(T)) << std::endl;

	Not_Evident_Euler(t0, T, tau);
	/*std::cout << "Precision of solution E-2 = " << fabs(Not_Evident_Euler(t0, T, tau/15)[0] - Reference_solution(T)) << std::endl;
	std::cout << "Precision of solution E-4 = " << fabs(Not_Evident_Euler(t0, T, tau/110)[0] - Reference_solution(T)) << std::endl;
	std::cout << "Precision of solution E-6 = " << fabs(Not_Evident_Euler(t0, T, tau / 10000)[0] - Reference_solution(T)) << std::endl;*/

	Runge_Kutta_2(t0, T, tau);
	//std::cout << "Precision of solution E-2 = " << fabs(Runge_Kutta_2(t0, T, tau/10)[0] - Reference_solution(T)) << std::endl;
	//std::cout << "Precision of solution E-4 = " << fabs(Runge_Kutta_2(t0, T, tau / 100)[0] - Reference_solution(T)) << std::endl;
	//std::cout << "Precision of solution E-7 = " << fabs(Runge_Kutta_2(t0, T, tau / 110)[0] - Reference_solution(T)) << std::endl;

	Symmetrical(t0, T, tau);
	/*std::cout << "Precision of solution E-2 = " << fabs(Symmetrical(t0, T, tau / 5)[0] - Reference_solution(T)) << std::endl;
	std::cout << "Precision of solution E-4 = " << fabs(Symmetrical(t0, T, tau / 100)[0] - Reference_solution(T)) << std::endl;
	std::cout << "Precision of solution E-7 = " << fabs(Symmetrical(t0, T, tau / 1000)[0] - Reference_solution(T)) << std::endl;*/

	for (size_t i = 0; i < 2; ++i) {
	std::cout << Runge_Kutta_4(t0, T, tau, 1)[i] << std::endl;}

	//std::cout << "Precision of solution = " << fabs(Runge_Kutta_4(t0, T, 1 * tau,1)[0] - Reference_solution(T)) << std::endl; 

	

	for (size_t i = 0; i < 2; ++i)
		std::cout << Evident_Adams(t0, T, tau)[i] << std::endl;

	for (size_t i = 0; i < 2; ++i)
		std::cout << Predictor_Corrector(t0, T, tau)[i] << std::endl;

	for (size_t i = 0; i < 2; ++i)
		std::cout << Runge_Rule(t0, T, tau)[i] << std::endl;

	return 0;
}


double* Evident_Euler(double t0, double T, double tau)
{
	double* y = new double[2];
	double* temp = new double[2];
	
	y[0] = 1.; y[1] = 0.;//начальное условие
	int N = round((T - t0) / tau) + 1;

	std::ofstream file;
	file.open("1.dat");
	if (file.is_open())
	{
		for (size_t i = 0; i < N; ++i) //i отвечает за шаг
		{
			for (size_t j = 0; j < 2; ++j)
				temp[j] = y[j];

			for (size_t j = 0; j < 2; ++j)
				y[j] = temp[j] + tau * F(j, i * tau, temp);

			file << y[0] << " " << y[1] << std::endl;
		}
	}
	file.close();
	delete[] temp;
	return y;
}

double* Not_Evident_Euler(double t0, double T, double tau)
{
	double* y = new double[2];
	double* temp = new double[2];

	y[0] = 1.; y[1] = 0.; //начальное условие
	int N = round((T - t0) / tau) + 1;

	std::ofstream file;
	file.open("2.dat");
	if (file.is_open())
	{
	  for(size_t i=0; i<N; ++i)
	  {
		y = Newton_System(y, tau, t0, T, 0);
		file << y[0] << " " << y[1] << std::endl;
	  }	
	}
	file.close();
	delete[] temp;
	return y;
}

double* Runge_Kutta_2(double t0, double T, double tau)
{
	double* y = new double[dim];
	double* temp = new double[dim];
	double* k1 = new double[dim];
	double* k2 = new double[dim];

	y[0] = 1.; y[1] = 0.;//начальное условие

	int N = (T - t0) / tau + 1;

	std::ofstream file;
	file.open("RK_2.dat");

	if (file.is_open()) {
		for (size_t i = 0; i < N; ++i) 
		{
			for (size_t j = 0; j < 2; ++j)
				temp[j] = y[j];

			for (size_t j = 0; j < 2; ++j)
				k1[j] = F(j, i * t0, temp);

			for (size_t j = 0; j < 2; ++j)
				y[j] = temp[j] + tau * k1[j];

			for (size_t j = 0; j < 2; ++j)
				k2[j] = F(j, i * t0 + tau, y);

			for (size_t j = 0; j < 2; ++j)
				y[j] = temp[j] + tau * (0.5 * k1[j] + 0.5 * k2[j]);

			file << y[0] << " " << y[1] << std::endl;
		}
	}
	file.close();

	return y;
}

double* Symmetrical(double t0, double T, double tau)
{
	double* y = new double[dim];
	double* temp = new double[dim];

	y[0] = 1.; y[1] = 0.; //начальное условие
	int N = round((T - t0) / tau) + 1;

	std::ofstream file;
	file.open("S.dat");
	if (file.is_open())
	{
		for (size_t i = 0; i < N; ++i)
		{
			y = Newton_System(y, tau, t0, T, 1);
			file << y[0] << " " << y[1] << std::endl;
		}
	}
	file.close();
	delete[] temp;
	return y;
}

double* Runge_Kutta_4(double t0, double T, double tau, bool record)
{
	double* y = new double[dim];

	double* temp = new double[dim];
	double* k1 = new double[dim];
	double* k2 = new double[dim];
	double* k3 = new double[dim];
	double* k4 = new double[dim];

	y[0] = 1.; y[1] = 0.;//начальное условие

	int N = round((T - t0) / tau) + 1;

	std::ofstream file;
	if(record)
	    file.open("RK_4.dat");
	else
		file.open("RK_4_copy.dat");

	if (file.is_open())
	{
	   for (size_t i = 0; i < N; ++i)
	   {
		for (size_t j = 0; j < dim; ++j)
			temp[j] = y[j];

		for (size_t j = 0; j < dim; ++j)
			k1[j] = F(j, i * t0, temp);

		for (size_t j = 0; j < dim; ++j)
			 y[j] = temp[j] + tau / 2 * k1[j];

		for (size_t j = 0; j < dim; ++j)
			k2[j] = F(j, i * t0 + tau / 2, y);

		for (size_t j = 0; j < dim; ++j)
			 y[j] = temp[j] + tau / 2 * k2[j];

		for (size_t j = 0; j < dim; ++j)
			k3[j] = F(j, i * t0 + tau / 2, y);

		for (size_t j = 0; j < dim; ++j)
			 y[j] = temp[j] + tau * k3[j];

		for (size_t j = 0; j < dim; ++j)
			k4[j] = F(j, i * t0 + tau, y);

		for (size_t j = 0; j < dim; ++j)
			y[j] = temp[j] + tau / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);

			file << y[0] << " " << y[1] << std::endl;

	   if (fabs(i*tau-0.02)<1e-15)
	   {
		   std::cout << "err rk4 5tau " << fabs(y[0] - Reference_solution(i*tau))<<"\n";

	   }
	   }



	}
	file.close();
	return y;
}

double* Evident_Adams(double t0, double T, double tau)
{
	double* y = new double[dim];
	double* temp = new double[dim];
	double* y_0 = new double[dim];//y_{n}
	double* y_1 = new double[dim];//y_{n-1} ...
	double* y_2 = new double[dim];
	double* y_3 = new double[dim];

	y[0] = 1.; y[1] = 0.;//начальное условие

	int N = round((T - t0) / tau) + 1;

	std::ofstream file;
	file.open("A.dat");

	if (file.is_open())
	{ 
		y_3[0] = 1.; y_3[1] = 0.;
		y_2 = Runge_Kutta_4(0, tau, tau, 0);
		y_1 = Runge_Kutta_4(0, 2 * tau, tau, 0);
		y_0 = Runge_Kutta_4(0, 3 * tau, tau, 0);

		file << y_3[0] << " " << y_3[1] << std::endl;
		file << y_2[0] << " " << y_2[1] << std::endl;
		file << y_1[0] << " " << y_1[1] << std::endl;
		file << y_0[0] << " " << y_0[1] << std::endl;

		for (size_t i = 4; i < N; ++i)
		{
			for (size_t j = 0; j < dim; ++j)
				y[j] = y_0[j] + (tau / 24) * (55 * F(j, tau, y_0) - 59 * F(j, tau, y_1) + 37 * F(j, tau, y_2) - 9 * F(j, tau, y_3));

			for (size_t k = 0; k < dim; ++k)
			{
				y_3[k] = y_2[k];
				y_2[k] = y_1[k];
				y_1[k] = y_0[k];
				y_0[k] = y[k];
			}
			
			file << y[0] << " " << y[1] << std::endl;
		}
				
	}
	file.close();
	return y;
}

double* Predictor_Corrector(double t0, double T, double tau)
{
	double* y = new double[dim];
	double* y_0 = new double[dim];//y_{n}
	double* y_1 = new double[dim];//y_{n-1} ...
	double* y_2 = new double[dim];
	double* y_3 = new double[dim];

	y[0] = 1.; y[1] = 0.;//начальное условие

	int N = round((T - t0) / tau) + 1;

	std::ofstream file;
	file.open("PC.dat");

	if (file.is_open())
	{
		y_3[0] = 1.; y_3[1] = 0.;
		y_2 = Runge_Kutta_4(0, tau, tau, 0);
		y_1 = Runge_Kutta_4(0, 2 * tau, tau, 0);
		y_0 = Runge_Kutta_4(0, 3 * tau, tau, 0);

		file << y_3[0] << " " << y_3[1] << std::endl;
		file << y_2[0] << " " << y_2[1] << std::endl;
		file << y_1[0] << " " << y_1[1] << std::endl;
		file << y_0[0] << " " << y_0[1] << std::endl;

		for (size_t i = 4; i < N; ++i)
		{
			for (size_t j = 0; j < dim; ++j)
			{
				y[j] = y_0[j] + (tau / 24) * (55 * F(j, tau, y_0) - 59 * F(j, tau, y_1) + 37 * F(j, tau, y_2) - 9 * F(j, tau, y_3));
				y[j] = y_0[j] + (tau / 24) * (9 * F(j, tau, y) + 19 * F(j, tau, y_0) - 5 * F(j, tau, y_1) + F(j, tau, y_2));
			}

			for (size_t k = 0; k < dim; ++k)
			{
				y_3[k] = y_2[k];
				y_2[k] = y_1[k];
				y_1[k] = y_0[k];
				y_0[k] = y[k];
			}

			file << y[0] << " " << y[1] << std::endl;
		}

	}
	file.close();
	return y;

}

//double* Runge_Rule(double t0, double T, double tau)
//{
//	double* y = new double[dim];
//	double* y1_2 = new double[dim];
//
//	y[0] = 1.; y[1] = 0.;//начальное условие
//
//	file << y[0] << " " << y[1] << std::endl;
//
//	int N = round((T - t0) / tau) + 1;
//
//	std::ofstream file;
//	file.open("RR.dat");
//
//	if (file.is_open())
//	{
//		for (size_t i = 1; i < N; ++i)
//		{
//			for (size_t j = 0; j < dim; ++j) 
//			{
//				y[j] = Runge_Kutta_4(t0, i*tau, tau,0);
//			    y1_2[j] = Runge_Kutta_4(t0, i * tau, tau/2, 0);
//			}
//			
//			if (norma(y, y1_2) / 15 < eps)
//			{
//				for (size_t j = 0; j < dim; ++i)
//					y[j] = y1_2[j];
//
//				file << y[0] << " " << y[1] << std::endl;
//			}
//			else 
//			{
//				tau = tau / 2;
//				--i;
//			}
//			//file << y[0] << " " << y[1] << std::endl;
//
//		}
//	}
//	file.close();
//	return y;
//
//}

double* Runge_Rule(double t0, double T, double tau)
{
	double* y = new double[dim];
	double* y1_2 = new double[dim];
	double* y_temp = new double[dim];
	double* temp = new double[dim];
	double* k1 = new double[dim];
	double* k2 = new double[dim];
	double* k3 = new double[dim];
	double* k4 = new double[dim];

	y[0] = 1.; y[1] = 0.;   //начальное условие
	y1_2[0] = 1.; y1_2[1] = 0.;
	double length = 0. ;

	int i = 0;
	int N = round((T - t0) / tau) + 1;

	std::ofstream err;
	err.open("err.dat");

	std::ofstream file;
	file.open("RR.dat");

	if (file.is_open())
	{
		while (length < T)
		{
				for (size_t j = 0; j < dim; ++j)
					temp[j] = y[j];

				for (size_t j = 0; j < dim; ++j)
					k1[j] = F(j, i * t0, temp);

				for (size_t j = 0; j < dim; ++j)
					y[j] = temp[j] + tau / 2 * k1[j];

				for (size_t j = 0; j < dim; ++j)
					k2[j] = F(j, i * t0 + tau / 2, y);

				for (size_t j = 0; j < dim; ++j)
					y[j] = temp[j] + tau / 2 * k2[j];

				for (size_t j = 0; j < dim; ++j)
					k3[j] = F(j, i * t0 + tau / 2, y);

				for (size_t j = 0; j < dim; ++j)
					y[j] = temp[j] + tau * k3[j];

				for (size_t j = 0; j < dim; ++j)
					k4[j] = F(j, i * t0 + tau, y);

				for (size_t j = 0; j < dim; ++j)
					y[j] = temp[j] + tau / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);

				std::cout << "y: ";
				Print_Vector(y);

			tau = tau / 2; 

			for (size_t k = 0; k < 2; ++k) 
			{
				for (size_t j = 0; j < dim; ++j)
					temp[j] = y1_2[j];

				for (size_t j = 0; j < dim; ++j)
					k1[j] = F(j, i * t0, temp);

				for (size_t j = 0; j < dim; ++j)
					y1_2[j] = temp[j] + tau / 2 * k1[j];

				for (size_t j = 0; j < dim; ++j)
					k2[j] = F(j, i * t0 + tau / 2, y1_2);

				for (size_t j = 0; j < dim; ++j)
					y1_2[j] = temp[j] + tau / 2 * k2[j];

				for (size_t j = 0; j < dim; ++j)
					k3[j] = F(j, i * t0 + tau / 2, y1_2);

				for (size_t j = 0; j < dim; ++j)
					y1_2[j] = temp[j] + tau * k3[j];

				for (size_t j = 0; j < dim; ++j)
					k4[j] = F(j, i * t0 + tau, y1_2);

				for (size_t j = 0; j < dim; ++j)
					y1_2[j] = temp[j] + tau / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
			}
		
			std::cout << "y1_2: ";
				Print_Vector(y1_2);

				if (err.is_open())
				{
					err << i << " " << (norma(y, y1_2) / 15) << std::endl;
				}


			if ((norma(y, y1_2) / 15) <= EPSILON)
			{
				std::cout << "if" << std::endl;

				/*if ((norma(y, y1_2) / 15) <= (1000 * EPSILON))
				{
					tau = tau * 2;
					std::cout << "if-if";
				}*/
				
				for (size_t j = 0; j < 2; j++)
					y[j] = y1_2[j];
				
				tau = tau * 2;

				length = length + tau;
				std::cout <<"length = "<< length << std::endl;

			}
			else
			{

				--i;
				std::cout << "else" << std::endl;
			}

			i++;
			
			file << y[0] << " " << y[1] << std::endl;
		}
	}
	file.close();
	err.close();
	return y;
}