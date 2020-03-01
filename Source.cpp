#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>
#include "Newton.h"

//#include "Header1.h"//test

double* Runge_Kutta_4 (double t0, double T, double tau);
double* Evident_Euler (double t0, double T, double tau);
double* Runge_Kutta_2 (double t0, double T, double tau);
double* Not_Evident_Euler (double t0, double T, double tau);
double* Symmetrical(double t0, double T, double tau);

int main(int argc, char* argv[])
{
	double t0, T, tau;
	t0 = 0.; T = 1.; tau = 0.01;

	for (size_t i = 0; i < 2; ++i)
		std::cout << Evident_Euler(t0, T, tau)[i] << std::endl;

	for (size_t i = 0; i < 2; ++i)
		std::cout << Runge_Kutta_2(t0, T, tau)[i] << std::endl;

	for (size_t i = 0; i < 2; ++i)
		std::cout << Not_Evident_Euler(t0, T, tau)[i] << std::endl;

	for (size_t i = 0; i < 2; ++i)
		std::cout << Symmetrical(t0, T, tau)[i] << std::endl;


	return 0;
}


double* Evident_Euler(double t0, double T, double tau)
{
	double* y;
	double* temp;
	for (size_t i = 0; i < 2; ++i) {
		temp = new double[2];
		y = new double[2];
	}

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
	double* y;
	double* temp;
	for (size_t i = 0; i < 2; ++i) 
	{
		temp = new double[2];
		y = new double[2];
	}

	y[0] = 1.; y[1] = 0.; //начальное условие
	int N = round((T - t0) / tau) + 1;

	std::ofstream file;
	file.open("2.dat");
	if (file.is_open())
	{
	for(size_t i=0; i<N; ++i)
	{
		//y = NewtonSys(t0, T, tau, y, 1);
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
	double* y;
	double* temp;
	double* k1;
	double* k2;

	for (size_t i = 0; i < 2; ++i) {
		y = new double[2];
		temp = new double[2];
		k1 = new double[2];
		k2 = new double[2];
	}

	y[0] = 1.; y[1] = 0.;//начальное условие

	int N = (T - t0) / tau + 1;

	std::ofstream file;
	file.open("RK_2.dat");

	if (file.is_open()) {
		for (size_t i = 0; i < N; ++i) {

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
	double* y;
	double* temp;

	for (size_t i = 0; i < 2; ++i)
	{
		temp = new double[2];
		y = new double[2];
	}

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

double* Runge_Kutta_4(double t0, double T, double tau)
{
	double* y;
	double* temp;
	double* k1;
	double* k2;
	double* k3;
	double* k4;
	for (size_t i = 0; i < 2; ++i) {
		y = new double[2];
		temp = new double[2];
		k1 = new double[2];
		k2 = new double[2];
		k3 = new double[2];
		k4 = new double[2];
	}

	y[0] = 0.; y[1] = 1.;//начальное условие

	int N = (T - t0) / tau + 1;

	for (size_t i = 0; i < N; ++i) {

		for (size_t j = 0; i < 2; ++j)
			temp[j] = y[j];

		for (size_t j = 0; i < 2; ++j)
			k1[j] = F(j, i * t0, y);

		for (size_t j = 0; i < 2; ++j)
			temp[j] = y[j] + tau / 2 * k1[j];

		for (size_t j = 0; i < 2; ++j)
			k2[j] = F(i, i * t0 + tau / 2, temp);

		for (size_t j = 0; i < 2; ++j)
			temp[j] = y[j] + tau / 2 * k2[j];

		for (size_t j = 0; i < 2; ++j)
			k3[j] = F(i, i * t0 + tau / 2, temp);

		for (size_t j = 0; i < 2; ++j)
			temp[j] = y[j] + tau * k3[j];

		for (size_t j = 0; i < 2; ++j)
			k4[j] = F(i, i * t0 + tau, temp);

		for (size_t j = 0; i < 2; ++j)
			y[j] = y[j] + tau / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);

	}

	return y;
}