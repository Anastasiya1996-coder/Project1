#pragma once
int dim = 2;
double eps = 1e-7;
int nomertesta = 1;
const double k = 20; // Коэф. жесткости.
const double m = 0.3;

double F(int i, double t, double* y) // Systema yravnenii razmera razm, i - nomer yravneniya v systeme
{
	double F;
	if (i == 0) // Pervoe yravnenie systemi
	{
		switch (nomertesta) {
		case 1:
			F = y[1];
			return F;
		case 2:
			F = 2 * y[0] + y[1] * y[1] - 1;
			return F;
		case 3:
			F = 1 - y[0] * y[0] - y[1] * y[1];
			return F;
		case 4:
			F = (0.4 + 0.3 * sin(0.5 * t)) * y[1] + 0.1 * y[0] - y[0] * y[0] * y[0] - y[0] * y[1] * y[1];
			return F;
		case 5:
			F = 10 * (y[1] - y[0]);
			return F;
		case 6:
			F = y[1];
			return F;
		default:
			std::cout << "Invalid number!";
			F = 1;
			return F;
		}
	}
	else if (i == 1) // Vtoroe yravnenie systemi
	{
		switch (nomertesta) {
		case 1:
			F = -k / m * y[0];
			return F;
		case 2:
			F = 6 * y[0] - y[1] * y[1] + 1;
			return F;
		case 3:
			F = 2 * y[0];
			return F;
		case 4:
			F = -(0.4 + 0.3 * sin(0.5 * t)) * y[0] + 0.2 * y[1] - y[1] * y[1] * y[1] - y[1] * y[0] * y[0];
			return F;
		case 5:
			F = y[0] * (28 - y[2]) - y[1];
			return F;
		case 6:
			F = -4 * y[1] - 5 * y[0];
			return F;
		default:
			std::cout << "Nevernii nomer testa!!!!!";
			F = 1;
			return F;
		}
	}
	else if (i == 2)
	{
		F = y[0] * y[1] - 8 / 3 * y[2];
		return F;
		//	cout << "Ne zadana systema pod tekyshii razmer = " << razm << endl;
	}

}

double Function_system(size_t i, double tau, double* y, double* temp)
{
	double Y = y[i] - temp[i] - tau * F(i, tau, y); //задача для маятника, t_n не учитывается функцией
	return Y;
}

double Function_system_for_symmetrical(size_t i, double tau, double* y, double* temp)
{
	double Y = y[i] - temp[i] - tau/2 *(F(i, tau, temp) + F(i, tau, y)); //задача для маятника, t_n не учитывается функцией
	return Y;
}

void Print_Vector(double* A)
{
	for (int i = 0; i < dim; i++)
		std::cout << A[i] << " ";
	std::cout << std::endl;
}

void Print_Matrix(double** A)
{
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
			std::cout << A[i][j] << " ";
		std::cout << std::endl;
	}
}

void Delete_Matrix(double** A)
{
	for (size_t i = 0; i < dim; ++i)
		delete[] A[i];
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

