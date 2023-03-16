///Program to solve system of linear equations using Cholesky Decomposition 
///				Kabir Lovero
#include <iostream>
#include <vector>


using namespace std;


vector<vector<double>> CholeskyDecomposition(vector<vector<double>> A, int n)
{
	
	n = A.size();
	//initialize lower diagonal array with zeros
	vector<vector<double>> L(n,vector<double>(n));
	
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j <= i ; j++)
		{
			double sum = 0.;
			if (i == j)
			{
				for (int k = 0; k < j; k++)
				{
			
					sum += pow(abs(L[j][k]), 2);
				}
				L[j][j] = sqrt(abs(A[j][j] - sum));
			}
			else
			{ 
			/// get value of L[i][j] 
				for (int k = 0; k < j; k++)
				{

					sum += L[i][k] * L[j][k];
				}
				L[i][j] = (A[i][j] - sum) / L[j][j];
			}
		}
		
	}
	return L;
}



/// <summary>
/// LU implementation
/// </summary>
vector<double> LUDecomp(vector<vector<double>> L, vector<vector<double>> U, vector<double> b)
{
	int n = L.size();

	//set y values to zero
	vector<double> y(n);
	vector<double> x(n);


	//forward substitution
	for (int i = 0; i < n; i++)
	{
		double sum = 0.;
		for (int j = 0; j < i; j++)
		{
			sum += L[i][j] * y[j];
		}
		y[i] = (b[i] - sum) / L[i][i];
	}

	//backward substitution
	for (int i = n - 1; i >= 0; i--)
	{
		double sum = 0.;
		for (int j = i + 1; j < n; j++)
		{
			sum += U[i][j] * x[j];
		}
		x[i] = (y[i] - sum) / U[i][i];
	}

	return x;
}

vector<vector<double>> transpose(vector<vector<double>> A)
{
	int n = A.size();
	vector<vector<double>> transpose(n, vector<double>(n));
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			transpose[j][i] = A[i][j];
		}
	}
	
	return transpose;
}



int main()
{
	vector<vector<double>> Arr ={{8.,	 3.22, 0.8,   0. ,   4.1},
								  {3.22, 7.76, 2.33,  1.91,  -1.03},
								  {0.8,	 2.33, 5.25,  1.00,  3.02},
								  {0.,   1.91, 1.00,  7.50,  1.03},
								  {4.1,  -1.03, 3.02, 1.03, 6.44} };

	vector<double> b = { 9.45,-12.2,7.78,-8.10,10.0 };

	int n = sizeof(Arr) / sizeof(double);


	vector<vector<double>> L = CholeskyDecomposition(Arr,n); // decompose into lowe

	cout << "Cholesky Decomposition of Matrix A is:  \n" << endl;
	for (int i = 0; i < L.size(); i++)
	{
		for (int j = 0; j < L[i].size(); j++)
		{
			cout << L[i][j] << "   " ;
		}
		cout << "\n" << endl;
	}

	cout << "Solution:  \n" << endl;
	vector<vector<double>> U = transpose(L);
	vector<double> x = LUDecomp(L, U, b);

	cout << "TRANSPOSE OF L  is:  \n" << endl;
	for (int i = 0; i < U.size(); i++)
	{
		for (int j = 0; j < U[i].size(); j++)
		{
			cout << U[i][j] << "   ";
		}
		cout << "\n" << endl;
	}


	cout << "Solution for system of equations is :  \n" << endl;
	/// print solution vector
	for (int i = 0; i < x.size(); i++)
	{
		cout << x[i] << " ";
	}

	return 0;
}





















