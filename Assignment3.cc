/*
	Assignment 3 - Jacobi, Gauss-Seidel, SOR Method
	Name: Isac Rajan
	Roll No: 14ME130
	Course: Applied Computational Methods in Mechanical Sciences
*/

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>

using std::vector;
using std::cout;
using std::cin;
using std::endl;
using std::ofstream;

typedef double mydoub;
typedef vector<double> vdouble1;
typedef vector<vector<double> > vdouble2;

class iterativeMethods
{
	private:
		int N; //No of equations
		mydoub es; //Specified Error Tolerance
		vdouble2 inputMatrix; //Read from a CSV file
		vdouble2 A; //Matrix A
		vdouble1 B; //Matrix B
		
	public:
		int itr;
		int MAX_ITR; //Max iteration
		iterativeMethods(); //constructor function
		void jacobi();
		void gaussSeidel();
		void sor();
};

int main()
{
	int choice;
	cout << "Enter 1 for Jacobi, 2 for Gauss, 3 for SOR: ";
	cin >> choice;
	
	iterativeMethods ans; //class object
	
	switch(choice)
	{
		case 1:
			ans.jacobi();
			break;
		case 2:
			ans.gaussSeidel();
			break;
		case 3:
			ans.sor();
			break;
		default:
			cout << "Wrong Choice::Try Again:(" << endl;
			return -1;
	}
	return 0;
}

iterativeMethods::iterativeMethods(void)
{
	cout << "Enter the number of equations to be solved, present in CSV: ";
	cin >> N;
	
	//Specified Tolerance
	es = 0.05/100;
	
	//Initializing inputMatrix
	inputMatrix.resize(N);
	for(int i=0; i<N; i++)
		inputMatrix[i].resize(N+1);
	
	//Initializing A & B
	A.resize(N);
	B.resize(N);
	for(int i=0; i<N; i++)
		A[i].resize(N);

	//CSV READER: Reading the CSV file containing matrix A & B
	std::ifstream file("inputMatrix.csv");
	for(int rows=0; rows<N; rows++)
	{
		std::string line;
		std::getline(file, line);
		line += ',';
		if(!file.good())
			break;
		
		std::stringstream iss(line);
		for(int cols=0; cols<N+1; cols++)
		{
			std::string val;
			std::getline(iss, val, ',');
			if(!iss.good())
				break;
			
			std::stringstream convertor(val);
			convertor >> inputMatrix[rows][cols];
		}
	}
	
	//Assigning A & B matrix
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			A[i][j] = inputMatrix[i][j];
	for(int i=0; i<N; i++)
			B[i] = inputMatrix[i][N];
		
	//Checking for diagonal dominance
	int c = 0;
	for(int i=0; i<N; i++)
	{
		mydoub sum = 0.0;
		for(int j=0; j<N; j++)
		{
			if (i!=j)
				sum += fabs(A[i][j]);
		}
		if(sum <= fabs(A[i][i]))
			c +=1;
	}
	if (c>=1)
		cout << "The Input Matrix is diagonally dominant\nThus iterative methods are valid\n";
	else
	{
		cout << "Matrix is not diagonally dominant\nProgram Terminated" << endl;
		exit(EXIT_FAILURE);
	}
	
}

void iterativeMethods::jacobi(void)
{
	ofstream result ("JacobiResult.txt");
	result << "Results:: JACOBI ITERATIVE METHOD\n";
	result << "Contents of Matrix A:\n";
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
			result << A[i][j] << "\t";
		result << endl;
	}
	result << "Contents of Matrix B:" << endl;
	for(int i=0; i<N; i++)
		result << B[i] << endl;
	result << "The Input Matrix is diagonally dominant\nThus iterative methods are valid\n";
	
	vdouble1 x_new(N, 0), x_old(N, 50);
	int MAX_ITR = 100;
	
	for(itr=1; itr<=MAX_ITR; ++itr)
	{
		mydoub maxerr = 0.0;
		for(int i=0; i<N; i++)
		{
			mydoub err, sum = 0.0;
			for(int j=0; j<N; j++)
			{
				if(i!=j)
					sum += A[i][j]*x_old[j]; 
			}
			x_new[i] = (B[i]-sum)/A[i][i];
			err = (fabs(x_new[i]-x_old[i])/fabs(x_new[i]))*100;
			if(err > maxerr)
        maxerr = err;				
		}
		result << "Iteration = " << itr << " & Max Error = " << maxerr << endl;
		if (maxerr < es)
		{
			result << "\nSolution converges at Iteration = " << itr << endl;
			for(int m=0; m<N; m++)
				result << "x[" << m+1 << "] = " << x_new[m] << endl;
			result << "Residual Errors:\n";
			for(int p=0; p<N; p++)
			{
				mydoub sum = 0.0;
				for(int q=0; q<N; q++)
					sum += A[p][q]*x_new[q];
				result << "\tError in Equation " << p+1 << " is " << B[p] - sum << endl;
			}
			result.close();
			exit(EXIT_FAILURE);
		}
		for(int n=0; n<N; n++)
			x_old[n] = x_new[n];
	}
	result << "\nSolution didn't converge::Iterations not sufficient" << endl;
	result.close();
}

void iterativeMethods::gaussSeidel(void)
{
	ofstream result ("GaussSeidelResult.txt");
	result << "Results:: GAUSS-SEIDEL ITERATIVE METHOD\n";
	result << "Contents of Matrix A:\n";
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
			result << A[i][j] << "\t";
		result << endl;
	}
	result << "Contents of Matrix B:" << endl;
	for(int i=0; i<N; i++)
		result << B[i] << endl;
	result << "The Input Matrix is diagonally dominant\nThus iterative methods are valid\n";
	
	vdouble1 x(N, 50);
	int MAX_ITR = 100;
	
	for(itr=1; itr<=MAX_ITR; ++itr)
	{
		mydoub maxerr = 0.0;
		for(int i=0; i<N; i++)
		{
			mydoub err, t, sum = 0.0;
			for(int j=0; j<N; j++)
			{
				if(i!=j)
					sum += A[i][j]*x[j]; 
			}
			t = (B[i]-sum)/A[i][i];
			err = (fabs(t-x[i])/fabs(t))*100;
			if(err > maxerr)
        maxerr = err;		
			x[i] = t;
		}
		result << "Iteration = " << itr << " & Max Error = " << maxerr << endl;
		if (maxerr < es)
		{
			result << "\nSolution converges at Iteration = " << itr << endl;
			for(int m=0; m<N; m++)
				result << "x[" << m+1 << "] = " << x[m] << endl;
			result << "Residual Errors:\n";
			for(int p=0; p<N; p++)
			{
				mydoub sum = 0.0;
				for(int q=0; q<N; q++)
					sum += A[p][q]*x[q];
				result << "\tError in Equation " << p+1 << " is " << B[p] - sum << endl;
			}
			result.close();
			exit(EXIT_FAILURE);
		}
	}
	result << "\nSolution didn't converge::Iterations not sufficient" << endl;
	result.close();
}

void iterativeMethods::sor(void)
{
	mydoub omega;
	cout << "Enter the value of of omega (1 < omega < 2) = ";
	cin >> omega;
	
	ofstream result ("SORResult.txt");
	result << "Results:: SOR (Successive Over-Relaxation) ITERATIVE METHOD\n";
	result << "Contents of Matrix A:\n";
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
			result << A[i][j] << "\t";
		result << endl;
	}
	result << "Contents of Matrix B:" << endl;
	for(int i=0; i<N; i++)
		result << B[i] << endl;
	result << "The Input Matrix is diagonally dominant\nThus iterative methods are valid\n";
	
	vdouble1 x_new(N, 0), x_old(N, 50);
	int MAX_ITR = 500;
	
	for(itr=1; itr<=MAX_ITR; ++itr)
	{
		mydoub maxerr = 0.0;
		for(int i=0; i<N; i++)
		{
			mydoub err, t, sum = 0.0;
			for(int j=0; j<N; j++)
			{
				if(i!=j)
					sum += A[i][j]*x_old[j]; 
			}
			t = (B[i]-sum)/A[i][i];
			x_new[i] = omega*t + (1.0-omega)*x_old[i];
			err = (fabs(x_new[i]-x_old[i])/fabs(x_new[i]))*100;
			if(err > maxerr)
        maxerr = err;
			x_old[i] = x_new[i];
		}
		result << "Iteration = " << itr << " & Max Error = " << maxerr << endl;
		if (maxerr < es)
		{
			result << "\nSolution converges at Iteration = " << itr << endl;
			for(int m=0; m<N; m++)
				result << "x[" << m+1 << "] = " << x_new[m] << endl;
			result << "Residual Errors:\n";
			for(int p=0; p<N; p++)
			{
				mydoub sum = 0.0;
				for(int q=0; q<N; q++)
					sum += A[p][q]*x_new[q];
				result << "\tError in Equation " << p+1 << " is " << B[p] - sum << endl;
			}
			result.close();
			exit(EXIT_FAILURE);
		}
		for(int n=0; n<N; n++)
			x_old[n] = x_new[n];
	}
	result << "\nSolution didn't converge::Iterations not sufficient" << endl;
	result.close();
	
}