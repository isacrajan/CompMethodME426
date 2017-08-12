/*
	Assignment 2 - Doo-Little & Cholesky Decomposition Method
	Name:	Isac Rajan
	Roll No:	14ME130
	Course: Applied Computational Methods in Mechanical Sciences
*/

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <fstream>

using std::vector;
using std::cout;
using std::cin;
using std::endl;
using std::ofstream;

typedef double mydoub;
typedef vector<double> vdouble1;
typedef vector<vector<double> > vdouble2;

int main()
{
	clock_t start;
	start = clock(); //Variable to measure CPU time
	ofstream result("Assignment2_result.txt");
	result << "Results: DOO-LITTLE DECOMPOSITION METHOD" << endl;
	/* DOO-LITTLE DECOMPOSITION METHOD */
	//Declaration & Initialization of vector A in Ax=B, L, U
	vdouble2 A, L, U;
	A.resize(3);
	L.resize(3);
	U.resize(3);
	for(int i=0; i<3; i++)
	{
		A[i].resize(3);
		L[i].resize(3);
		U[i].resize(3);
	}
	A[0][0] = 1;	A[0][1] = 4; A[0][2] = 1;
	A[1][0] = 1;	A[1][1] = 6; A[1][2] = -1;
	A[2][0] = 2;	A[2][1] = -1; A[2][2] = 2;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
		{
			if(i == j)
				L[i][j] = 1.0; //Diagonal Elements
			if(i < j)
				L[i][j] = 0.0; //Upper Triangular Matrix
			if(i > j)
				U[i][j] = 0.0; //Lower Triangular Matrix
		}
	//Declaration & Initialization of vector x, z, B
	vdouble1 x(3, 0), z(3, 0), B(3);
	B[0] = 7;
	B[1] = 13;
	B[2] = 5;
	result << "Contents of Matrix A:\n";
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
			result << A[i][j] << "\t";
		result << endl;
	}
	result << "Contents of Matrix B:" << endl;
	for(int i=0; i<3; i++)
		result << B[i] << endl;

	//1. Getting L and U matrix
	for(int k=0; k<3; k++)
	{
		for(int m=k; m<3; m++)
		{
			mydoub sum;
			sum = 0.0;
			for(int j=0; j<=(k-1); j++)
				sum += L[k][j]*U[j][m];
			U[k][m] = A[k][m] - sum;
		}
		for(int i=k+1; i<3; i++)
		{
			mydoub sum;
			sum = 0.0;
			for(int j=0; j<=(k-1); j++)
				sum += L[i][j]*U[j][k];
			L[i][k] = (A[i][k] - sum)/U[k][k];
			L[i][i] = 1;
		}
	}
	result << "Contents of Matrix U:\n";
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
			result << U[i][j] << "\t";
		result << endl;
	}
	result << "Contents of Matrix L:\n";
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
			result << L[i][j] << "\t";
		result << endl;
	}
	//3. Solve LZ=B to get Z
	for(int i=0; i<3; i++)
	{
		mydoub sum;
		sum = 0.0;
		for(int j=0; j<=(i-1); j++)
			sum += L[i][j]*z[j];
		z[i] = B[i] - sum;
	}

	result << "Contents of Matrix Z:" << endl;
	for(int i=0; i<3; i++)
		result << z[i] << endl;

	//4. Solve UX=B to get X
	x[2] = z[2]/U[2][2];
	for(int i=1; i>=0; i--)
	{
		mydoub sum;
		sum = 0.0;
		for(int j=i+1; j<3; j++)
			sum += U[i][j]*x[j];
		x[i] = (z[i]-sum)/U[i][i];
	}
	result << "After solving, the values of x1, x2, and x3 are: " << endl;
	for(int i=0; i<3; i++)
		result << "x" << i+1 << " = " << x[i] << endl;
	result << "CPU time: " << (double)((-start + clock())/CLOCKS_PER_SEC) << " seconds" << endl;
	result << "\nVerification: Substituting results into original equation:\nError is:" << endl;
	for(int i=0; i<3; i++)
		result << "\tError in Equation " << i+1 << " is " << B[i] - (A[i][0]*x[0]+A[i][1]*x[1]+A[i][2]*x[2]) << endl;

	/* --------------------------------------------------- */
	/* CHOLESKY DECOMPOSITION METHOD */
	result << "\n\n---------------------------------------" << endl;
	result << "Results: CHOLESKY DECOMPOSITION METHOD" << endl;
	start = clock();
	//Initialization of Vector A
	A[0][0] = 10;	A[0][1] = -1; A[0][2] = 2;
	A[1][0] = -1;	A[1][1] = 12; A[1][2] = -1;
	A[2][0] = 2;	A[2][1] = -1; A[2][2] = 20;
	//Initialization of vector B
	B[0] = 11;
	B[1] = 10;
	B[2] = 21;
	result << "Contents of Matrix A:\n";
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
			result << A[i][j] << "\t";
		result << endl;
	}
	result << "Contents of Matrix B:" << endl;
	for(int i=0; i<3; i++)
		result << B[i] << endl;

	//Checking for symmetry in Matrix A
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			if(i!=j)
				if(A[i][j] != A[j][i])
				{
					result << "\nMatrix A is not symmetric" << endl;
					return -1;
				}
	result << "\nCheck Result: Matrix A is symmetric" << endl;
	//Checking positive definitivity of Matrix A
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			if(i==j)
				if(A[i][i] <= 0)
				{
					result << "\nMatrix A in not Positive Definite" << endl;
					return -1;
				}
	result << "Check Result: Matrix A is Positive Definite\n" << endl;

	//1. Get U Matrix and L=U^T
	for(int i=0; i<3; i++)
	{
		mydoub sum;
		sum = 0.0;
		for(int k=0; k<=i-1; k++)
			sum += U[k][i]*U[k][i];
		U[i][i] = sqrt(A[i][i] - sum);

		for(int j=i+1; j<3; j++)
		{
			sum = 0.0;
			for(int k=0; k<=i-1; k++)
				sum += U[k][i]*U[k][j];
			U[i][j] = (A[i][j] - sum)/U[i][i];
		}
	}
	result << "Contents of Matrix U:\n";
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
			result << U[i][j] << "\t";
		result << endl;
	}

	//2. Solving U^TZ = B
	for(int i=0; i<3; i++)
	{
		mydoub sum;
		sum = 0.0;
		for(int j=0; j<=i-1; j++)
			sum += U[j][i]*z[j];
		z[i]=(B[i]-sum)/U[i][i];
	}
	result << "Contents of Matrix Z:" << endl;
	for(int i=0; i<3; i++)
		result << z[i] << endl;

	//3. Solving UX=Z
	for(int i=2; i>=0; i--)
	{
		mydoub sum;
		sum = 0.0;
		for(int j=i+1; j<3; j++)
			sum += U[i][j]*x[j];
		x[i]=(z[i]-sum)/U[i][i];
	}
	result << "After solving, the values of x1, x2, and x3 are: " << endl;
	for(int i=0; i<3; i++)
		result << "x" << i+1 << " = " << x[i] << endl;
	result << "CPU time: " << (double)((-start + clock())/CLOCKS_PER_SEC) << " seconds" << endl;
	result << "\nVerification: Substituting results into original equation:\nError is:" << endl;
	for(int i=0; i<3; i++)
		result << "\tError in Equation " << i+1 << " is " << B[i] - (A[i][0]*x[0]+A[i][1]*x[1]+A[i][2]*x[2]) << endl;
	result.close();
	return 0;
}
