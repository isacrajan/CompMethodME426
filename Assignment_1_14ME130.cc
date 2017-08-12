/*
	Assignment 1 - Gauss Elimination
	Name:	Isac Rajan
	Roll No:	14ME130
	Course: Applied Computational Methods in Mechanical Sciences
*/

#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>

using std::vector;
using std::cout;
using std::cin;
using std::endl;
using std::ofstream;

typedef double mydoub;
typedef int myint;
typedef vector<double> vdouble1;
typedef vector<vector<double> > vdouble2;

int main()
{
	mydoub N;
	cout << "Enter the number of equation to be solved: ";
	cin >> N;
	clock_t start; //Variable to measure CPU time
	
	//Input Time and velocities
	vdouble1 t(N), v(N);
	cout << "Enter " << N << " times and velocities:\n" << endl;
	for(int i=0; i<N; i++)
		cin >> t[i] >> v[i];
	
	//Declaration & Initialization of Vector A
	vdouble2 A;
	A.resize(N);
	for(int i=0; i<N; i++)
		A[i].resize(N);
	for(int i=0; i<N; i++)
	{
		A[i][0] = t[i]*t[i];
		A[i][1] = t[i];
		A[i][2] = 1;
	}	
	
	//Declaration & Initialization of Vector x & B
	vdouble1 x(N, 0), B(N);
	for(int i=0; i<N; i++)
	{
		B[i] = v[i];
	}
	
	ofstream result("Assignment1_result.txt");
	result << "The matrix A is as shown below:" << endl;
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
			result << A[i][j] << "\t";
		result << endl;
	}
	result << "The matrix B is as shown below:" << endl;
	for(int i=0; i<N; i++)
		result << B[i] << endl;
	
	start = clock();
	mydoub ratio, sum;
	//Forward Elimination
	for(int i=0; i<N-1; i++)
		for(int k=i+1; k<N; k++)
		{
			ratio = A[k][i]/A[i][i];
			for(int j=0; j<N; j++)
				A[k][j]-=ratio*A[i][j];
			B[k]-=ratio*B[i];
		}
	
	//Backward Substitution
	x[N-1] = B[N-1]/A[N-1][N-1];
	for(int i=N-2; i>=0; i--)
	{
		sum = 0.0;
		for(int j=i+1; j<N; j++)
			sum+=A[i][j]*x[j];
		x[i] = (B[i]-sum)/A[i][i];
	}
	
	result << "After solving, the values of a1, a2, and a3 are: " << endl;
	for(int i=0; i<N; i++)
		result << "a" << i+1 << " = " << x[i] << endl;
	
	result << "\nThe velocity at t = 6 sec is " << x[0]*6*6+x[1]*6+x[2] << endl;
	result << "CPU time: " << (double)((-start + clock())/CLOCKS_PER_SEC) << " seconds" << endl;
	
	result << "\nVerification: Substituting results into original equation:\nError is:" << endl;
	for(int i=0; i<N; i++)
		result << "\tError in Equation " << i+1 << " is " << v[i] - (x[0]*t[i]*t[i]+x[1]*t[i]+x[2]) << endl;
	
	result.close();
	return 0;
}

