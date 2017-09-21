"""
Assignment 6 - Newton's Divided Difference & Lagrange Interpolation Method
Name: Isac Rajan
Roll No: 14ME130
Course: Applied Computational Methods in Mechanical Sciences
"""

import numpy as np
import sys
import timeit

#Input Data is read from a file
T, alpha = np.loadtxt("inputData.txt", unpack=True)

file = open('Results6.txt', 'w')

def newtonDividedDiff():
	file.write("Results: NEWTON'S DIVIDED DIFFERENCE METHOD\n\n")
	file.write("Divided Difference Table:\n")
	t = float(input("Enter the Temperature at which alpha is to be found:"))
	start_time = timeit.default_timer()
	n = T.size #Number of data points
	dd = np.zeros((n,n))
	#Evaluating the Divided Difference Table
	for i in range(n):
		dd[i,0] = alpha[i]
	for j in range(1,n):
		for i in range(0,n-j):
			dd[i,j] = (dd[i+1,j-1]-dd[i,j-1])/(T[i+j]-T[i])
	for i in range(n):
		file.write("{}\t" .format(T[i]))
		for j in range(n-i):
			file.write("{:.5e}\t".format(dd[i,j]))
		file.write("\n")
	#Calculation of the interpolated Value
	sum = dd[0,0]
	for i in range(1,n):
		xterm = 1.0
		for j in range(0,i):
			xterm = xterm * (t - T[j])
		sum += dd[0,i]*xterm	
	file.write("\nValue of alpha at T = {} is {:5e}".format(t,sum))
	file.write("\nTime Taken = {:3f} seconds".format(timeit.default_timer() - start_time))
	file.close()

def lagrange():
	file.write("Results: LAGRANGE INTERPOLATION METHOD\n\n")
	t = float(input("Enter the Temperature at which alpha is to be found:"))
	start_time = timeit.default_timer()
	n = T.size #Number of data points
	sum = 0.0
	for i in range(n):
		mul = 1.0
		for j in range(n):
			if j!=i:
				mul *= (t-T[j])/(T[i]-T[j])
		sum += mul*alpha[i]
	file.write("\nValue of alpha at T = {} is {:5e}".format(t,sum))
	file.write("\nTime Taken = {:3f} seconds".format(timeit.default_timer() - start_time))
	file.close()

#Main Program Starts Here
choice = int(input("Enter 1 for Newton's Divided Difference & 2 for Lagrange Interpolation:"))
if choice == 1:
	newtonDividedDiff()
elif choice == 2:
	lagrange()
else:
	print("Invalid Choice::Program Terminated")
	sys.exit()
