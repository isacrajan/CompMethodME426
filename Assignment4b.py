"""
Assignment 4(b) -Fixed-Point Iteration Method 
Name: Isac Rajan
Roll No: 14ME130
Course: Applied Computational Methods in Mechanical Sciences
"""
import numpy as np

ho = 5 #Initial Guess of H
hn = 0.0 #New Value of H- Initalization

es = 0.05/100 #Prespecified Tolerance
ea = 0.0 #Relative Absolute error - Initialization

MAX_ITR = 500
f = open('Results4b.txt', 'w')
f.write("Results: FIXED-POINT ITERATION METHOD\n")

for i in range(1,MAX_ITR+1):
	f.write("Iteration {}\n".format(i))
	hn = ((0.0135*ho**2 + 0.27*ho + 1.35)/9.05097) ** 0.2 #Formulation
	#hn = ((9.05097*ho**5 - 0.27*ho - 1.35)/0.0135) ** 0.5
	#hn = (9.05097*ho**5 - 0.0135*ho**2 - 1.35)/0.27
	ea = 100 * (hn-ho)/hn
	ho = hn
	f.write("\tRelative Aprrox Error(%) = {0:.6f}\n".format(abs(ea)))
	if abs(ea) < es:
		f.write("Solution is {:.5f}: Reached at Iteration = {}".format(hn,i))
		break

f.write("\n\nSubtituting the solution in equation of Flow Rate:\n")
f.write("We get: Q = {}".format( ((0.0002**0.5) * (20 * hn) ** (5./3))/(0.03 * (20 + 2*hn) ** (2./3)) ))		
f.close()	