"""
Assignment 4(a) - Newton Raphson's Method 
Name: Isac Rajan
Roll No: 14ME130
Course: Applied Computational Methods in Mechanical Sciences
"""

import numpy as np

ho = -1.0 #Old value of x - INITIAL GUESS
hn = 0.0 #New value of x - Initialization

es = 1e-5 #Prespecified Tolerance Limit
ea = 0.0 #Relative Absolute Error - Initialization

MAX_ITR = 500
f = open('Results4a.txt', 'w')
f.write("Results: NEWTON-RAPHSON METHOD\n")

for i in range(1,MAX_ITR):
	f.write("Iteration {}\n".format(i))
	hn = (2*ho**3 - 9*ho**2 - 28.64789)/(3*ho**2 - 18*ho) #Formulation
	ea = (hn-ho)/hn
	ho = hn
	f.write("\tRelative Aprrox Error = {0:.6f}\n".format(abs(ea)))
	if abs(ea) < es:
		f.write("Solution is {:.5f}: Reached at Iteration = {}".format(hn,i))
		break

f.write("\n\nSubtituting the solution in equation of Volume:\n")
f.write("We get: V = {}".format((np.pi * hn**2 * (3*3 - hn) )/3))		
f.close()	