"""
Assignment 5 - Golden Section Search Method & Newton's Method
Name: Isac Rajan
Roll No: 14ME130
Course: Applied Computational Methods in Mechanical Sciences
"""
import numpy as np

file = open('Results5.txt', 'w')
#Prespecified Tolerance Limit
es = 0.05/100

#Function to be optimized
def f(x):
    return 4*np.sin(x)*(1 + np.cos(x))

#First Derivative
def f1(x):
    return 4*(np.cos(x)+np.cos(2*x))

#Second Derivative
def f2(x):
    return -4*(np.sin(x)+2*np.sin(2*x))
    
#Module for GSSM
def goldenSection():
    file.write("Results: GOLDEN SECTION SEARCH METHOD\n\n")
    #Intial Bounds
    xl = 0
    xu = np.pi/2
    MAX_ITR = 500 #Max iterations, incase of divergence
    for i in range(MAX_ITR):
        file.write("Iteration {}\n".format(i+1))
        d = ((np.sqrt(5)-1)/2)*(xu - xl)
        x1 = xl + d
        x2 = xu - d
        if f(x1) > f(x2):
            xl = x2
        elif f(x1) < f(x2):
            xu = x1
        file.write("\tError = {:.5f}\n".format(xu-xl))
        #Check if it's converged
        if (xu - xl) < es:
            file.write("Solution converged at iteration = {}\n".format(i+1))
            file.write("Final Optimal Solution = {:.4f}\n".format((xu+xl)/2))
            file.write("Max Value of Function = {:.4f}".format(f((xu+xl)/2)))
            break

#Module for Newton's Method
def newtonMethod():
    file.write("\n\nResults: NEWTON'S METHOD\n\n")
    xo = np.pi/4 #Initial Guess
    MAX_ITR = 500
    for i in range(MAX_ITR):
        file.write("Iteration {}\n".format(i+1))
        xn = xo - f1(xo)/f2(xo)
        ea = (xn - xo)/xn
        file.write("\tRelative Approx Error = {:.5f}\n".format(ea))
        if abs(ea) < es and f1(xn) < 1e-6:
            file.write("Solution converged at iteration = {}\n".format(i+1))
            file.write("Final Optimal Solution = {:.4f}\n".format(xn))
            file.write("Max Value of Function = {:.4f}".format(f(xn)))
            break
        xo = xn

#GSSM Method Called        
goldenSection()
#Newton's Method Called
newtonMethod()        

file.close()