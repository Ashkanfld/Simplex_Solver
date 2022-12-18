# Simplex Solver

This program Solves linear programs using the simplex algorithm (general formulas)

# How to use 

1. Download a Python file (Simplex.py)

2. Create a CSV file for the linear problem in the standard form (standard form is defined below)

3. Determine a file path for the CSV file in the " importData " function.

4. Determine a file path for reporting the minimum of the cost function and optimal point in the " reportAnswer " function. 

Note: Currently, only linear programs in standard form are supported.


# Standard form 

Z=CX

AX=b 

X>0

Where: 

Z = Cost Function 

C = vector of coefficients 

A = the constraint matrix 

b = the right-hand side vector 

Note: all components of the vector 'b' must be Non-negative.

# Example (Igor Griva, Stephen G. Nash, Ariela Sofer, Linear and Nonlinear Optimization, p. 126.)

1) basic form 

Minimize      Z = -X1 - 2X2

subject to    -2X1 + X2 <= 2
	       -X1 + 2X2<= 7
	             X1 <= 3
             X1,X2>=0

2) standard form 

Minimize   	Z = -X1 - 2X2

subject to     -2X1 + X2 + X3 = 2
	       -X1 + 2X2 + X4 = 7
	              X1 + X5 = 3
             X1,X2,X3,X4,X5>=0              
So; 

A = [[-2,1,1,0,0],[-1,2,0,1,0],[1,0,0,0,1]]

C = [-1,-2,0,0,0]

b=[2,7,3]

3) the CSV file for importing the problem is Example.csv

4) the results could be seen in Results.csv


# Acknowledgments

The entire program is written by Ashkan Fouladi (fooladiashkang@gmail.com)

