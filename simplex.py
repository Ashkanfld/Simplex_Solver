''' This code has been developed by Ashkan Fouladi (fooladiashkang@gmail.com)'''
''' 1) This program Solves linear programs using the simplex algorithm 

2) General Formulas have been adopted as the main method. 

3)To run this program, the user must input the coefficients
 of the standard form of the linear problem in CSV format  

4) In the setting section user could define
the maximum cycle number for finding the optimal point

5) there are three possibilities for the output; first, if the problem is not bounded,
the program prints 'the problem is unbounded'. Second, if there is no feasible solution
for the given constraints, the program prints 'This problem has no feasible solution'.
Third, if the linear problem has an optimal solution, the program prints 
'The result has been reported in the given file path' and reports both the minimum value for
the cost function and the optimal point. 

 '''

# Modules
import numpy as np
from numpy import linalg
import itertools
import csv

#Setting
maxCycle = 100

#Functions
def importData(file_path):
    
    ' This Function imports data from a CSV file'
    with open(file_path,'r',newline='') as f:
        content=csv.reader(f)
        next(content,None)
        c=[]                                # C is the vector of the ceofficients for X 
                                            #in Z=CX
        A=[]                                # A is the constraint matrix 
        b=[]                                # b is the right-hand side vector in AX=b
        general=[]
        for row in content:
            general.append(row[1:])
        m=len(general)-1                    # Number of constraints 
        n=len(general[1])-1                 # Number of Variables
    for i in range(0,n):
        c.append(float(general[0][i]))
    for j in range(1,m+1):
        b.append(float(general[j][-1]))
        A.append(general[j][:-1])
    for r in range(0,m):
        for t in range(0,n):
            A[r][t]=float(A[r][t])
    return A,b,c,m,n

def initialAssign(m,n,A,b):
    ' This Function attempt to find a non-negative pivot.'
    guess=[]
    x=[]
    xn=[]                                   # Nonbasic variables
    xb=[]                                   # Basic Variables
    for i in range(0,n):
        x.append(i)
    seq=itertools.combinations(x, n-m)
    seq_list=[]
    for f in seq:
        seq_list.append(list(f))
    for k in range(0,len(seq_list)):
        terminator=0
        guess=seq_list[k]
        A2=np.delete(A,guess,1)
        Ans=np.linalg.solve(A2,b)
        for w in Ans:
            if w<0:
                terminator+=1
                break
        if terminator==0:
            break
    xn=guess
    for i in x:
        if i not in xn:
            xb.append(i)        
    return xn,xb,terminator

def calculate(m,n,A,b,c,xn,xb):    
    ''' This function seperates the features of variables into 
    nonbasic and basic variables.
    '''
    B=[]                                    # B is the basis matrix (the matrix of the ceofficients
                                            # ceofficients for the basic variables in constraints)
    N=[]                                    # N is the matrix of the ceofficients 
                                            # for the nonbasic variables in constraints                                                   
    for i in range(0,m):
        e=[]
        cb=[]                               # Cb is the vector of the ceofficients 
                                            #for the basic variables in Z function
        for j in range(0,len(xb)):
            e.append(A[i][xb[j]])
            cb.append(c[xb[j]])
        B.append(e)
    for i in range(0,m):
        f=[]
        cn=[]
        for j in range(0,len(xn)):
            f.append(A[i][xn[j]])
            cn.append(c[xn[j]])             #Cn is the vector of the ceofficients
                                            #for the nonbasic variables in Z function
        N.append(f)
    return B,N,cb,cn

def simplex(B,b,cb,cn,N):
    ''' this function operates simplex algorithm '''
    B_inverse=np.linalg.inv(B)              # Calculates the inverse matrix of B
    bHat=np.matmul(B_inverse,b)             # bHat is B(-1)*b
    z=np.matmul(cb,bHat)                    # z is the value for the objective function  
    y=np.matmul(cb,B_inverse)               # y is the vector of simplex mulipliers
    cHat=cn-np.matmul(y,N)                  # cHat is the reduced cost of x
    return B_inverse,bHat,z,cHat
def enteringVariable(cHat):
    ''' This function determines the entering varaible '''
    index=-1
    for j in range(0,len(cHat)):
        if cHat[j]<0:
            index=j                         # Choosing the entering variable
            break
    return index
def leavingVariable(N,m,B_inverse,A,bHat,xt,xn):
    ''' This function determines the leaving varaible '''
    At=np.array([])
    for i in range(0,m):
        At=np.append(At,A[i][xn[xt]])
    At=At.reshape(m,1)
    AHat=np.matmul(B_inverse,At)
    s=1e10
    index=-1
    for i in range(0,len(AHat)):
        if AHat[i]>0 and (bHat[i]/AHat[i])<s:     
            index=i                                 # Choosing the leaving variable
            s=(bHat[i]/AHat[i])
    return index
def update(i,s,xn,xb):
    ''' this function updates the basic and nonbasic varaibles '''
    auxiliary1=xn[i]
    auxiliary2=xb[s]
    xn[i]=auxiliary2
    xb[s]=auxiliary1
    
def reportAnswer(file_path,xb,xn,A,b,z,n,m):
    ''' This function reports the minimum for the cost function 
    and also the optimal basis as a CSV file
    '''
    with open(file_path,'w',newline='') as f:
        content=csv.writer(f)
        a=[]
        for i in range(0,n):
            a.append('X'+str(i+1))
        content.writerow(a+['z'])
        A2=np.delete(A,xn,1)
        aux=np.linalg.solve(A2,b)
        x=[]
        counter=0
        for i in range(0,n):
            if i in xb:
                x.append(aux[counter])
                counter+=1
            if i in xn:
                x.append(0)    
        content.writerow(x+[z])
        
def main():
    
    for p in range(0,maxCycle):
        if p ==0:
            A,b,c,m,n=importData('')
            nonbasic,basic,switch=initialAssign(m,n,A,b)
            if switch !=0:
                print('This problem has no feasible solution')
                break
            B,N,Cb,Cn=calculate(m,n,A,b,c,nonbasic,basic)
            B_inverse,bHat,valueFunction,reducedCost=simplex(B,b,Cb,Cn,N)
            i=enteringVariable(reducedCost)
            if i==-1:
                reportAnswer('',basic,nonbasic,A,b,valueFunction,n,m)
                print('The result has been reported in the given file path')
                break
            s=leavingVariable(N, m, B_inverse, A, bHat, i,nonbasic)
            if s==-1:
                print(' the problem is unbounded')
                break
            update(i,s,nonbasic,basic)
        if p!=0:
            B,N,Cb,Cn=calculate(m,n,A,b,c,nonbasic,basic)
            B_inverse,bHat,valueFunction,reducedCost=simplex(B,b,Cb,Cn,N)
            i=enteringVariable(reducedCost)
            if i==-1:
                reportAnswer('',basic,nonbasic,A,b,valueFunction,n,m)
                print('The result has been reported in the given file path')
                break
            s=leavingVariable(N, m, B_inverse, A, bHat, i,nonbasic)
            if s==-1:
                print(' the problem is unbounded')
                break
            update(i,s,nonbasic,basic)

main()