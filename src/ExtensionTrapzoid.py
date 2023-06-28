# -*- coding: utf-8 -*-
"""
Created on Wed May 24 16:05:10 2023

@author: A.L. Termaat
"""

# For this python file, 
# Define non linear system


#dy1/dt = -a/N * y1*y2           
#dy2/dt = a/N * y1*y2 - b* y2    
#dy3/dt = b*y2                   

import numpy as np 
from numpy.linalg import inv
import matplotlib.pyplot as plt 

# Time Interval
t_start = 0
t_end = 500
n = 10000
dt = (t_end - t_start)/n
dt = 1

iteration = 0

#constants
alpha = 0.0804
beta = 0.0201
gamma = 5.0
delta = 0.8
epsilon = 0.0
phi = 0.0

# Initial Guess Value in millions
N = 3.444 * 10**9
S = 3*10**9
Si = 10**6
Ii = 10**6
Ri = 10**6
R = N - S - Si - Ii - Ri



    
def eq1(S,Si,Ii,Ri,R):
    return -alpha/N * S*(Si+Ii+Ri)

def eq2(S,Si,Ii,Ri,R):
    return alpha/N*S*(Si+Ii+Ri)-gamma/(Si+Ii+Ri) *Si*Ii-beta*Si

def eq3(S,Si,Ii,Ri,R):
    return delta * gamma/(Si+Ii+Ri) *Si*Ii-beta*Ii

def eq4(S,Si,Ii,Ri,R):
    return (1-delta) * gamma/(Si+Ii+Ri) *Si*Ii-beta*Ri

def eq5(S,Si,Ii,Ri,R):
    return beta*(Si+Ii+Ri)

# Equation
#BE : y old + dt * g - y
#FE : y old + dt * g old - y
                                  #TZ : y old + dt/2 (g old + g) - y
#ME : w = y old + dt * g old
    # and y old +dt/2 (g old + g(w)) - y


############

def f1(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt):
	return S_old + (eq1(S,Si,Ii,Ri,R)+eq1(S_old,Si_old,Ii_old,Ri_old,R_old))*dt/2 - S

def f2(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt):
	return Si_old + (eq2(S,Si,Ii,Ri,R)+eq2(S_old,Si_old,Ii_old,Ri_old,R_old))*dt/2 - Si

def f3(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt):
	return Ii_old + (eq3(S,Si,Ii,Ri,R)+eq3(S_old,Si_old,Ii_old,Ri_old,R_old))*dt/2 - Ii

def f4(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt):
    return Ri_old + (eq4(S,Si,Ii,Ri,R)+eq4(S_old,Si_old,Ii_old,Ri_old,R_old))*dt/2 - Ri

def f5(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt):
    return R_old + (eq5(S,Si,Ii,Ri,R)+eq5(S_old,Si_old,Ii_old,Ri_old,R_old))*dt/2 - R

# Jacobian Matrix
def jacobian(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt):
    J = np.zeros((5,5))
    h = 1e-5
    # we will use the limit definition of a derivative, but then approximated by taking h very small.

	# First Row
    J[0,0] = (f1(S_old,Si_old,Ii_old,Ri_old,R_old,S+h,Si  ,Ii  ,Ri  ,R  ,dt)-f1(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[0,1] = (f1(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si+h,Ii  ,Ri  ,R  ,dt)-f1(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[0,2] = (f1(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii+h,Ri  ,R  ,dt)-f1(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[0,3] = (f1(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii  ,Ri+h,R  ,dt)-f1(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[0,4] = (f1(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii  ,Ri  ,R+h,dt)-f1(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    
	# Second Row 
    J[1,0] = (f2(S_old,Si_old,Ii_old,Ri_old,R_old,S+h,Si  ,Ii  ,Ri  ,R  ,dt)-f2(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[1,1] = (f2(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si+h,Ii  ,Ri  ,R  ,dt)-f2(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[1,2] = (f2(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii+h,Ri  ,R  ,dt)-f2(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[1,3] = (f2(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii  ,Ri+h,R  ,dt)-f2(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[1,4] = (f2(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii  ,Ri  ,R+h,dt)-f2(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    
	# Third Row
    J[2,0] = (f3(S_old,Si_old,Ii_old,Ri_old,R_old,S+h,Si  ,Ii  ,Ri  ,R  ,dt)-f3(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[2,1] = (f3(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si+h,Ii  ,Ri  ,R  ,dt)-f3(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[2,2] = (f3(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii+h,Ri  ,R  ,dt)-f3(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[2,3] = (f3(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii  ,Ri+h,R  ,dt)-f3(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[2,4] = (f3(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii  ,Ri  ,R+h,dt)-f3(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h

    # Fourth row
    J[3,0] = (f4(S_old,Si_old,Ii_old,Ri_old,R_old,S+h,Si  ,Ii  ,Ri  ,R  ,dt)-f4(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[3,1] = (f4(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si+h,Ii  ,Ri  ,R  ,dt)-f4(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[3,2] = (f4(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii+h,Ri  ,R  ,dt)-f4(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[3,3] = (f4(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii  ,Ri+h,R  ,dt)-f4(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[3,4] = (f4(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii  ,Ri  ,R+h,dt)-f4(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
	
    # Fifth row
    J[4,0] = (f5(S_old,Si_old,Ii_old,Ri_old,R_old,S+h,Si  ,Ii  ,Ri  ,R  ,dt)-f5(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[4,1] = (f5(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si+h,Ii  ,Ri  ,R  ,dt)-f5(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[4,2] = (f5(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii+h,Ri  ,R  ,dt)-f5(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[4,3] = (f5(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii  ,Ri+h,R  ,dt)-f5(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    J[4,4] = (f5(S_old,Si_old,Ii_old,Ri_old,R_old,S  ,Si  ,Ii  ,Ri  ,R+h,dt)-f5(S_old,Si_old,Ii_old,Ri_old,R_old,S,Si,Ii,Ri,R,dt))/h
    return J



# Newton Rhaphson Method
def Newton_Rhaphson(S,Si,Ii,Ri,R,S_guess,Si_guess,Ii_guess,Ri_guess,R_guess,dt):
	# Vector Equation
    global iteration
    Y_old = np.zeros((5,1))
    Y_old[0] = S_guess
    Y_old[1] = Si_guess
    Y_old[2] = Ii_guess
    Y_old[3] = Ri_guess
    Y_old[4] = R_guess
	#print(Y_old)

    F = np.zeros((5,1))

	# Solver Parameter
    error = 9e9
    tol = 1e-4
    alfa = 1
	#iter = 0

    while(error>=tol):
		# Jacobi set for the equation
        iteration = iteration + 1
        J = jacobian(S,Si,Ii,Ri,R,Y_old[0],Y_old[1],Y_old[2],Y_old[3],Y_old[4],dt)

        F[0] = f1(S,Si,Ii,Ri,R,Y_old[0],Y_old[1],Y_old[2],Y_old[3],Y_old[4],dt)
        F[1] = f2(S,Si,Ii,Ri,R,Y_old[0],Y_old[1],Y_old[2],Y_old[3],Y_old[4],dt)
        F[2] = f3(S,Si,Ii,Ri,R,Y_old[0],Y_old[1],Y_old[2],Y_old[3],Y_old[4],dt)
        F[3] = f4(S,Si,Ii,Ri,R,Y_old[0],Y_old[1],Y_old[2],Y_old[3],Y_old[4],dt)
        F[4] = f5(S,Si,Ii,Ri,R,Y_old[0],Y_old[1],Y_old[2],Y_old[3],Y_old[4],dt)

        Y_new = Y_old - alfa*(np.matmul(inv(J),F))
        error = np.max(np.abs(Y_new-Y_old))
		#print(error,tol)
        Y_old = Y_new
		#print(iter)
        #log_message = "hoi"#'iteration = {0} y1 = {1} y2 = {2} y3 = {3}\'.format(iteration,Y_new[0],Y_new[1],Y_new[2])
		#print(log_message)
		#iter = iter + 1

    return [Y_new[0],Y_new[1],Y_new[2],Y_new[3],Y_new[4]]

# Implicit Euler

def implicit(S,Si,Ii,Ri,R,t_start,t_end,dt):
    t = np.arange(t_start,t_end,dt)
    time = len(t)

    Y_1 = np.zeros(time)
    Y_2 = np.zeros(time)
    Y_3 = np.zeros(time)
    Y_4 = np.zeros(time)
    Y_5 = np.zeros(time)

    Y_1[0] = S
    Y_2[0] = Si
    Y_3[0] = Ii
    Y_4[0] = Ri
    Y_5[0] = R

    S_guess = 1
    Si_guess = 1
    Ii_guess = 1
    Ri_guess = 1
    R_guess = 1

    for i in range (1,time):
        Y_1[i],Y_2[i],Y_3[i], Y_4[i], Y_5[i] = Newton_Rhaphson(Y_1[i-1],Y_2[i-1],Y_3[i-1],Y_4[i-1],Y_5[i-1],S_guess,Si_guess,Ii_guess,Ri_guess,R_guess,dt)
        S_guess = Y_1[i]
        Si_guess = Y_2[i]
        Ii_guess = Y_3[i]
        Ri_guess = Y_4[i]
        R_guess = Y_5[i]
        

    return [t,Y_1,Y_2,Y_3,Y_4,Y_5]

t_lst, S_lst, Si_lst, Ii_lst, Ri_lst, R_lst = implicit(S,Si,Ii,Ri,R,t_start,t_end,dt)

plt.plot(t_lst,S_lst)
plt.plot(t_lst,Si_lst)
plt.plot(t_lst,Ii_lst)
plt.plot(t_lst,Ri_lst)
plt.plot(t_lst,R_lst)
plt.xlabel('Months')
plt.ylabel('Number of people in millions')
plt.legend(['S','Si','Ii','Ri','R'])
title = "Trapezoidal method for δ = " + str(delta) + " and γ = " + str(gamma)
plt.title(title)
plt.grid('on')
plt.show()