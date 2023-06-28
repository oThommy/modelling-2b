# Program      : Euler's method for a system
# Author       : D.E. Bos, T. Berkhout, and A.L. Termaat
# Created      : 2023

import numpy as np
import matplotlib.pyplot as plt




Dt = 0.01                   # timestep Delta t
S_init = 5000000000                       # initial population of S
I_init = 100000000                      # initial population of I
R_init = 0                          # initial population of R
t_init = 0                        # initial time
t_end = 150                      # final time

alpha = 0.00000000002
beta = 0.02
gamma = 0 / 212

n_steps = int(round((t_end-t_init)/Dt)) # total number of timesteps

X = np.zeros(3)                   # create space for current X=[P,G]^T
dXdt = np.zeros(3)                # create space for current derivative
t_arr = np.zeros(n_steps + 1)     # create a stor  age array for t
X_arr = np.zeros((3,n_steps+1))   # create a storage array for X=[P,G]^T
t_arr[0] = t_init                 # add the initial t to the storage array
X_arr[0,0] = S_init               # add the initial P to the storage array
X_arr[1,0] = I_init               # add the initial G to the storage array
X_arr[2,0] = R_init

# Euler's method

for i in range (1, n_steps + 1):
    t = t_arr[i-1]                 # load the time
    S = X_arr[0,i-1]               # load the value of P
    I = X_arr[1,i-1]               # load the value of G
    R = X_arr[2,i-1]
    X[0] = S                       # fill current state vector X=[P,G]^T
    X[1] = I
    X[2] = R
    dSdt = - alpha * S * I + gamma * R
    dIdt = alpha * S * I - beta * I                 # calculate the derivative dG/dt
    dRdt = beta * I - gamma * R
    dXdt[0] = dSdt                 # fill derivative vector dX/dt
    dXdt[1] = dIdt         
    dXdt[2] = dRdt
    Xnew = X + Dt*dXdt             # calculate X on next time step
    X_arr[:,i] = Xnew              # store Xnew 
    t_arr[i] = t + Dt              # store new t-value 

# Plot the results

fig = plt.figure()
plt.plot(t_arr, X_arr[0,:], linewidth = 4, label="S") 
plt.plot(t_arr, X_arr[1,:], linewidth = 4, label="I")  
plt.plot(t_arr, X_arr[2,:], linewidth = 4, label="R")

  # set title
plt.xlabel('t', fontsize = 20)   # name of horizontal axis
plt.ylabel('S, I and R', fontsize = 20) # name of vertical axis

plt.xticks(fontsize = 15)               # adjust the fontsize
plt.yticks(fontsize = 15)               # adjust the fontsize
                # set the range of the axes

plt.legend(fontsize=15)                 # show the legend
plt.show()                              # necessary for some platforms

