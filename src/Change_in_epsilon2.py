# Program      : Euler's method for a system
# Author       : D.E. Bos, T. Berkhout, en A.L. Termaat
# Created      : 2023

import numpy as np
import matplotlib.pyplot as plt




deltalist = [[0.5,0.2],[0.8,5]]
Iimaxlst1 = []
Iimaxlst2 = []


for k in deltalist:
    ebsilonlst = []
    for j in range(101):
        
        
        N = 3.444 * 10 **(9)
        Dt = 0.05                  # timestep Delta t
        S_init = 3000000000                      # initial population of S
        SI_init = 1000000
        II_init = 1000000                     # initial population of I
        RI_init = 1000000
        R_init = N - S_init - SI_init - II_init - RI_init                     # initial population of R
        t_init = 0                        # initial time
        t_end = 500                      # final time
        
        alpha = 0.02748682735904053
        beta = 0.00458084878685135
        gamma = k[1]
        delta = k[0]
        ebsilon = j/100
        eta = 0.0
    
        n_steps = int(round((t_end-t_init)/Dt)) # total number of timesteps
        
        X = np.zeros(5)                   # create space for current X=[P,G]^T
        dXdt = np.zeros(5)                # create space for current derivative
        t_arr = np.zeros(n_steps + 1)     # create a stor  age array for t
        X_arr = np.zeros((5,n_steps+1))   # create a storage array for X=[P,G]^T
        t_arr[0] = t_init                 # add the initial t to the storage array
        X_arr[0,0] = S_init               # add the initial P to the storage array
        X_arr[1,0] = SI_init               # add the initial G to the storage array
        X_arr[2,0] = II_init
        X_arr[3,0] = RI_init
        X_arr[4,0] = R_init
        
        # Euler's method
        
        for i in range (1, n_steps + 1):
            t = t_arr[i-1]                 # load the time
            S = X_arr[0,i-1]               # load the value of P
            SI = X_arr[1,i-1]               # load the value of G
            II = X_arr[2,i-1]
            RI = X_arr[3,i-1]
            R = X_arr[4,i-1]
            I = SI + II + RI
            X[0] = S                       # fill current state vector X=[P,G]^T
            X[1] = SI
            X[2] = II
            X[3] = RI
            X[4] = R
            dSdt = - (alpha/N) * S * I + eta * R
            dSIdt = (alpha/N) * S * I - (gamma / I) * SI * II - beta * SI                 # calculate the derivative dG/dt
            dIIdt = delta * (gamma / I) * SI * II - beta * II - ebsilon * II
            dRIdt = (1 - delta) * (gamma / I) * SI * II - beta * RI + ebsilon * II
            dRdt = beta * I - eta * R
            dXdt[0] = dSdt                 # fill derivative vector dX/dt
            dXdt[1] = dSIdt         
            dXdt[2] = dIIdt
            dXdt[3] = dRIdt
            dXdt[4] = dRdt
            Xnew = X + Dt*dXdt             # calculate X on next time step
            X_arr[:,i] = Xnew              # store Xnew 
            t_arr[i] = t + Dt              # store new t-value 
            
            
        if k == deltalist[0]:
            Iimaxlst1.append(max(X_arr[2,:]))
        if k == deltalist[1]:
            Iimaxlst2.append(max(X_arr[2,:]))
        
        ebsilonlst.append(ebsilon)




# Plot the results

fig = plt.figure()
plt.plot(ebsilonlst, Iimaxlst1, linewidth = 1, label = 'γ = 0.5, δ = 0.2')
plt.plot(ebsilonlst, Iimaxlst2, linewidth = 1, label = 'γ = 5, δ = 0.8')


  # set title
plt.xlabel('ε', fontsize = 20)   # name of horizontal axis
plt.ylabel('Max. fake news believers', fontsize = 20) # name of vertical axis

plt.xticks(fontsize = 15)               # adjust the fontsize
plt.yticks(fontsize = 15)               # adjust the fontsize
                # set the range of the axes

plt.legend(fontsize=15)                 # show the legend
plt.grid('on')
plt.show()                              # necessary for some platforms

