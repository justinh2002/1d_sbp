import numpy as np
import matplotlib.pyplot as plt
import timeit


# Implement the finite difference approximation for the first derivatives du/dx

def d_x(u, nx, dx, order):
    # finite difference approximation for the first derivatives du/dx
    
    # initialise the derivative vector
    ux = np.zeros(nx) #0*u #np.zeros(u) #0*u
    # second order accurate case
    if order==2:
        
        # consider first the endpoints
        ux[0] = (u[1] -  u[0])/dx                               # at x_0 = 0
        ux[nx-1] = (u[nx-1] -  u[nx-2])/dx                      # at x_N = L  

        # interior points x_j (j = 1, 2, ... N-1)
        for j in range(1, nx-1):
            ux[j] = (u[j+1] -  u[j-1])/(2*dx)
            
            
    # fourth order accurate case        
    if order==4:
        ################################################# 
        # calculate partial derivatives on the boundaries:(0,1,2,3, : nx-4, nx-3, nx-2, nx-1)
        # with one-sided difference operators
        
        ux[0] = -24./17*u[0] + 59./34*u[1]  - 4./17*u[2] - 3./34*u[3]
        ux[1] = -1./2*u[0] + 1./2*u[2] ;
        ux[2] = 4./43*u[0] - 59./86*u[1]  + 59./86*u[3] - 4./43*u[4]
        ux[3] = 3./98*u[0] - 59./98*u[2]  + 32./49*u[4] - 4./49*u[5]


        ux[nx-1] = 24./17*u[nx-1] - 59./34*u[nx-2]  + 4./17*u[nx-3] + 3./34*u[nx-4]
        ux[nx-2] = 1./2*u[nx-1] - 1./2*u[nx-3] ;
        ux[nx-3] = -4./43*u[nx-1] + 59./86*u[nx-2]- 59./86*u[nx-4]+ 4./43*u[nx-5]
        ux[nx-4] = -3./98*u[nx-1] + 59./98*u[nx-3]- 32./49*u[nx-5]+ 4./49*u[nx-6]
    
        # interior points x_j (j = 4, ... nx-5)     
        #------------------------------------------------------------------------------------------------------------------------------
        for j in range(4, nx-4):
            ux[j] = 1./12*u[j-2] - 2./3*u[j-1] + 2./3*u[j+1] - 1./12*u[j+2]

        ux[:] = ux/dx
            
            
    return ux