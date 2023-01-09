import numpy as np
import matplotlib.pyplot as plt
import timeit
from derivative_sbp import * 
from boundary import * 
from sat import *

def advect_test(order, cfl,tend,nx):
    L = 10.0                               # length of the domain (km)
    c = 1.0                                # wavespeed (m/s)
    t = 0.0                                # initial time
    tend =  tend #1 #6.5                              # final time

    # set discretisation parameterscfl
    nx = nx #101                                # number of gridpoints

    # discretise the domain into nx grid points
    x = np.linspace(0, L, nx)              # discrete domain
    dx = x[1]-x[0]                                                

    order = float(order)   
    #Initialise the arrays 
    u =  np.zeros(nx) 

    U =  np.zeros(nx)
                           # order of accuracy: 2 or 4.    

# Time stepping parameters
    cfl = float(cfl)                      # CFL number
    dt = cfl/c*dx                     # Time step
    nt = int(tend/dt)          # number of time steps
    n = 0                             # counter
    t=0 

    print("cfl is: ", cfl)
    print("time step is: ", dt)
    print("number of time steps is:", nt)

    method = input("Euler, RK2 or RK4: ")
    input('press enter to begin simulation')

    # initialise timer
    start = timeit.default_timer()

    # loop in time
    for t in np.arange(0, tend, dt):
        n = n+1
        
        #forward Euler
        if method == "Euler":
            u = u + dt*rate(u, nx, dx, dt, order, c, t, x)
            
        # RK2
        if method == "RK2":
            # compute predictor
            k1 = u + dt/2*rate(u, nx, dx, dt, order, c, t, x)
            
            # corrector
            u = u + dt*rate(u+dt/2*k1, nx, dx, dt, order, c, t+0.5*dt, x)
            
        
        # RK4
        if method == "RK4":
            # compute predictor (RK-stages)
            k1 = rate(u, nx, dx, dt, order, c, t, x)
            k2 = rate(u+dt/2*k1, nx, dx, dt, order, c, t+0.5*dt, x)
            k3 = rate(u+dt/2*k2, nx, dx, dt, order, c, t+0.5*dt, x)
            k4 = rate(u+dt*k3, nx, dx, dt, order, c, t+dt, x)
            
            # corrector
            u = u + dt/6*(k1 + 2*k2 + 2*k3 + k4)

            import csv
        
            with open('csv_datafile{}.csv'.format(round(t,1)), 'w') as f:
            # create the csv writer
                writer = csv.writer(f)

                # write a row to the csv file
                for w in range(len(u)):
                    writer.writerow([x[w],u[w]])

            print('outputting file for time: ', round(t,2), ', iteration:', n   )
            
        
        # Exact solution
        #U=np.exp(-(((x-c*(t+dt)-x0)/delta)**2))ls

        for j in range(nx):
            U[j]=g(t+dt-x[j]/c) 
        



    # Simulation end time
    stop = timeit.default_timer()


    print('total simulation time = ', stop - start, 's')            # print the time required for simulation
    print('spatial order  of accuracy = ', order)                   # print the polynomial degree used
    print('number of grid points = ', nx)                           # print the degree of freedom                    # print the spatial step
    print('total number of time steps',nt )
    print('maximum relative error = ', np.max(np.abs(U-u))/np.max(np.abs(U)))             # print the max. relative error
    print('log2 max relative error = ', np.log2(np.max(np.abs(U-u))/np.max(np.abs(U)))) 
    
    return   


import argparse
# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='SBP.code.advection.')
    parser.add_argument("-simulation", "--simulation",
                        choices=['advect_test'],
                        default='advect_test', help="The simulation setup/test (default: %(default)s).")
    parser.add_argument("-order", "--order", choices=['2', '4'], default = 4,
                        help="Finite difference scheme for advection (default: %(default)s).")
    parser.add_argument("-cfl", "--cfl", type=float,default = 1.0, 
                        help="CFL timestep safety factor (default: %(default)s).")
    parser.add_argument("-tend", "--tend", type=float,default = 1.0, 
                        help="Simulation end time (default: %(default)s).")
    parser.add_argument("-nx", "--nx", type=int,default = 101, 
                        help="Number of grid points (default: %(default)s).")                                             

    args = parser.parse_args()

    # new hydro class object

    # advection test
    if args.simulation == "advect_test":
        advect_test(args.order,args.cfl,args.tend,args.nx)



