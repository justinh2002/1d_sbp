import numpy as np
from derivative_sbp import * 
from boundary import * 
# Compute the RHS:  -cdu/dx + SAT
def rate(u, nx, dx, dt, order, c, t, x):
    tau0 = 3

    # set penalty parameter
    if order == 2:
        tau = tau0/(0.5 * dx)
    if order == 4:
        tau =  tau0/((17.0 / 48.0) * dx)
        #print(tau)
        
    #  compute numerical derivative of u 
    r = -c*d_x(u, nx, dx, order)
    
    # SAT term is - H^-1 tau0 e1 ( boundary terms) where H is the positive definite, symmetric discrete norm
    # penalize the boundary term
    r[0] -=  tau*(u[0]-g(t))
    
    return r