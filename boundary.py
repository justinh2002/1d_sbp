# consider the boundary forcing
def g(t):

    import numpy as np
    
    g0 = 0.0

    if t <= 2.0 and t >= 0.0:
        g0 = (np.sin(np.pi/2 * t)) ** 4
        
    
    return g0