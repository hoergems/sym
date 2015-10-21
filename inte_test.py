from libintegrate import *
import time
import numpy as np

class InteTest:
    def __init__(self):
        integrate = Integrate(2.0)
        thetas_star = v_double()
        dot_thetas_star = v_double()
        rho_star = v_double()
        
        current_state = v_double()
        
        
        thetas_star[:] = [np.pi / 2.0, 0.0]
        dot_thetas_star[:] = [0.0, 0.0]
        rho_star[:] = [0, -0.1]
        current_state[:] = [np.pi / 2.0, 0.0, 0.0, 0.0]
        
        t0 = 0.0
        te = 0.03
        delt = 0.03
        
        int_times = v_double()
        int_times[:] = [t0, te, delt]
        
        integrate.setup(thetas_star, dot_thetas_star, rho_star)        
        t0 = time.time()
        integrate.doIntegration(current_state, int_times)
        print "integration took " + str(time.time() - t0) + " seconds"
        
if __name__ == "__main__":
    InteTest()