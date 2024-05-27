import numpy as np
from scipy.optimize import minimize

def rosen(x):
    """The Rosenbrock function"""
    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
for method in ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP']:
    for i in range(10):
        res = minimize(rosen, x0, method=method,
                options={'xatol': 1e-8, 'disp': True})
        if res.success:
            print(method, res.x)
            break

# returns x- answer, success 
# - bool if optimization was successful, 
# status - int, 
# message - str, 
# fun, jas, hess - value of the function at the minimum, 
# nfev - number of function evaluations, 
# nit - number of iterations. 
# maxcv - max constraint violation