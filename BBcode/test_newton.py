from scipy.integrate import quad 
from scipy.optimize import newton
import numpy as np

th = np.sqrt(np.pi)/2

def gaussian(x, sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-x**2/(2*sigma)**2)

def p_corr(sigma):
    return quad(gaussian, -th, th, args=(sigma,))[0]    

def root_func(sigma, p):
    return p_corr(sigma) + p - 1

def find_sigma(p):
    sol = newton(root_func, x0=1, args=(p,))
    return sol


sigma = find_sigma(0.1)
print(sigma)
print(p_corr(sigma))