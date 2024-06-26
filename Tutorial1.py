# Finding the roots using Newton Raphson
from scipy.misc import derivative
import numpy as np

# Initial conditions
p1 = 1
rho1 = 1
u1 = 0

p6 = 0.1
rho6 = 0.125
u6 = 0

gamma = 1.4
cs1 = np.sqrt(gamma*p1/rho1)
cs6 = np.sqrt(gamma*p6/rho6)

A6 = (2/(rho6*(gamma+1)))
B6 = p6*((gamma-1)/(gamma+1))

def eqn(x):
    pow = (gamma-1)/(2*gamma)
    y = cs6*(x - p6) * np.sqrt(A6/(x + B6)) + ((2*cs1)/(gamma-1))*((x/p1)**(pow) - 1)
    return y

def d_deqn(x):
    dy_dx = cs6 * np.sqrt(A6/(x + B6)) - cs6 * (1/2)*((x - p6)/(x+B6)) * np.sqrt(A6/(x + B6)) + (cs1/gamma)*(x/p1)**(-(gamma+1)/(2*gamma))
    return dy_dx

def new_raph(x0):
    x = x0
    while abs(eqn(x))>1e-4:
        fx = eqn(x)
        dfx = d_deqn(x)
        x = x - (fx/dfx)
        print('value of x: ', x)
    return x

p3 = new_raph(1)
print('Value of pressure in region 3: ', p3)

# Calculating the value of lambda
# lamb = 0.78*((0.4/2.8)+((2.4*0.36)/(2.8*0.1)))**0.5
# print('Lambda: ', lamb)

# # Calculating the value of u3 
# u3 = 0.78 * ((p3/0.1)-1) * (1.43/(2.4*(p3/0.1)+0.4))**0.5
# print('Value of u3: ', u3)

# # Calculating the shock position
# x0 = 0.5
# t = 0.414
# xsh = x0 + lamb*t

