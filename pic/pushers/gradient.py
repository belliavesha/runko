import numpy as np

Nx = 100
Ny = 100

grid_x = np.linspace(0,1,Nx)
grid_y = np.linspace(0,1,Nx)
vals = np.zeros((Nx,Ny))

def F(x,y):
    return x**2*np.exp(y), np.cos(x)/(1+y)

def Fdx(x,y):
    return 2*x*np.exp(y), -np.sin(x)/(1+y) 

def Fdy(x,y):
    return x**2*np.exp(y), -np.cos(x)/(1+y)/(1+y)

def comp_grads(vals):
    return dxx,dyy,dxy,dyx
