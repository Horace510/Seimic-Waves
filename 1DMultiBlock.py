import numpy as np
import matplotlib.pyplot as plt
import timeit


'''
Notations:


System of equations:
{
p*u_t +p_x = 0, 1/K*p_t + u_x = F(x,t), xin [0,L], t>=0

p(o,t) = 0, 1/2(Zu(L,t)-p(L,t)) = 0, Z = sqrt(rho K)

u(x,0) = 0, p(x,0) = 0
}


Layer   roh[g/cm]   K[GPa]      Domain[Km]
1       1.5         6.4813      [0,10]
2       1.9         8.3538      [10,20]
3       2.1         10.3719     [20,30]
4       3           9.0         [30,40]

u: prticle velocity [m/s]
p : pressure [MPa]
ro>0: density [g/cm]
K>0: bulk modulus [GPa]

x in {10,20,30}

# DP-SBP-SAT scheme using 4th order DP-SBP and RK4 
'''

rohs = [1.5,1.9,2.1,3]
Ks = [6.4813,8.3538,10.3719,9]
Domains = [(0,10),(10,20),(20,30),(30,40)]
X = [10,20,30]

dimension = 3

B = np.zeros((3,3))
B[0][0] = -1
B[dimension-1][dimension-1] = 1
print(B)




def d_x(u, nx, dx, order):
    # finite difference approximation for the first derivatives du/dx
    
    # initialise the derivative vector
    ux = 0*u
    # second order accurate case
    if order==2:
        
        # consider first the endpoints
        ux[0] = (u[1] -  u[0])/dx                               # at x_0 = 0
        ux[nx-1] = (u[nx-1] -  u[nx-2])/dx                      # at x_N = L  

        # interior points x_j (j = 1, 2, ... N-1)
        for j in range(1, nx-1):
            ux[j] = (u[j+1] -  u[j-1])/(2.0*dx)
            
            
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


# boundary forcing 
def g(t):
    
    import numpy as np
    
    g0 = 0.0

    if t <= 2.0 and t >= 0.0:
        g0 = (np.sin(np.pi/2 * t)) ** 4
        
    
    return g0
