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

D = np.zeros()