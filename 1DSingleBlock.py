# Parameters initialization and plotting the simulation
# Import necessary routines
import numpy as np
import matplotlib.pyplot as plt
import time_integrator
import rate
import utils
import timeit

#plt.switch_backend("TkAgg")          # plots in external window
# plt.switch_backend("nbagg")           # plots within this notebook
'''
rho*u_t + p_x = 0, 1/K*p_t + u_x = 0, x in [0,L], t>=0
p(0,t) = 0, 1/2(Z*u(L,t)-p(L,t))=0, Z = sqrt(rho*K)
u(x,0)=f(x), p(x,0)=g(x)



'''






# Initializations
L = 10.0         # length of the domain (km)
t = 0.0          # initial time
tend = 1.45      # final time
nx = 501        # grid points in x                                                                                                                       
dx = L/(nx-1)    # grid increment in x
cs = 3.464       # velocity (km/s) (can be an array)                                                                                                             
iplot = 5       # snapshot frequency
rho = 1    # density [g/cm^3]
mu = rho*cs**2   # shear modulus [GPa]
Zs = rho*cs      # shear impedance 

order = 4        # order of accuracy


#Initialize the domain
y = np.zeros((nx, 1))

# Initial particle velocity perturbation and discretize the domain
for j in range(0, nx):
    y[j, :] = j*dx                                             # discrete domain



# Time stepping parameters
cfl = 1.0                         # CFL number
dt = (cfl/cs)*dx                  # Time step
nt = int(round(tend/dt))          # number of time steps
n = 0                             # counter

# Boundary condition reflection coefficients 
r0 = 1                            # r=0:absorbing, r=1:free-surface, r=-1: clamped 
r1 = 1                            # r=0:absorbing, r=1:free-surface, r=-1: clamped

# penalty parameters
tau_11 = 1 
tau_12 = 1
tau_21 = 1 
tau_22 = 1

# Initialize: particle velocity (v); and shear stress (s)

##############################
v = np.zeros((nx, 1))       ##
s = np.zeros((nx, 1))       ##
                            ##
U = np.zeros((nx, 1))       ##
V = np.zeros((nx, 1))       ##
##############################
U_t = np.zeros((nx, 1))
V_t = np.zeros((nx, 1))
U_x = np.zeros((nx, 1))
V_x = np.zeros((nx, 1))

                                

# Difference between analyticla and numerical solutions
EV = [0]                                 # initialize errors in V (velocity)
EU = [0]                                 # initialize errors in U (stress)
T = [0]                                  # later append every time steps to this



# Computation and plotting

# Initialize animated plot for velocity and stress
# fig1 = plt.figure(figsize=(10,10))
# ax1 = fig1.add_subplot(4,1,1)
# line1 = ax1.plot(y, v, 'r', y, U, 'k--')
# plt.title('numerical vs exact')
# plt.xlabel('x [km]')
# plt.ylabel('velocity [m/s]')

# ax2 = fig1.add_subplot(4,1,2)
# line2 = ax2.plot(y, s, 'r', y, V, 'k--')
# plt.title('numerical vs exact')
# plt.xlabel('x[km]')
# plt.ylabel('stress [MPa]')

# # Initialize error plot (for velocity and stress)
# ax3 = fig1.add_subplot(4,1,3)
# line3 = ax3.plot(T, EV, 'r')
# plt.title('relative error in particle velocity')
# plt.xlabel('time [s]')
# ax3.set_ylim([10**-5, 1])
# plt.ylabel('error')

# ax4 = fig1.add_subplot(4,1,4)
# line4 = ax4.plot(T, EU, 'r') 
# plt.ylabel('error')
# plt.xlabel('time[t]')
# ax4.set_ylim([10**-5, 1])
# plt.title('relative error in stress')

# plt.tight_layout()
# plt.ion()
# plt.show()


t=0   # initial time

forcing = 1.0  # forcing function, forcing = 1,  and no forcing function, forcing = 0

# type of initial data: Gaussian or Sinusoidal
type_0 = 'Gaussian'
#type_0 = 'Sinusoidal'


if type_0 in ('Sinusoidal'):
    forcing = 1.0  # we must use forcing for Sinusoidal initial condition



# L2-norm normalizer
# Generate  conditions for normalization
rate.mms(v, s, U_t, V_t, U_x, V_x, y, 0.65, type_0)
A =  (np.linalg.norm(v)) 
B =  (np.linalg.norm(s))


# Loop through time and evolve the wave-fields using ADER time-stepping scheme of N+1 order of accuracy
start = timeit.default_timer()

# Generate initial conditions
rate.mms(v, s, U_t, V_t, U_x, V_x, y, t, type_0)

for t in utils.drange (0.0, tend+dt,dt):
    n = n+1
    
    # compute numerical solution 
    time_integrator.elastic_RK4(v, s, v, s, rho, mu, nx, dx, order, y, t, dt, r0, r1,  tau_11,\
                    tau_21, tau_12, tau_22, type_0, forcing)
 
    # Analytical solution
    rate.mms(U, V, U_t, V_t, U_x, V_x, y, t+dt, type_0)
    
    # compute error and append to the error array
    EU.append(np.linalg.norm(U-v)/A)
    EV.append(np.linalg.norm(V-s)/B)
    
    
    T.append(t)









    # Updating plots
    if n % iplot == 0: 
        for l in line1:
            l.remove()
            del l               
        for l in line2:
            l.remove()
            del l
        for l in line3:
            l.remove()
            del l               
        for l in line4:
            l.remove()
            del l 

        # Display lines
        line1 = ax1.plot(y, v, 'r', y, U, 'k--')
        ax1.legend(iter(line1),('Numerical', 'Analytical'))
        line2 = ax2.plot(y, s, 'r', y, V, 'k--')
        ax2.legend(iter(line2),('Numerical', 'Analytical'))
        line3 = ax3.plot(T, EU, 'k--')
        ax3.set_yscale("log")#, nonposx='clip')
        line4 = ax4.plot(T, EV, 'k--')
        ax4.set_yscale("log")#, nonposx='clip')
        plt.gcf().canvas.draw()
       
plt.ioff()
plt.show()

# Simulation end time
stop = timeit.default_timer()
print('total simulation time = ', stop - start)                   # print the time required for simulation
print('spatial order  of accuracy = ', order)                                  # print the polynomial degree used
print('number of grid points = ', nx)                     # print the degree of freedom
print('maximum relative error in particle velocity = ', max(EU))  # max. relative error in particle velocity
