##accoustic wave 
import first_derivative_sbp_operators
import numpy as np
import boundarycondition        

def elastic_rate(
    hv,
    hs,
    v,
    s,
    rho,
    mu,
    nx,
    dx,
    order,
    t,
    y,
    r0,
    r1,
    tau0_1,
    tau0_2,
    tauN_1,
    tauN_2,
    type_0,
    forcing,
):

    ## hv,hs,y,r0,r1 is still unknown, 
    ## v is v, 
    ## s is sigma, 
    ## rho,mu is constant
    ## nx is the number of grids in space line
    ## dx is the distance betweeen consecutive space grid
    ## oreder is the order of time intergration
    ## t is the time 
    ## tau is the SAT penalty coefficients
    ## type_0 is the type of initial guess of our fuctions v and sigma.
    ## forcing is the forcing term as the name suggested.

    # we compute rates that will be used for Runge-Kutta time-stepping

    V = np.zeros((nx, 1)) ## This is the V function.
    S = np.zeros((nx, 1))  ## This is the \sigma function.
    Vt = np.zeros((nx, 1)) ## This is the LHS: \frac{d V}{d t}.
    St = np.zeros((nx, 1))  ## This is the LHS: \frac{d \sigma}{d t}.
    Vx = np.zeros((nx, 1))  ## This is the discritised spatial derivative of v, that will be obtain using SBP operators.
    Sx = np.zeros((nx, 1))  ## This is the discritised spatial derivative of \sigma, that will be obtain using SBP operators.

    mms(V, S, Vt, St, Vx, Sx, y, t, type_0) ###

    # initialize arrays for computing derivatives
    vx = np.zeros((nx, 1))
    sx = np.zeros((nx, 1))

    # compute first derivatives for velocity and stress fields
    '''
    dx() is a function to perform the matrix multiplication of D \cdot *, where * is the thing you wanna take spacial derivative.
    '''
    first_derivative_sbp_operators.dx(vx, v, nx, dx, order) 
    first_derivative_sbp_operators.dx(sx, s, nx, dx, order)

    # compute the elastic rates
    hv[:, :] = (1.0 / rho) * sx + forcing * (Vt - (1.0 / rho) * Sx)
    hs[:, :] = mu * vx + forcing * (St - mu * Vx)

    # impose boundary conditions using penalty: SAT
    impose_bc(
        hv,
        hs,
        v,
        s,
        rho,
        mu,
        nx,
        dx,
        order,
        forcing * V,
        forcing * S,
        r0,
        r1,
        tau0_1,
        tau0_2,
        tauN_1,
        tauN_2,
    )





def mms(V, S, V_t, S_t, V_x, S_x, y, t, type_0):
    '''
    This is  to initialise V,S,V_t,S_t,V_x,S_x.
    '''
    import numpy as np

    if type_0 in ("Gaussian"):

        delta = 0.015 * (y[-1, 0] - y[0, 0])

        cs = 3.464
        rho = 2.6702

        Zs = rho * cs

        x0 = 0.5 * (y[-1, 0] - y[0, 0])

        V[:, :] = (
            1
            / np.sqrt(2.0 * np.pi * delta ** 2)
            * 0.5
            * (
                np.exp(-(y + cs * (t) - x0) ** 2 / (2.0 * delta ** 2))
                + np.exp(-(y - cs * (t) - x0) ** 2 / (2.0 * delta ** 2))
            )
        )
        S[:, :] = (
            1
            / np.sqrt(2.0 * np.pi * delta ** 2)
            * 0.5
            * Zs
            * (
                np.exp(-(y + cs * (t) - x0) ** 2 / (2.0 * delta ** 2))
                - np.exp(-(y - cs * (t) - x0) ** 2 / (2.0 * delta ** 2))
            )
        )

        V_t[:, :] = 0
        S_t[:, :] = 0
        V_x[:, :] = 0
        S_x[:, :] = 0

    if type_0 in ("Sinusoidal"):

        delta = y[-1, 0] - y[0, 0]

        ny = 20.5 / delta * np.pi
        nt = 2.5 * np.pi
        fs = 9.33

        V[:, :] = np.cos(nt * t) * np.sin(ny * y + fs)

        S[:, :] = ny * np.sin(nt * t) * np.cos(ny * y - fs)

        V_t[:, :] = -nt * np.sin(nt * t) * np.sin(ny * y + fs)
        S_t[:, :] = nt * ny * np.cos(nt * t) * np.cos(ny * y - fs)

        V_x[:, :] = ny * np.cos(nt * t) * np.cos(ny * y + fs)
        S_x[:, :] = -ny * ny * np.sin(nt * t) * np.sin(ny * y - fs)




# To update penalty term h11 based on dx and the
# required order

def penaltyweight(h11, dx, order):
    '''
    This is to calculate the penalyty weight, and it is done by updating the h11, with given dx and order
    '''
    if order == 2:
        h11[:] = 0.5 * dx
    if order == 4:
        h11[:] = (17.0 / 48.0) * dx

    if order == 6:
        h11[:] = 13649.0 / 43200.0 * dx

def impose_bc(
    hv,
    hs,
    v,
    s,
    rho,
    mu,
    nx,
    dx,
    order,
    V,
    S,
    r0,
    r1,
    tau0_1,
    tau0_2,
    tauN_1,
    tauN_2,
):
    '''
    Impose bc on to hv and hs, which in this content can be think of the H^-1 bit and the forcing bit...
    So dv/dt (LHS term) = Dv + hv (penalty term + forcing term)...
    '''
    # impose boundary conditions
    import numpy as np
    import boundarycondition

    # penalty weights
    h11 = np.zeros((1, 1))
    penaltyweight(h11, dx, order)

    mv = np.zeros((1, 1)) ## mu
    ms = np.zeros((1, 1))

    pv = np.zeros((1, 1))
    ps = np.zeros((1, 1))

    v0 = v[0, :]
    s0 = s[0, :]

    vn = v[nx - 1, :]
    sn = s[nx - 1, :]

    # boundary forcing
    V0 = V[0, :] ## These are the bondary stencile 
    S0 = S[0, :]

    Vn = V[nx - 1, :]
    Sn = S[nx - 1, :]

    # compute SAT terms
    boundarycondition.bcm(mv, ms, v0, s0, V0, S0, rho, mu, r0)
    boundarycondition.bcp(pv, ps, vn, sn, Vn, Sn, rho, mu, r1)

    # penalize boundaries with the SAT terms
    hv[0, :] = hv[0, :] - tau0_1 / h11 * mv
    hs[0, :] = hs[0, :] - tau0_2 / h11 * ms

    hv[nx - 1, :] = hv[nx - 1, :] - tauN_1 / h11 * pv

# def advection_rate(hv, v, nx, dx, order, t, y, tau):
'''
Not relevant, this is for aadvection
'''
#     # we compute rates that will be used for Runge-Kutta time-stepping
#     #
#     import first_derivative_sbp_operators
#     import numpy as np

#     # initialize arrays for computing derivatives
#     vx = np.zeros((nx, 1))

#     # compute first derivatives of the advected field v
#     first_derivative_sbp_operators.dx(vx, v, nx, dx, order)

#     # compute rates
#     hv[:, :] = -vx

#     # impose boundary conditions using penalty: SAT

#     # penalty weights
#     h11 = np.zeros((1, 1))
#     penaltyweight(h11, dx, order)

#     V0 = np.zeros((1, 1))
#     # boundary forcing
#     g(V0, t)

#     # penalize boundaries with the SAT terms (Change)
#     hv[0, :] = hv[0, :] - tau / h11 * (v[0, :] - V0)
        
# def g(V, t):
'''
Not used

'''
#     import numpy as np

#     V[:, :] = 0.0

#     if t <= 1.0 and t >= 0.0:
#         V[:, :] = (np.sin(np.pi * t)) ** 4