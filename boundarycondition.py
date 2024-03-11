import numpy as np


def bcm(mv, ms, v, s, V, S, rho, mu, r): #K):
    # Z = np.sqrt(rho*K)
    # #########################
    # #########################
    # This is how we use mv and ms
                        # hv[0, :] = hv[0, :] - tau0_1 / h11 * mv
                        # hs[0, :] = hs[0, :] - tau0_2 / h11 * ms


    cs = np.sqrt(mu / rho)
    zs = rho * cs

    p = 0.5 * (zs * v - s)
    q = 0.5 * (zs * v + s)

    P = 0.5 * (zs * V - S)
    Q = 0.5 * (zs * V + S)
    g = P - r * Q

    mv[:] = 1.0 / rho * ((p - r * q) - (P - r * Q))
    ms[:] = -mu / zs * ((p - r * q) - (P - r * Q))


def bcp(pv, ps, v, s, V, S, rho, mu, r):
    ### hv[nx - 1, :] = hv[nx - 1, :] - tauN_1 / h11 * pv

    cs = np.sqrt(mu / rho)
    zs = rho * cs

    p = 0.5 * (zs * v + s)
    q = 0.5 * (zs * v - s)

    P = 0.5 * (zs * V + S)
    Q = 0.5 * (zs * V - S)

    g = P - r * Q

    pv[:] = 1.0 / rho * ((p - r * q) - (P - r * Q))
    ps[:] = -mu / zs * ((p - r * q) - (P - r * Q))
