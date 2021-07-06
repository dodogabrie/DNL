"""
Calculate the Lyapunov exponents for a set of ODEs using the method described 
in Sandri (1996), through the use of the variational matrix.
"""

import numpy as np
from numba import jit, njit, prange

@njit
def complete_motion(t, t_inner, x, a, r):
    #Ssol_temp = np.append(x[: , 0], Phi.flatten())
    for i,(t1,t2) in enumerate(zip(t[:-1], t[1:])):
        x[i+1] = x[i] + RK4(func, x[i], t1, t2, a, r)
    x0 = x[-1]
    D = 4
    Phi = np.eye(D, dtype=np.float64).flatten() # Identità di tipo float
    Ssol = np.append(x0, Phi.flatten())
    for i,(t1,t2) in enumerate(zip(t_inner[:-1], t_inner[1:])):
        Ssol_temp = Ssol + RK4(dSdt, Ssol, t1, t2, a, r)
        Ssol = Ssol_temp
        
    return x, np.reshape(Ssol[D:], (D,D))

@njit
def motion(t, x, a, r):
    for i,(t1,t2) in enumerate(zip(t[:-1], t[1:])):
        x[i+1] = x[i] + RK4(func, x[i], t1, t2, a, r)
    return x

@njit
def func(t, val, a, r):
    """
    Questa funzione definisce gli incrementi secondo il formato utile a 'solve_ivp' di scykit_learn.
    """
    ss = [np.sum(val * i) for i in a]
    diff = [r[i] * val[i] * (1 - ss[i] ) for i in range(len(r))]
    return np.array(diff)

@njit
def dSdt(t, S, a, r):
    """
    Differential equations for combined state/variational matrix
    propagation. This combined state is called S.
    """
    D = 4
    x = S[:D]
    Phi = S[D:]
    rPhi = np.reshape(Phi, (D, D))
    rdPhi = np.dot(fjac(t, x, a, r), rPhi)
    return np.append(func(t,x,a, r), rdPhi.flatten())

@njit
def fjac(t, val, a, r):
    s = [np.sum(val * a[i]) for i in range(len(a))]
    J = np.empty(np.shape(a))
    for i in range(len(r)):
        for j in range(len(r)):
            if i == j: 
                J[i][j] = r[i] * (1 - s[i]) - r[i] * val[i] * a[i][j]
            else:
                J[i][j] = - r[i] * val[i] * a[i][j]
    return J


@njit
def RK4(f, x, t1, t2, a, r):
    """
    Fourth-order, 4-step RK routine.
    Returns the step, i.e. approximation to the integral.
    If x is defined at time t_1, then stim should be an array of
    stimulus values at times t_1, (t_1+t_2)/2, and t_2 (i.e. at t1 and t2, as
    well as at the midpoint).
    Alternatively, stim may be a function pointer.
    """
    tmid = (t1 + t2)/2
    dt = t2 - t1

    K1 = f(t1, x, a, r)
    K2 = f(tmid, x + dt*K1/2, a, r)
    K3 = f(tmid, x + dt*K2/2, a, r)
    K4 = f(t2, x + dt*K3, a, r)

    return dt * (K1/2 + K2 + K3 + K4/2) / 3

@jit
def break_cond(f, xi, t1, t2, a, r, lim_dead):
    c1 = np.min(xi) <= lim_dead 
    #c2 = np.max(RK4(f, xi, t1, t2, a, r)) < 0.0000001
    #c3 = np.sum(xi) > 1.3
    return c1# or c2 or c3

@njit
def computeLE(x0, t, a, r, ttrans=None, lim_dead = 0.00001, early_stopping = False):
    """
    Computes the global Lyapunov exponents for a set of ODEs.
    f - ODE function. Must take arguments like f(t, x, p) where x and t are 
        the state and time *now*, and p is a tuple of parameters. If there are 
        no model paramters, p should be set to the empty tuple.
    x0 - Initial position for calculation. Integration of transients will begin 
         from this point.
    t - Array of times over which to calculate LE.
    p - (optional) Tuple of model parameters for f.
    fjac - Jacobian of f.
    ttrans - (optional) Times over which to integrate transient behavior.
             If not specified, assumes trajectory is on the attractor.
    method - (optional) Integration method to be used by scipy.integrate.ode.
    """
    block = False
    D = len(x0)
    N = len(t)
    if ttrans is not None:
        Ntrans = len(ttrans)
    dt = t[1] - t[0]

    # integrate transient behavior
    Phi0 = np.eye(D, dtype=np.float64).flatten() # Identità di tipo float

    if ttrans is not None:
        xi = x0
        for i,(t1,t2) in enumerate(zip(ttrans[:-1], ttrans[1:])):
            xip1 = xi + RK4(func, xi, t1, t2, a, r)
            xi = xip1
            if break_cond(func, xi, t1, t2, a, r, lim_dead):
                block = True
                return np.zeros((len(t)-1, 4)), block
        x0 = xi
    # start LE calculation
    if early_stopping: 
        early = 0
    else:
        early = N-1
    temp = 0
    L1 = np.inf
    LE = np.zeros((N-1, D), dtype=np.float64)
    final_LE = np.zeros((N-1, D), dtype=np.float64)
    LE_aux = np.zeros((N-1, D), dtype=np.float64)
    Ssol = np.zeros((N, D*(D+1)), dtype=np.float64)
    Ssol[0] = np.append(x0, Phi0)
    for i,(t1,t2) in enumerate(zip(t[:-1], t[1:])):
        Ssol_temp = Ssol[i] + RK4(dSdt, Ssol[i], t1, t2, a, r)

        if break_cond(func, Ssol_temp[:D], t1, t2, a, r, lim_dead):
            block = True
            return np.zeros((len(t)-1, 4)), block 

        # perform QR decomposition on Phi
        rPhi = np.reshape(Ssol_temp[D:], (D, D))
        Q,R = np.linalg.qr(rPhi)
        Ssol[i+1] = np.append(Ssol_temp[:D], Q.flatten())
        LE[i] = np.abs(np.diag(R))
        logLE = np.log(LE[i])
        if i > 0:
            LE_aux[i, :] = LE_aux[i-1, :] + logLE
        else:
            LE_aux[i, :] = logLE
        final_LE[i] = LE_aux[i]/(t2)
        if early_stopping:
            early = early + 1
            if np.abs(np.max(final_LE[i]) - L1) < 0.0001:
                temp = temp + 1
            else: 
                temp = 0
            if temp >= N/5:
                return final_LE[:early], block
            L1 = np.max(final_LE[i])
    return final_LE, block
               

@njit
def update_init(a, r, radi):
    rr = np.random.rand(4,4)
    a = rr - (np.diag(rr))*np.eye(4)
    r = np.zeros(4)
    r[1:] = np.random.rand(3)
    return a/radi, r/radi

@njit
def patience_gestion(pat, patience, a, r, a_aux, r_aux):
    if pat >= patience:
        rr = np.random.rand(4,4)
        a = rr - (np.diag(rr) - 1 )*np.eye(4)
        rr = np.random.rand(3)
        r[0] = 1
        r[1:] = rr
        pat = 0
        a_aux = a
        r_aux = r
    return pat, a, r, a_aux, r_aux

@njit
def store_val(LE, dead, save, a_save, radi, r_save, a, r, a_aux, r_aux, LE_save, pat):
    if not dead:
        if np.max(LE[-1]) > np.max(LE_save):
            a_save = a
            r_save = r
            LE_save = LE
            save = save + 1
            radi = radi + radi
            a_aux = a_save
            r_aux = r_save
            pat = 0
            #print(LE_save[-1])
        else:
            a = a_aux
            r = r_aux
            radi = max(2, radi/2)
    else: 
        a = a_aux
        r = r_aux
        radi = max(2, radi/2)
    return save, a_save, r_save, a, r, a_aux, r_aux, LE_save, radi, pat

@njit
def find_best(x0, t, ttrans=None, n = 10, seed = None, patience = 10, early_stopping = False):
    np.random.seed(seed)
    rr = np.random.rand(4,4)
    a = rr - (np.diag(rr) - 1 )*np.eye(4)
    r = np.zeros(4)
    rr = np.random.rand(3)
    r[0] = 1
    r[1:] = rr
    a_save = a
    r_save = r
    a_aux = a
    r_aux = r
    LE_save = -1*np.ones((len(t)-1, 4))
    pat = 0
    save = 0
    tot_dead = 0
    radi = 2
    for i in range(n):
        pat, a, r, a_aux, r_aux = patience_gestion(pat, patience, a, r, a_aux, r_aux)
        LE, dead = computeLE(x0, t, a, r, ttrans=ttrans, early_stopping = early_stopping)
        save, a_save, r_save, a, r, a_aux, r_aux, LE_save, radi, pat = store_val(LE, dead, save, a_save, radi,
                                                                      r_save, a, r, a_aux, r_aux, LE_save, pat)
        da, dr = update_init(a, r, radi)
        if pat == 0:
            radi = 2
        a = a + da
        r = r + dr
        tot_dead = tot_dead + dead 
        pat = pat + 1
    return a_save, r_save, LE_save, save, tot_dead

@njit
def maximal_le(x0, t, ttrans=None, n = 10, seed = None, early_stopping = False):
    np.random.seed(seed)
    LE_save = np.empty(n)
    tot_dead = 0
    i = 0
    while i < n:
        rr = np.random.rand(4,4)
        a = rr - (np.diag(rr) - 1 )*np.eye(4)
        r = np.zeros(4)
        rr = np.random.rand(3)
        r[0] = 1
        r[1:] = rr
        
        LE, dead = computeLE(x0, t, a, r, ttrans=ttrans, early_stopping = early_stopping)
        da, dr = update_init(a, r, radi=1)
        tot_dead = tot_dead + dead
        if not dead:
            LE_save[i] = np.max(LE[-1])
            i = i + 1
    return LE_save, tot_dead