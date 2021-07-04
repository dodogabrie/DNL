import numpy as np
import os
import importlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from numba import jit
import nolds
import sys
sys.path.insert(0,'./fortran_package/') # carico la directory con i pacchetti
import maps
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import plotly.graph_objects as go
import plotly.io as pio
import ComputeLE

def confronto(my_method, show = False):
    """
    Main per integrare con il metodo my_method (da passare a solve_ivp).
    Restitiusce il plot (da gestire esternamente con plt.show(), ad esempio).
    """
    n = 200000
    t = np.linspace(0, 1000, n)
    t1 = np.linspace(0, 5000, n)
    # Parameters
    a = np.array([[1,    1.09, 1.52, 0   ],
                  [0,    1,    0.44, 1.36],
                  [2.33, 0,    1,    0.47],
                  [1.21, 0.51, 0.35, 1   ]])
    r = np.array([1, 0.72,     1.53, 1.27])
    param = [a, r]
    
    # Inital Values
    init = np.array([0.1, 0.1, 0.1, 0.1])
    x = np.ones((n, 4))*init
    x = ComputeLE.motion(t, x, a, r)
    x[0] = x[-1]
    x = ComputeLE.motion(t1, x, a, r)
    xx = x[:,0]
    yy = x[:,1]
    zz = x[:,2]
    ww = x[:,3]
    #plt.plot(xx[:100], alpha=0.4)
    if show:
        marker_data = go.Scatter3d(
                        x=xx[::5],
                        y=yy[::5], 
                        z=zz[::5], 
                        marker=go.scatter3d.Marker(size=1, color='blue'), 
                        opacity=0.1, 
                        mode='markers',)
        fig=go.Figure(data=marker_data)
        fig.show()
confronto('DOP853', show=True)
