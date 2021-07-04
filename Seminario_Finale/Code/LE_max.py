import numpy as np
import plotly.express as px
x = np.loadtxt('data/MaxLE.txt')
fig = px.histogram(x, nbins=1000)
fig.update_layout(
    title="Distribuzione del massimo esponente di Lyapunov",
    xaxis_title="Max LE",
    yaxis_title="Count",
)
fig.show()
