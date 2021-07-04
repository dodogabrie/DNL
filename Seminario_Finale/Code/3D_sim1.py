# Import data
import time
import numpy as np
import plotly.graph_objects as go

nb_frames = 100
mu = np.linspace(1.2, 1.65, nb_frames) 
my_frames = []
x0, y0, z0 = np.loadtxt('data/frames/sim0.txt')
for k in range(nb_frames):
    xx, yy, zz = np.loadtxt(f'data/frames/sim{k}.txt')
    f = go.Frame(
        data=go.Scatter3d(
                x=xx,
                y=yy,
                z=zz,
                marker=go.scatter3d.Marker(size=1, color='blue'),
                opacity=0.1,
                mode='markers',
        ),
        name=str(k) # you need to name the frame for the animation to behave properly
    )

    my_frames.append(f)
fig = go.Figure(
        frames=my_frames)

# Add data to be displayed before animation starts
fig.add_trace(go.Scatter3d(
               x=x0,
               y=y0,
               z=z0,
              marker=go.scatter3d.Marker(size=1, color='blue'),
              opacity=0.1,
              mode='markers',)
            )

def frame_args(duration):
    return {
            "frame": {"duration": duration},
            "mode": "immediate",
            "fromcurrent": True,
            "transition": {"duration": duration, "easing": "linear"},
        }

sliders = [
            {
                "pad": {"b": 10, "t": 60},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": [
                    {
                        "args": [[f.name], frame_args(0)],
                        "label": 'mu = ' + str(np.round(mu[k], 3)),
                        "method": "animate",
                    }
                    for k, f in enumerate(fig.frames)
                ],
            }
        ]

# Layout
fig.update_layout(
         title='Simulazione con transiente',
         width=1200,
         height=700,
         scene=dict(
                    zaxis=dict(range=[-0.1, 6.8], autorange=False),
                    aspectratio=dict(x=1, y=1, z=1),
                    ),
         updatemenus = [
            {
                "buttons": [
                    {
                        "args": [None, frame_args(0)],
                        "label": "&#9654;", # play symbol
                        "method": "animate",
                    },
                    {
                        "args": [[None], frame_args(0)],
                        "label": "&#9724;", # pause symbol
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 70},
                "type": "buttons",
                "x": 0.1,
                "y": 0,
            }
         ],
         sliders=sliders
)

fig.show()
