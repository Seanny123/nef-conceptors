import nengo
import numpy as np

from idmp import InverseDMP
from preprocess import *

# load some sigs
# see if Eric's point eval thing can save me some time
# YES IT CAN, BUT I CAN'T FIGURE IT OUT

# create two oscillators with the same start point
# see if the position decoded from them is the same

# create a shifted oscillator or at least one with a different start
# check that it's position is different

# see if forcing functions are truly necessary


dt = 0.001
trange = (0, np.round(2*np.pi, decimals=3))
filename = "processed/out_61_pn_2.npz"

model = nengo.Network()
tau = 0.1

with model:
    osc = nengo.Ensemble(n_neurons=1500, dimensions=3, radius=1.4)

    def cycle(x):
        """makes a speed controlled oscillator"""
        a = 1.0
        b = 2.0 * np.pi * x[2]
        r = np.sqrt(x[0]**2.0 + x[1]**2.0)
        theta = np.arctan2(x[1], x[0])
        dr = 10.0*(-r**3.0 + a*r)
        dtheta = b

        dx = dr*np.cos(theta) - r*np.sin(theta)*dtheta
        dy = dr*np.sin(theta) + r*np.cos(theta)*dtheta

        return [x[0] + tau*dx, x[1] + tau*dy]

    nengo.Connection(osc, osc[:2], synapse=tau, function=cycle)

    rate = nengo.Node([1])
    nengo.Connection(rate, osc[2])

    bump = nengo.Node(lambda t: 1 if t < 0.5 else 0)
    nengo.Connection(bump, osc[0])

    ramp = nengo.Ensemble(200, 1)
    nengo.Connection(osc[:2], ramp,
        function=lambda x: np.arctan2(x[0], x[1]))

    ramp_idmp = InverseDMP(600, lambda x: np.array([x]))
    nengo.Connection(ramp, ramp_idmp.input)

    out_ramp = nengo.Node(size_in=1)
    nengo.Connection(ramp_idmp.state[0], out_ramp)

    sin_ramp = nengo.Node(size_in=1)
    nengo.Connection(ramp, sin_ramp,
        function=lambda x: x * np.sin(x))

    # none type stuff seems to happen after creating an object and manipulating it?
    # sin ramp is definitely not working, GET HELP
    sin_idmp = InverseDMP(600, lambda x: x * np.sin([x]))
    nengo.Connection(sin_ramp, sin_idmp.input)

    out_sin = nengo.Node(size_in=1)
    nengo.Connection(sin_idmp.state[0], out_sin)

    sec_osc = nengo.Ensemble(n_neurons=1500, dimensions=3, radius=1.4)

    sec_bump = nengo.Node(lambda t: 1 if t < 0.55 else 0)
    nengo.Connection(sec_bump, sec_osc[0])
    nengo.Connection(rate, sec_osc[2], transform=[-1])
    nengo.Connection(sec_osc, sec_osc[:2], synapse=tau, function=cycle)

    del_ramp = nengo.Node(size_in=1)
    nengo.Connection(sec_osc[:2], del_ramp,
        function=lambda x: np.arctan2(x[0], x[1]))
