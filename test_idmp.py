import nengo
import numpy as np

from idmp import InverseDMP


def heart_func(x):
    theta = np.arctan2(x[1], x[0])
    r = 2 - 2 * np.sin(theta) + np.sin(theta)*np.sqrt(np.abs(np.cos(theta)))/(np.sin(theta)+1.4)
    return -r*np.cos(theta), r*np.sin(theta)

def negheart_func(x):
    theta = np.arctan2(x[1], x[0])
    r = 1 - 1 * np.sin(theta) + np.sin(theta)*np.sqrt(np.abs(np.cos(theta)))/(np.sin(theta)+1.4)
    return r*np.cos(theta), -r*np.sin(theta)

def figure_8(x):
    return x[0]*x[1], x[1]

# create an oscillator and a figure-8 oscillator
# see if the position decoded from them is the same

# create a heart oscillator with a different oscillator start
# check that it's position is different

# see if forcing functions are truly necessary

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

    heart = nengo.Ensemble(200, 2)
    nengo.Connection(osc[:2], heart, function=heart_func)

    neg_heart = nengo.Ensemble(200, 2)
    nengo.Connection(osc[:2], neg_heart, function=negheart_func)

    sec_osc = nengo.Ensemble(n_neurons=1500, dimensions=3, radius=1.4)

    sec_bump = nengo.Node(lambda t: 1 if t < 0.55 else 0)
    nengo.Connection(sec_bump, sec_osc[0])
    nengo.Connection(rate, sec_osc[2], transform=[-1])
    nengo.Connection(sec_osc, sec_osc[:2], synapse=tau, function=cycle)

    fig8 = nengo.Ensemble(200, 2)
    nengo.Connection(sec_osc[:2], fig8, function=figure_8)

