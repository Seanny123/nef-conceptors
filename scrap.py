import nengo
import numpy as np
import nengo.dists
from nengo.utils.builder import default_n_eval_points

np.random.seed(3)
# maps from input value (in this case, theta) to output value

model = nengo.Network()
tau = 0.1
model.config[nengo.Ensemble].neuron_type = nengo.Direct()
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

    n_neurons = 1000
    readout = nengo.Ensemble(n_neurons=500, dimensions=1)
    nengo.Connection(osc[:2], readout, 
                     function=lambda x: np.arctan2(x[1], x[0]))

    act1 = nengo.Ensemble(n_neurons=100, dimensions=1, radius=np.sqrt(2))
    act2 = nengo.Ensemble(n_neurons=100, dimensions=1, radius=np.sqrt(2))
    act3 = nengo.Ensemble(n_neurons=100, dimensions=1, radius=np.sqrt(2))
    
    nengo.Connection(readout, act1, function=lambda x: x**2)
    nengo.Connection(readout, act2, function=lambda x: -x)
    nengo.Connection(readout, act3, function=lambda x: np.sqrt(abs(x)) * np.sign(x))
    
    weights = nengo.Node(output=[1,1,1])
    
    mult_encs = nengo.dists.Choice([[1,1],[1,-1],[-1,-1],[-1,1]])
    mult1 = nengo.Ensemble(n_neurons=500, dimensions=2,
                           encoders=mult_encs)
    mult2 = nengo.Ensemble(n_neurons=500, dimensions=2,
                           encoders=mult_encs)
    mult3 = nengo.Ensemble(n_neurons=500, dimensions=2,
                           encoders=mult_encs)

    nengo.Connection(act1, mult1[0])
    nengo.Connection(act2, mult2[0])
    nengo.Connection(act3, mult3[0])
    nengo.Connection(weights[0], mult1[1])
    nengo.Connection(weights[1], mult2[1])
    nengo.Connection(weights[2], mult3[1])
    
    output = nengo.Ensemble(n_neurons=100, dimensions=1)
    nengo.Connection(mult1, output, function=lambda x: x[0]*x[1])
    nengo.Connection(mult2, output, function=lambda x: x[0]*x[1])
    nengo.Connection(mult3, output, function=lambda x: x[0]*x[1])
