import nengo
import numpy as np
import nengo.dists
from nengo.utils.builder import default_n_eval_points
import ipdb

T = 1.0
dt = 0.001

t = np.arange(1000)*dt

pattern1 = 2 * t -1
pattern2 = np.sin(t*2*np.pi)
pattern3 = np.where(t<0.5, 2*t, 2-2*t)
pattern4 = -pattern3

def make_pattern():
    sig = nengo.processes.WhiteSignal(T, high=2)
    # this is returning a flat signal and I don't know why
    return sig.run(1000)

np.random.seed(3)
#patterns = [make_pattern()[:,0] for i in range(4)]
patterns = [pattern1, pattern2, pattern3, pattern4]
# maps from input value (in this case, theta) to output value
patterns = np.array(patterns)

model = nengo.Network()
#model.config[nengo.Ensemble].neuron_type = nengo.Direct()
#model.config[nengo.Connection].solver = nengo.solvers.LstsqL2(reg=0.001)
tau = 0.1
omega = 2*np.pi
with model:
    system = nengo.Ensemble(n_neurons=1000, dimensions=3, radius=1.4, label='cycle')

    def cycle(x):
        """makes a speed controlled oscillator"""
        a = 1
        b = 2 * np.pi * x[2]
        r = np.sqrt(x[0]**2 + x[1]**2)
        theta = np.arctan2(x[1], x[0])
        dr = 10*(-r**3 + a*r)
        dtheta = b

        dx = dr*np.cos(theta) - r*np.sin(theta)*dtheta
        dy = dr*np.sin(theta) + r*np.cos(theta)*dtheta

        return [x[0] + tau*dx, x[1] + tau*dy]

    nengo.Connection(system, system[:2], synapse=tau, function=cycle)

    rate = nengo.Node([1])
    nengo.Connection(rate, system[2])

    bump = nengo.Node(lambda t: 1 if t < 0.5 else 0)
    nengo.Connection(bump, system[0])

    n_neurons = 1000
    n_points = default_n_eval_points(n_neurons, 2)
    e_p1 = nengo.dists.UniformHypersphere(surface=True).sample(n_points, 2)
    n_points = default_n_eval_points(n_neurons, 1)
    e_p2 = nengo.dists.Uniform(0, 1).sample(n_points, 1)
    readout = nengo.Ensemble(n_neurons=n_neurons, dimensions=3,
                             eval_points=np.concatenate((e_p1, e_p2), axis=1),
                             label="readout")
    def read_func(x):
        """based off angle, read out pattern"""
        theta_p = np.arctan2(x[1], x[0])
        r_p = np.sqrt(x[0]**2 + x[1]**2)
        if theta_p < 0: 
            theta_p += 2 * np.pi
        index_p = int(len(patterns[0]) * theta_p / (2 * np.pi))
        #ipdb.set_trace()
        return patterns[int(np.round(x[2]))][index_p] * r_p

    nengo.Connection(system[:2], readout[:2])

    switch = nengo.Node([0])
    nengo.Connection(switch, readout[2])

    output = nengo.Ensemble(n_neurons=30, dimensions=1, label="output")
    nengo.Connection(readout, output, function=read_func)