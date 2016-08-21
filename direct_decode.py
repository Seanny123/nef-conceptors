import nengo
import numpy as np

model = nengo.Network()
tau = 0.1
ENS_SEED = 1
#model.config[nengo.Ensemble].neuron_type = nengo.Direct()

with model:
    sig = nengo.Node(lambda t: (np.cos(t), np.sin(t)))
    readout = nengo.Ensemble(n_neurons=100, dimensions=2, seed=ENS_SEED)
    nengo.Connection(sig, readout)
    p_spikes = nengo.Probe(readout.neurons, synapse=0.01)

target = 
solver = nengo.solvers.LstsqL2(reg=0.02)
print("getting decoders")
# Try a communication channel
decoders, info = solver(sim_train.data[p_spikes], target)
print("got decoders")


test_model = nengo.Network()
with test_model:
    osc = nengo.Ensemble(n_neurons=1, dimensions=3, neuron_type=nengo.Direct())

    def cycle(x):
        '''makes a speed controlled oscillator'''
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

    readout = nengo.Ensemble(n_neurons=100, dimensions=2, seed=ENS_SEED)
    nengo.Connection(osc[:2], readout)

    output = nengo.Node(size_in=1)
    p_out = nengo.Probe(output, synapse=0.01)