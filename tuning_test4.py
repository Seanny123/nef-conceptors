from tanh_neuron import TanhWithBias

import numpy as np
import nengo
import ipdb
import matplotlib.pyplot as plt

SEED = 0
n_neurons = 5
neuron_type = TanhWithBias(seed=SEED)
sig_dims = 1


with nengo.Network() as model:
    in_nd = nengo.Node(lambda t: np.ones(n_neurons) * (1.0 - t))
    # should use np.eye if you actually want orthogonal stuff
    sig_reserv = nengo.Ensemble(n_neurons, n_neurons, neuron_type=neuron_type, seed=SEED)
    nengo.Connection(in_nd, sig_reserv, synapse=None, seed=SEED)

    p_rate = nengo.Probe(sig_reserv.neurons, synapse=None)
    p_in = nengo.Probe(in_nd)

with nengo.Simulator(model) as conc_sim:
    conc_sim.run(2)

plt.plot(conc_sim.data[p_rate])
plt.plot(conc_sim.data[p_in])
plt.show()

ipdb.set_trace()
