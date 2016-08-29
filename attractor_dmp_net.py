import nengo
import numpy as np
from nengo.utils.connection import target_function

from constants import *


def make_attrac_net(target_func, period, n_neurons=100, seed=0, label=None):

    t_val = np.arange(0, period, dt)

    with nengo.Network(label=label) as ad_net:
        ad_net.input = nengo.Node(size_in=1)
        # ideally this should be using a config object
        readout = nengo.Ensemble(n_neurons=n_neurons, dimensions=2, neuron_type=nengo.LIFRate(), seed=seed)
        nengo.Connection(ad_net.input, readout, synapse=None)

        ad_net.output = nengo.Node(size_in=1)
        nengo.Connection(readout, ad_net.output,
                         **target_function(np.array([np.cos(t_val), np.sin(t_val)]).T, target_func(t_val)))

    return ad_net
