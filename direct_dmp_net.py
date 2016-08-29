import nengo
import numpy as np
from nengo.utils.connection import target_function

from constants import *


def make_dd_net(target_func, period, n_neurons=100, seed=0, label=None):

    t_val = np.arange(0, period, dt)

    with nengo.Network(label=label) as dd_net:
        dd_net.input = nengo.Node(size_in=2)
        # ideally this should be using a config object
        readout = nengo.Ensemble(n_neurons=n_neurons, dimensions=2, neuron_type=nengo.LIFRate(), seed=seed)
        nengo.Connection(dd_net.input, readout, synapse=None)

        dd_net.output = nengo.Node(size_in=1)
        nengo.Connection(readout, dd_net.output,
                         **target_function(
                             np.array([np.cos(2*np.pi*t_val/period), np.sin(2*np.pi*t_val/period)]).T,
                             target_func(t_val)))

    return dd_net
