import nengo
import numpy as np
import matplotlib.pyplot as plt
import ipdb

from direct_dmp_net2 import get_direct_decoders
from constants import *


def bump_func(t):
    return 1 if t < 0.1 else 0

period = 0.5
t_blend = 1000 * dt  # same as Matlab
t_start = period * 3


class Switch:
    def __init__(self):
        self.return_val = 1
        self.dec = 1 / t_blend * dt

    def step(self, t):
        """start at 1 and go to 0"""
        if t < t_start:
            return 1
        elif self.return_val > 0:
                self.return_val -= self.dec
        return self.return_val

sw = Switch()


def scalar_blend(sig_a, sig_b, seed=0):
    with nengo.Network() as dd_model:
        osc = nengo.Network()
        osc.config[nengo.Ensemble].neuron_type = nengo.LIFRate()
        osc.config[nengo.Ensemble].seed = seed
        nengo.networks.Oscillator(0.1, 2*np.pi/period, 300, net=osc)

        n_prod = 250
        prod_dim = 2
        mult_config = nengo.Config(nengo.Ensemble)
        mult_config[nengo.Ensemble].encoders = nengo.dists.Choice([[1, 1], [1, -1], [-1, -1], [-1, 1]])
        mult_config[nengo.Ensemble].neuron_type = nengo.LIFRate()
        mult_config[nengo.Ensemble].radius = np.sqrt(2)

        with mult_config:
            prod_a = nengo.Ensemble(n_prod, prod_dim)
            prod_b = nengo.Ensemble(n_prod, prod_dim)

        scale_control = nengo.Node(sw.step)
        in_start = nengo.Node(lambda t, x: 1 - x, size_in=1)  # start at 0 and go to 1
        nengo.Connection(scale_control, in_start, synapse=None)

        output = nengo.Node(size_in=1)

        dda = get_direct_decoders(sig_a, period, osc, bump_func)
        ddb = get_direct_decoders(sig_b, period, osc, bump_func)

        bump = nengo.Node(bump_func)
        nengo.Connection(bump, osc.ensemble[0])

        nengo.Connection(osc.ensemble.neurons, prod_a[0], transform=dda.T)
        nengo.Connection(osc.ensemble.neurons, prod_b[0], transform=ddb.T)

        nengo.Connection(in_start, prod_a[1])
        nengo.Connection(scale_control, prod_b[1])

        conn_a = nengo.Connection(prod_a, output, function=lambda x: x[0]*x[1])
        conn_b = nengo.Connection(prod_b, output, function=lambda x: x[0]*x[1])

        p_in_start = nengo.Probe(in_start)
        p_inhib = nengo.Probe(scale_control)
        p_a = nengo.Probe(conn_a, synapse=0.01)
        p_b = nengo.Probe(conn_b, synapse=0.01)
        p_out = nengo.Probe(output)

    with nengo.Simulator(dd_model) as dd_sim:
        dd_sim.run((t_start + t_blend + 2*period))

    return dd_sim.data[p_a], dd_sim.data[p_b], dd_sim.data[p_out], dd_sim.data[p_inhib], dd_sim.data[p_in_start]

# direct DMPs also have a frequency limit, but it's not as severe as attractor DMPs
cos_per = (2 * np.pi * 5) / period
sin_per = (2 * np.pi * 2) / period


def sin_targ(t):
    return np.sin(sin_per*t)


def cos_targ(t):
    return 0.5 * np.cos(cos_per*t)

res = scalar_blend(sin_targ, cos_targ)

plt.plot(res[0][int(2*period/dt):])
plt.plot(res[1][int(2*period/dt):])
plt.show()

plt.plot(res[2][int(2*period/dt):])
plt.show()
ipdb.set_trace()
