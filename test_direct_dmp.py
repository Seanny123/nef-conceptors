import nengo
import numpy as np
import matplotlib.pyplot as plt

from direct_dmp_net2 import get_direct_decoders

period = 0.5
sin_per = (2 * np.pi * 10)


def target_func(t):
    return np.sin(sin_per*t)


def bump_func(t):
    return 1 if t < 0.1 else 0


with nengo.Network() as dd_model:
    osc = nengo.Network()
    osc.config[nengo.Ensemble].neuron_type = nengo.LIFRate()
    osc.config[nengo.Ensemble].seed = 0
    nengo.networks.Oscillator(0.1, 2*np.pi/period, 300, net=osc)

    output = nengo.Node(size_in=1)

    dd = get_direct_decoders(target_func, period, osc, bump_func)

    bump = nengo.Node(bump_func)
    nengo.Connection(bump, osc.ensemble[0])

    nengo.Connection(osc.ensemble.neurons, output, transform=dd.T)

    p_out = nengo.Probe(output)

with nengo.Simulator(dd_model) as dd_sim:
    dd_sim.run(3*period)

plt.plot(dd_sim.data[p_out][int(2*period/dt):])
plt.show()