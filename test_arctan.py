import nengo
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import ipdb

from direct_dmp_net2 import get_direct_decoders
from constants import *

period = 0.5
seed = 0


def bump_func(t):
    return 1 if t < 0.1 else 0

xv = np.linspace(0, 2*np.pi, 100)
inter_arc = interpolate.interp1d(xv, np.arctan2(np.cos(xv), np.sin(xv)))


def arc_func(x):
    return inter_arc(2*np.pi*x/period % (2*np.pi))

test_val = np.linspace(0, period, 100)
plt.plot(test_val, arc_func(test_val))
plt.show()

with nengo.Network() as ad_model:
    bump = nengo.Node(bump_func)

    osc = nengo.Network()
    osc.config[nengo.Ensemble].neuron_type = nengo.LIFRate()
    osc.config[nengo.Ensemble].seed = 0
    nengo.networks.Oscillator(0.1, 2*np.pi/period, 300, net=osc)

    arctan = nengo.Ensemble(250, 1, radius=np.pi)

    dd = get_direct_decoders(arc_func, period, osc, bump_func)

    nengo.Connection(bump, osc.ensemble[0])
    nengo.Connection(osc.ensemble.neurons, arctan, transform=dd.T)

    p_arc = nengo.Probe(arctan, synapse=0.01)

with nengo.Simulator(ad_model) as ad_sim:
    ad_sim.run(4*period)

ipdb.set_trace()
gd = ad_sim.data[p_arc][int(2*period/dt):]
plt.plot(gd)
plt.show()
