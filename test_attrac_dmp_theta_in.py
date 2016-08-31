import nengo
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

from attractor_dmp_net import make_attrac_net
from direct_dmp_net2 import get_direct_decoders
from constants import *

period = 0.5
sin_per = (2 * np.pi * 10)
seed = 0


def target_func(t):
    return np.sin(sin_per*t)


def bump_func(t):
    return 1 if t < 0.1 else 0

pre_dat = target_func(np.linspace(0, period, 100))
xv = np.linspace(-np.pi, np.pi, pre_dat.shape[0])
proc_func = interpolate.interp1d(xv, pre_dat)

xv = np.linspace(0, 2*np.pi, 100)
inter_arc = interpolate.interp1d(xv, np.arctan2(np.cos(xv), np.sin(xv)))


def arc_func(x):
    return inter_arc(2*np.pi*x/period % (2*np.pi))

with nengo.Network() as ad_model:
    bump = nengo.Node(bump_func)

    osc = nengo.Network()
    osc.config[nengo.Ensemble].neuron_type = nengo.LIFRate()
    osc.config[nengo.Ensemble].seed = seed
    nengo.networks.Oscillator(0.1, 2*np.pi/period, 300, net=osc)
    dd = get_direct_decoders(arc_func, period, osc, bump_func)

    dmp, conn_func = make_attrac_net(proc_func, 300, dd=dd, seed=seed)

    nengo.Connection(bump, osc.ensemble[0])
    nengo.Connection(osc.ensemble.neurons, dmp.input, function=conn_func)

    p_arc = nengo.Probe(arctan, synapse=0.01)

    p_out = nengo.Probe(dmp.output, synapse=0.01)

with nengo.Simulator(ad_model) as ad_sim:
    ad_sim.run(4*period)

plt.plot(ad_sim.data[p_out][int(2*period/dt):])
plt.show()
