import nengo
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

from attractor_dmp_net import make_attrac_net
from constants import *
from process import d3_scale

period = 0.5
sin_per = (2 * np.pi * 10)


def target_func(t):
    return np.sin(sin_per*t)


def bump_func(t):
    return 1 if t < 0.1 else 0

pre_dat = target_func(np.linspace(0, period, 100))
xv = np.linspace(-np.pi, np.pi, pre_dat.shape[0])
proc_func = interpolate.interp1d(xv, pre_dat)

with nengo.Network() as ad_model:
    bump = nengo.Node(bump_func)

    osc = nengo.Network()
    osc.config[nengo.Ensemble].neuron_type = nengo.LIFRate()
    osc.config[nengo.Ensemble].seed = 0
    nengo.networks.Oscillator(0.1, 2*np.pi/period, 300, net=osc)

    dmp, conn_func = make_attrac_net(proc_func, 300)

    nengo.Connection(bump, osc.ensemble[0])
    nengo.Connection(osc.ensemble, dmp.input, function=conn_func)

    p_out = nengo.Probe(dmp.output, synapse=0.01)

with nengo.Simulator(ad_model) as ad_sim:
    ad_sim.run(4*period)

g_dat = ad_sim.data[p_out][int(2*period/dt):]
plt.plot(d3_scale(g_dat))
plt.show()

