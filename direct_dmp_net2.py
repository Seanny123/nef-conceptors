import nengo
import numpy as np

from constants import *


def get_direct_decoders(target_func, period, arg_osc, bump_func, n_neurons=100, seed=0, label=None):

    with nengo.Network(label=label) as train_net:
        bump = nengo.Node(bump_func)
        osc = arg_osc.copy()
        nengo.Connection(bump, osc.ensemble[0])

        p_neur = nengo.Probe(osc.ensemble.neurons, synapse=0.01)

    with nengo.Simulator(train_net) as sim_train:
        sim_train.run(period*3)

    solver = nengo.solvers.LstsqL2(reg=0.02)
    start = int(2*period/dt)
    end = int(3*period/dt)
    decoders, info = solver(sim_train.data[p_neur][start:end], target_func(np.arange(0, (end-start)*dt, dt))[:, None])

    return decoders


