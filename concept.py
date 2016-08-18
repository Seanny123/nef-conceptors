from utils import gen_w_rec, get_conceptors

import nengo
import nengo.solvers

import numpy as np
import matplotlib.pyplot as plt

import ipdb

SEED = 0

dt = 0.001
t_scale = 0.5
t_period = 1.0
n_neurons = 1000
sig_dims = 1


def conv_t(t):
    return t % (t_period*t_scale)


def funky_sig(t):
    """rise with three different linear slopes"""
    t_step = conv_t(t)

    # getting a weird signal could have been done way more efficiently...
    if t_step < 0.3 * t_scale:
        return t_step
    elif 0.3*t_scale <= t_step < 0.6*t_scale:
        return t_step * 0.3
    elif 0.6 * t_scale <= t_step:
        return t_step * 1.3
    else:
        return 0

# TODO: make sure these are actually periodic
sigs = [
    lambda t: np.sin(20*t),
    lambda t: 0.5*np.cos(80*t),
    lambda t: funky_sig(t)
]

n_sigs = len(sigs)
apert = np.ones(n_sigs) * 10
t_len = t_period*t_scale
t_steps = int(t_len / dt)

rate_data = np.zeros((n_sigs, t_steps, n_neurons))
pat_data = np.zeros((n_sigs, t_steps))

w_rec = gen_w_rec(n_neurons)

for i_s, sig in enumerate(sigs):
    # get the rate data
    with nengo.Network() as rate_acc:
        in_sig = nengo.Node(sig)
        sig_reserv = nengo.Ensemble(n_neurons, 1, neuron_type=nengo.LIFRate(), seed=SEED)

        nengo.Connection(sig_reserv.neurons, sig_reserv.neurons, transform=w_rec)
        nengo.Connection(in_sig, sig_reserv, synapse=None)
        p_rate = nengo.Probe(sig_reserv.neurons)
        p_pat = nengo.Probe(in_sig)

    with nengo.Simulator(rate_acc) as rate_sim:
        rate_sim.run(t_len*2)

    rate_data[i_s] = rate_sim.data[p_rate][t_steps:]
    pat_data[i_s] = rate_sim.data[p_pat][t_steps:, sig_dims-1]

# get the output weights using the usual decoder math
solver = nengo.solvers.LstsqL2(reg=0.02)
w_out, _ = solver(
    rate_data.reshape((t_steps*n_sigs, -1)),
    pat_data.reshape((t_steps*n_sigs, -1)),

)
# TODO: compare RMSE

# do SVD on the neuron data to get Conceptors
# slowest part of the process and appears to only be using one core
conceptors = get_conceptors(rate_data, n_sigs, t_steps, apert, n_neurons)

# this process gives a n_neuron dim conceptor. Is this what Matlab gives?


def conc_func(t, x):
    """change periodically between output patterns for testing"""
    # is this supposed to be a dot product or a multiplication?
    return np.dot(x, conceptors[int(t % t_len)])


# use it for Conceptor testing
with nengo.Network() as conc_model:
    reserv = nengo.Ensemble(n_neurons, 1, neuron_type=nengo.LIFRate(), seed=SEED)
    conc_node = nengo.Node(conc_func, size_in=n_neurons)
    output = nengo.Node(size_in=1)

    nengo.Connection(reserv.neurons, conc_node, transform=w_rec)
    nengo.Connection(reserv.neurons, output, transform=w_out[None, :, sig_dims-1])

    p_res = nengo.Probe(output)

with nengo.Simulator(conc_model) as conc_sim:
    conc_sim.run(t_len*n_sigs)

plt.plot(conc_sim.data[p_res])
plt.show()
ipdb.set_trace()
