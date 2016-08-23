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
t_len = t_period*t_scale
n_neurons = 600  # same as Matlab
sig_dims = 1


def conv_t(t):
    return t % t_len


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

t_steps = int(t_len / dt)

sin_per = (2 * np.pi * 10) / t_len
cos_per = (2 * np.pi * 20) / t_len
sigs = [
    lambda t: np.sin(sin_per*t),
    lambda t: 0.5*np.cos(cos_per*t),
    #lambda t: funky_sig(t)
]

n_sigs = len(sigs)
apert = np.ones(n_sigs) * 10

rate_data = np.zeros((n_sigs, t_steps, n_neurons))
pat_data = np.zeros((n_sigs, t_steps))
init_x = np.zeros((n_sigs, n_neurons))

w_rec = gen_w_rec(n_neurons)

for i_s, sig in enumerate(sigs):
    # get the rate data
    with nengo.Network() as rate_acc:
        in_sig = nengo.Node(sig)
        sig_reserv = nengo.Ensemble(n_neurons, 1, neuron_type=nengo.LIFRate(), seed=SEED)

        nengo.Connection(sig_reserv.neurons, sig_reserv.neurons, transform=w_rec, synapse=0)
        nengo.Connection(in_sig, sig_reserv, synapse=None)
        p_rate = nengo.Probe(sig_reserv.neurons, synapse=None)
        p_pat = nengo.Probe(in_sig)

    with nengo.Simulator(rate_acc) as rate_sim:
        rate_sim.run(t_len*2)

    rate_data[i_s] = rate_sim.data[p_rate][t_steps:]
    init_x[i_s] = rate_sim.data[p_rate][t_steps]
    pat_data[i_s] = rate_sim.data[p_pat][t_steps:, sig_dims-1]

# get the output weights using the usual decoder math
# previously verified by decoding original activities
solver = nengo.solvers.LstsqL2(reg=0.02)
w_out, d_info = solver(
    rate_data.reshape((t_steps*n_sigs, -1)),
    pat_data.reshape((t_steps*n_sigs, -1)),

)

# do SVD on the neuron data to get Conceptors
# slowest part of the process and appears to only be using one core
conceptors = get_conceptors(rate_data, n_sigs, t_steps, apert, n_neurons)


def get_idx(t):
    return int(np.floor(t % (t_len * n_sigs) / t_len))


def conc_func(t, x):
    """change periodically between output patterns for testing"""
    return np.dot(conceptors[get_idx(t)], x)


def kick_func(t):
    """initialise the neurons"""
    tt = conv_t(t)
    if tt < 0.1 * t_scale:
        return init_x[get_idx(t)]
    else:
        return 0

# use it for Conceptor testing
with nengo.Network() as conc_model:
    kick = nengo.Node(kick_func)
    reserv = nengo.Ensemble(n_neurons, 1, neuron_type=nengo.LIFRate(), seed=SEED)
    conc_node = nengo.Node(conc_func, size_in=n_neurons)
    output = nengo.Node(size_in=1)
    sanity = nengo.Node(get_idx)

    nengo.Connection(kick, reserv.neurons, synapse=None)
    nengo.Connection(reserv.neurons, conc_node, transform=w_rec, synapse=0)
    nengo.Connection(conc_node, reserv.neurons, synapse=None)
    nengo.Connection(reserv.neurons, output, transform=w_out.T)

    p_rate = nengo.Probe(reserv.neurons)
    p_res = nengo.Probe(output)
    p_kick = nengo.Probe(kick)
    p_conc = nengo.Probe(conc_node)
    p_san = nengo.Probe(sanity)

with nengo.Simulator(conc_model) as conc_sim:
    conc_sim.run(t_len*n_sigs)

plt.plot(conc_sim.trange(), conc_sim.data[p_res])
plt.show()

ipdb.set_trace()
