from utils import gen_w_rec, get_w_out, nrmse

import nengo
import nengo.solvers

import numpy as np
import matplotlib.pyplot as plt
from tanh_neuron import TanhWithBias

import ipdb

SEED = 0

dt = 0.001
t_scale = 0.5
t_period = 1.0
t_len = t_period*t_scale
n_neurons = 600
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
]

n_sigs = len(sigs)
apert = np.ones(n_sigs) * 10

rate_data = np.zeros((n_sigs, t_steps, n_neurons))
pat_data = np.zeros((n_sigs, t_steps))
init_x = np.zeros((n_sigs, n_neurons))

w_rec = gen_w_rec(n_neurons)
#neuron_type = nengo.LIFRate()
neuron_type = TanhWithBias(seed=SEED)

for i_s, sig in enumerate(sigs):
    # get the rate data
    with nengo.Network(seed=SEED) as rate_acc:
        in_sig = nengo.Node(sig)
        sig_reserv = nengo.Ensemble(n_neurons, sig_dims, neuron_type=neuron_type, seed=SEED)

        nengo.Connection(sig_reserv.neurons, sig_reserv.neurons, transform=w_rec, synapse=0, seed=SEED)
        nengo.Connection(in_sig, sig_reserv, synapse=None, seed=SEED)
        p_rate = nengo.Probe(sig_reserv.neurons, synapse=None)
        p_pat = nengo.Probe(in_sig)

    with nengo.Simulator(rate_acc) as rate_sim:
        rate_sim.run(t_len*2)

    rate_data[i_s] = rate_sim.data[p_rate][t_steps:]
    init_x[i_s] = rate_sim.data[p_rate][t_steps]
    pat_data[i_s] = rate_sim.data[p_pat][t_steps:, sig_dims-1]

rate_data = rate_data.reshape((t_steps*n_sigs, -1))
pat_data = pat_data.reshape((t_steps*n_sigs, -1))

w_out = get_w_out(rate_data.T, pat_data.T, n_neurons=n_neurons)

print("NRMSE: %s" % np.mean(nrmse(rate_data.T, pat_data.T)))
plt.plot(np.dot(rate_data, w_out))
plt.show()

# use it for Conceptor testing
with nengo.Network(seed=SEED) as conc_model:
    in_sig = nengo.Node(sig)
    reserv = nengo.Ensemble(n_neurons, sig_dims, neuron_type=neuron_type, seed=SEED)
    conc_node = nengo.Node(size_in=n_neurons)
    output = nengo.Node(size_in=sig_dims)

    nengo.Connection(in_sig, reserv, synapse=None)
    nengo.Connection(reserv.neurons, reserv.neurons, transform=w_rec, synapse=0, seed=SEED)
    nengo.Connection(reserv.neurons, output, transform=w_out.T, synapse=None, seed=SEED)

    p_res = nengo.Probe(output, synapse=0.01)

with nengo.Simulator(conc_model) as conc_sim:
    conc_sim.run(t_len*n_sigs)

plt.plot(conc_sim.trange(), conc_sim.data[p_res])
plt.show()

ipdb.set_trace()
