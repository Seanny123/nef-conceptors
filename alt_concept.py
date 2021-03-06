from utils import gen_w_rec, get_conceptors_w_solver, check_w_out
from tanh_neuron import TanhWithBias

import nengo
import nengo.solvers
import nengo.utils.numpy as npext
import nengolib

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
wash = 50


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

old_rate_data = np.zeros((n_sigs, t_steps-wash, n_neurons))
rate_data = np.zeros((n_sigs, t_steps-wash, n_neurons))
w_targ_data = np.zeros((n_sigs, t_steps-wash, n_neurons))
pat_data = np.zeros((n_sigs, t_steps-wash))

init_x = np.zeros((n_sigs, n_neurons))

# TODO: just initialise this randomly
w_rec = gen_w_rec(n_neurons)
enc_dist = nengo.dists.UniformHypersphere(surface=True)
encoders = enc_dist.sample(n_neurons, sig_dims)

neuron_type = nengo.LIFRate()

for i_s, sig in enumerate(sigs):
    # get the rate data
    with nengo.Network() as rate_acc:
        in_sig = nengo.Node(sig)
        w_targ_out = nengo.Node(size_in=n_neurons)
        sig_reserv = nengo.Ensemble(n_neurons, sig_dims,
                                    neuron_type=neuron_type,
                                    encoders=encoders, seed=SEED)

        nengo.Connection(in_sig, sig_reserv, synapse=0)
        nengo.Connection(sig_reserv.neurons, sig_reserv.neurons,
                         transform=w_rec, synapse=0)
        nengo.Connection(sig_reserv, sig_reserv,
                         transform=0.6, synapse=0)
        nengo.Connection(sig_reserv.neurons, w_targ_out, transform=w_rec, synapse=None)

        p_w_targ = nengo.Probe(w_targ_out, synapse=None)
        p_rate = nengo.Probe(sig_reserv.neurons, synapse=None)
        p_pat = nengo.Probe(in_sig)

    with nengo.Simulator(rate_acc) as rate_sim:
        rate_sim.run(t_len*2)

    rate_data[i_s] = rate_sim.data[p_rate][t_steps+wash:]
    old_rate_data[i_s] = rate_sim.data[p_rate][t_steps+wash-1:-1]
    w_targ_data[i_s] = rate_sim.data[p_w_targ][t_steps+wash:]
    pat_data[i_s] = rate_sim.data[p_pat][t_steps+wash:, sig_dims-1]

    init_x[i_s] = rate_sim.data[p_rate][t_steps+wash]

old_x_val = old_rate_data.reshape(((t_steps-wash)*n_sigs, -1)).T
x_val = rate_data.reshape(((t_steps-wash)*n_sigs, -1)).T
sig_val = pat_data.reshape(((t_steps-wash)*n_sigs, -1)).T
w_targ = w_targ_data.reshape(((t_steps-wash)*n_sigs, -1)).T

solver = nengo.solvers.LstsqL2(reg=0.02)
w_out = solver(x_val.T, sig_val.T)[0]

# demonstrates how w_out is functional
#plt.plot(np.dot(rate_data.reshape(((t_steps-wash)*n_sigs, -1)), w_out))
#plt.show()

opt_w_rec = solver(old_x_val.T, w_targ.T)[0]
print(npext.rmse(np.dot(opt_w_rec, old_x_val), w_targ))

# do SVD on the neuron data to get Conceptors
# slowest part of the process and appears to only be using one core
#check_w_out(x_val, w_out, sig_val)


def get_idx(t):
    return int(np.floor(t % (t_len * n_sigs) / t_len))


def conc_func(t, x):
    """change periodically between output patterns for testing"""
    return np.dot(conceptors[get_idx(t)], x)


# this doesn't actually matter
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
    reserv = nengo.Ensemble(n_neurons, sig_dims, neuron_type=neuron_type,
                            encoders=encoders, seed=SEED)
    conc_node = nengo.Node(conc_func, size_in=n_neurons)
    output = nengo.Node(size_in=1)

    nengo.Connection(kick, reserv.neurons, synapse=None)
    nengo.Connection(reserv.neurons, conc_node,
                     transform=opt_w_rec, synapse=None)
    nengo.Connection(conc_node, reserv, synapse=0)
    nengo.Connection(reserv.neurons, output, transform=w_out.T, synapse=None)

    p_rate = nengo.Probe(reserv.neurons)
    p_res = nengo.Probe(output, synapse=0.01)
    p_kick = nengo.Probe(kick)
    p_conc = nengo.Probe(conc_node)

with nengo.Simulator(conc_model) as conc_sim:
    conc_sim.run(t_len*n_sigs)

plt.plot(conc_sim.trange(), conc_sim.data[p_res])
plt.show()

ipdb.set_trace()
