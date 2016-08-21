from utils import gen_w_rec, get_conceptors

import nengo
import nengo.solvers

import numpy as np

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

rate_data = np.zeros((n_sigs, t_len))
pat_data = np.zeros((n_sigs, t_len))
w_rec_base = gen_w_rec(n_neurons)


def make_div_net(origin, w_func):
    sig_ens = []
    with nengo.Network as d_n:
        for s_i in range(n_sigs):
            sig_ens.append(nengo.Ensemble(n_neurons, 1, seed=SEED))
            nengo.Connection(origin.neurons, sig_ens[-1].neurons, transform=w_rec_base)
            nengo.Connection(sig_ens[-1].neurons, origin.neurons, transform=w_func(n_neurons))
    return d_n, sig_ens

for i_s, sig in enumerate(sigs):
    # get the rate data
    with nengo.Network() as rate_acc:
        in_sig = nengo.Node(sig)
        sig_reserv = nengo.Ensemble(n_neurons, 1, neuron_type=nengo.LIFRate(), seed=SEED)

        div_net, sig_div = make_div_net(sig_reserv)

        nengo.Connection(in_sig, sig_reserv, synapse=None)
        p_rate = nengo.Probe(sig_reserv.neurons)
        p_pat = nengo.Probe(in_sig)

    with nengo.Simulator(rate_acc) as rate_sim:
        rate_sim.run(t_period*t_scale*2)

    rate_data[i_s, :] = (rate_sim.data[p_rate][t_steps:])
    pat_data[i_s, :] = (rate_sim.data[p_pat][t_steps:])

# get the output weights using the usual decoder math
solver = nengo.solvers.LstsqL2(reg=0.02)
w_out, _ = solver(rate_data, pat_data)

# do SVD on the neuron data to get Conceptors
conceptors = get_conceptors(rate_data, n_sigs, t_steps, apert, n_neurons)


# use it for Conceptor testing
with nengo.Network() as conc_model:
    reserv = nengo.Ensemble(n_neurons, 1, neuron_type=nengo.LIFRate(), seed=SEED)

    class ConcGen:
        def __init__(self):
            self.count = 0

        def w_conc_gen(self, n):
            w = np.dot(conceptors[self.count], gen_w_rec(n))
            self.count += 1
            return w

    conc_gen = ConcGen()
    conc_net, conc_div = make_div_net(reserv, conc_gen.w_conc_gen)

    output = nengo.Node(size_in=1)

    nengo.Connection(reserv.neurons, output, transform=w_out[None, :, sig_dims-1])

    p_res = nengo.Probe(output)
