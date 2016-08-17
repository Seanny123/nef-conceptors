import nengo
import nengo.solvers

import numpy as np

SEED = 0

dt = 0.001
t_scale = 0.5
t_period = 1.0


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

rate_data = np.zeros((n_sigs, t_len))  # this should be a numpy array
pat_data = np.zeros((n_sigs, t_len))
n_neurons = 1000
w_rec = np.random.uniform(-0.5, 0.5, size=(n_neurons, n_neurons)) / n_neurons
w_rec /= np.max(np.abs(np.linalg.eigvals(w_rec)))
w_rec /= n_neurons

for i_s, sig in enumerate(sigs):
    # get the rate data
    with nengo.Network as rate_acc:
        in_sig = nengo.Node(sig)
        sig_reserv = nengo.Ensemble(n_neurons, 1, neuron_type=nengo.LIFRate(), seed=SEED)

        nengo.Connection(sig_reserv.neurons, sig_reserv.neurons, transform=w_rec)
        nengo.Connection(in_sig, sig_reserv, synapse=None)
        p_rate = nengo.Probe(sig_reserv.neurons)
        p_pat = nengo.Probe(in_sig)

    with nengo.Simulator(rate_acc) as rate_sim:
        rate_sim.run(t_period*t_scale*2)

    rate_data[i_s, :] = (rate_sim.data[p_rate][t_period*t_scale:])
    pat_data[i_s, :] = (rate_sim.data[p_pat][t_period*t_scale:])

# get the output weights using the usual decoder math
solver = nengo.solvers.LstsqL2(reg=0.02)
w_out, _ = solver(rate_data, pat_data)

# do SVD on the neuron data to get Conceptors
conceptors = []
for i_s in range(n_sigs):
    r_dat = rate_data[i_s, :]
    rate_corr = np.dot(r_dat, r_dat.T) / t_len
    unit, sing_vals, _ = np.linalg.svd(rate_corr)
    assert sing_vals.shape == np.eye(n_neurons)
    s_new = np.dot(sing_vals, np.linalg.pinv(sing_vals + apert[i_s]**-2 * np.eye(n_neurons)))
    conceptors.append(np.dot(np.dot(unit, s_new), unit.T))


def conc_func(t, x):
    """change periodically between output patterns for testing"""
    return x * conceptors[int(t % (t_period*t_scale))]


# use it for Conceptor testing
with nengo.Network as conc_model:
    reserv = nengo.Ensemble(n_neurons, 1, neuron_type=nengo.LIFRate(), seed=SEED)
    conc_node = nengo.Node(conc_func)
    output = nengo.Node(size_in=1)

    nengo.Connection(reserv.neurons, conc_node, transform=w_rec)
    nengo.Connection(reserv.neurons, output, transform=w_out)

    p_res = nengo.Probe(output)
