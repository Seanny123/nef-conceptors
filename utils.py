import numpy as np


def d3_scale(dat, out_range=(-1, 1), in_range=None):
    if in_range is None:
        domain = [np.min(dat, axis=0), np.max(dat, axis=0)]
    else:
        domain = in_range

    def interp(x):
        return out_range[0] * (1.0 - x) + out_range[1] * x

    def uninterp(x):
        b = 0
        if (domain[1] - domain[0]) != 0:
            b = domain[1] - domain[0]
        else:
            b =  1.0 / domain[1]
        return (x - domain[0]) / b

    return interp(uninterp(dat))


def gen_w_rec(n_neurons):
    w_rec = np.random.uniform(-0.5, 0.5, size=(n_neurons, n_neurons)) / n_neurons
    w_rec /= np.max(np.abs(np.linalg.eigvals(w_rec)))
    w_rec /= n_neurons
    return w_rec


def get_conceptors(rate_data, n_sigs, t_steps, apert, n_neurons):
    conceptors = []
    for i_s in range(n_sigs):
        r_dat = rate_data[i_s]
        rate_corr = np.dot(r_dat.T, r_dat) / t_steps
        unit, sing_vals, _ = np.linalg.svd(rate_corr)
        s_new = np.dot(sing_vals, np.linalg.pinv(sing_vals + apert[i_s] ** -2 * np.eye(n_neurons)))
        conceptors.append(np.dot(np.dot(unit, s_new), unit.T))
    return conceptors
