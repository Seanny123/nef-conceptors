import numpy as np
import scipy
import matplotlib.pyplot as plt
import nengo
import nengo.utils.numpy as npext

#import ipdb


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
            b = 1.0 / domain[1]
        return (x - domain[0]) / b

    return interp(uninterp(dat))


def gen_w_rec(n_neurons):
    #w_rec = scipy.sparse.random(n_neurons, n_neurons, density=10/n_neurons).A
    w_rec = np.random.randn(n_neurons, n_neurons)
    w_rec /= np.max(np.abs(np.linalg.eigvals(w_rec)))
    print("max weight:%s" % np.max(w_rec))
    return w_rec


def get_w(rate_data, target, n_neurons=1, reg=0):
    """The solution method from the Matlab implementation"""
    return np.dot(
             np.dot(
               np.linalg.pinv(
                  np.dot(rate_data, rate_data.T) + reg * np.eye(n_neurons)
               ), rate_data),
           target.T)


def get_conceptors(rate_data, n_sigs, t_steps, apert, n_neurons):
    """this is nearly identical to
     https://github.com/nengo/nengo/blob/master/nengo/utils/least_squares_solvers.py#L250"""
    conceptors = []
    for i_s in range(n_sigs):
        r_dat = rate_data[i_s]
        # TODO: replace t_steps with just te length of r_dat
        rate_corr = np.dot(r_dat.T, r_dat) / t_steps
        print("max correlation:%s" % np.max(rate_corr))
        unit, sing_vals, _ = np.linalg.svd(rate_corr)
        print("max sing vals:%s" % np.max(sing_vals))
        #plt.plot(sing_vals[:100])
        #plt.show()
        sing_vals = np.diag(sing_vals)
        # This is equivalent to the SVD solver
        # si = s / (s**2 + m * sigma**2)
        s_new = np.dot(
            sing_vals,
            np.linalg.pinv(
                sing_vals + (apert[i_s] ** -2) * np.eye(n_neurons)
            )
        )
        assert s_new.shape == sing_vals.shape
        # this is the one thing that's different
        # because unlike the SVD solver, there's no target
        # np.dot(V.T, si[:, None] * np.dot(U.T, Y))
        # you're just trying to reproduce the input
        # so if you give the reference signal it is the
        # same thing
        c_res = np.dot(np.dot(unit, s_new), unit.T)
        assert c_res.shape == s_new.shape
        assert c_res.shape == (n_neurons, n_neurons)
        conceptors.append(c_res)
    return conceptors

def get_conceptors_w_solver(rate_data, sig_val, apert):
    conceptors = []
    svd_solver = nengo.utils.least_squares_solvers.SVD()
    n_sigs = rate_data.shape[0]

    sigs = sig_val.reshape((n_sigs, -1))
    for s_i in range(n_sigs):
        conceptors.append(
            svd_solver(
                rate_data[s_i],
                sigs[s_i],
                1.0 / apert[s_i]
            )[0]
        )
    return conceptors

def check_w_out(x_val, w_out, sig_val):
    approx = np.dot(x_val.T, w_out)
    print(npext.rmse(approx, sig_val.T))
    
    plt.plot(approx, alpha=0.5)
    plt.plot(sig_val.T, alpha=0.5)
    plt.show()