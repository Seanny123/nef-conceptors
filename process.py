import nengo
import numpy as np
import scipy
import scipy.io
from scipy import interpolate

# TODO: remove "d3_scale" and "pre", then rename the file to postprocess


def d3_scale(dat, out_range=(-1, 1), in_range=None):
    """scale function for mapping from one range to another stolen from D3.js"""

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


def pre(output_dims, pattern_file_names):
    min_maxs = np.zeros((output_dims, 2))
    min_maxs[:, 0] = np.inf
    min_maxs[:, 1] = -np.inf
    raw_dats = []

    # get the actual maximum and minimums for each dimension
    for nm in pattern_file_names:
        name = nm[5:]
        raw_dats.append(scipy.io.loadmat("section2.3_demoMotionCapture/nnData/%s.mat" %(nm))["nnRawData"+name].T)
        for o_i in range(output_dims):
            assert raw_dats[-1][o_i].shape != (61,)
            min_val = np.min(raw_dats[-1][o_i])
            if min_val < min_maxs[o_i, 0]:
                min_maxs[o_i, 0] = min_val

            max_val = np.max(raw_dats[-1][o_i])
            if max_val > min_maxs[o_i, 1]:
                min_maxs[o_i, 1] = max_val

    assert np.all(min_maxs[:, 0] != np.inf)
    assert np.all(min_maxs[:, 1] != -np.inf)

    function_list = []

    for n_i, nm in enumerate(pattern_file_names):
        # make each pattern values normalised between -1, 1
        # and temporally squash them between -pi and pi too
        function_list.append([])
        raw_dat = raw_dats[n_i]
        xv = np.linspace(-np.pi, np.pi, raw_dat.shape[1])
        #assert raw_dat.shape[0] == output_dims
        normed_data = np.zeros_like(raw_dat)

        for o_i in range(output_dims):
            assert min_maxs[o_i][0] <= np.min(raw_dat[o_i, :]) 
            assert min_maxs[o_i][1] >= np.max(raw_dat[o_i, :]) 
            normed_data[o_i, :] = d3_scale(raw_dat[o_i, :], in_range=min_maxs[o_i])
            assert np.max(normed_data) <= 1.5
            assert np.min(normed_data) >= -1.5
            function_list[-1].append(interpolate.interp1d(xv, normed_data[o_i, :]))
    return function_list, min_maxs


def post(sim, p_out, min_maxs, final_only=True):
    tmp_out = sim.data[p_out]
    reg_out = np.zeros_like(tmp_out)
    reg_min = np.min(tmp_out)
    reg_max = np.max(tmp_out)
    out_scaled = np.zeros_like(tmp_out)
    compressed = []
    final_out = []

    for t_i in range(tmp_out.shape[1]):
        # additional filtering for reg_out to stop the crazy jitters
        tmp_scaled = d3_scale(tmp_out[:, t_i], out_range=min_maxs[t_i],
                              in_range=(reg_min, reg_max))
        print("Initial min: %s Initial max: %s" %(reg_min, reg_max))
        print("Scaled min: %s Scaled max: %s" %(np.min(tmp_scaled), np.max(tmp_scaled)))
        print("Intended min: %s Intended max %s" %(min_maxs[t_i][0], min_maxs[t_i][1]))

        good_start = 1000
        total_len = tmp_scaled.shape[0] - good_start
        # this doesn't do a perfect translation, because some actions are longer/shorter than others
        # but that doesn't really matter in the end
        signal_len = 291
        sample_interval = int( total_len / signal_len)
        re_size = good_start + total_len % sample_interval
        compressed.append(nengo.Lowpass(0.005).filt(
                          tmp_scaled[re_size:][::sample_interval], dt=0.001))
        print("Compressed min: %s, Compressed max: %s" %(np.min(compressed[t_i]), np.max(compressed[t_i])))
        final_out.append(d3_scale(compressed[t_i],
                                  in_range=(np.min(compressed[t_i]), np.max(compressed[t_i])),
                                  out_range=min_maxs[t_i]))
        final_max = np.max(final_out[t_i])
        final_min = np.min(final_out[t_i])
        print("Final min: %s, Final max: %s\n" %(final_min, final_max))
        print("Final diff: %s\n\n" %(
              np.abs(min_maxs[t_i][0]-final_min) + np.abs(min_maxs[t_i][1]-final_max)))

        if not final_only:
            out_scaled[:, t_i] = tmp_scaled
            reg_out[:, t_i] = nengo.Lowpass(0.1).filt(tmp_scaled, dt=0.001)

    # try running the patterns in Matlab to see if they're legit
    if not final_only:
        scipy.io.savemat("pattern_out.mat", {"out_scaled": out_scaled, "reg_out": reg_out, "compressed": compressed})
    scipy.io.savemat("final_pattern.mat", {"final_out": final_out})
