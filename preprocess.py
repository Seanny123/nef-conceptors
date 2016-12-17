import numpy as np
import scipy.io
from scipy import interpolate

from utils import d3_scale

pattern_file_names = (
    "nnRawExaStride",
    "nnRawSlowWalk",
    "nnRawWalk",
    "nnRawRunJog",
    "nnRawCartWheel",
    "nnRawWaltz",
    "nnRawCrawl",
    "nnRawStandup",
    "nnRawGetdown",
    "nnRawSitting",
    "nnRawGetSeated",
    "nnRawStandupFromStool",
    "nnRawBox1",
    "nnRawBox2",
    "nnRawBox3",
)


def preprocess(output_dims, pattern_num, save_res=False):
    """load the patterns from matlab and save the results"""

    ptrn_file_names = pattern_file_names[:pattern_num]
    min_maxs = np.zeros((output_dims, 2))
    min_maxs[:, 0] = np.inf
    min_maxs[:, 1] = -np.inf
    raw_dats = []

    final_normed = {}

    # get the actual maximum and minimums for each dimension
    for nm in ptrn_file_names:
        name = nm[5:]
        raw_dats.append(scipy.io.loadmat("section2.3_demoMotionCapture/nnData/%s.mat" % (nm))["nnRawData" + name].T)
        for o_i in range(output_dims):
            assert raw_dats[-1][o_i].shape != (61,)
            min_val = np.min(raw_dats[-1][o_i])
            if min_val < min_maxs[o_i, 0]:
                min_maxs[o_i, 0] = min_val

            max_val = np.max(raw_dats[-1][o_i])
            if max_val > min_maxs[o_i, 1]:
                min_maxs[o_i, 1] = max_val

    for n_i, nm in enumerate(ptrn_file_names):
        # make each pattern values normalised between -1, 1
        # and temporally squash them between -pi and pi too
        raw_dat = raw_dats[n_i]
        # assert raw_dat.shape[0] == output_dims
        normed_data = np.zeros_like(raw_dat)

        for o_i in range(output_dims):
            assert min_maxs[o_i][0] <= np.min(raw_dat[o_i, :])
            assert min_maxs[o_i][1] >= np.max(raw_dat[o_i, :])
            normed_data[o_i, :] = d3_scale(raw_dat[o_i, :], in_range=min_maxs[o_i])
            assert np.max(normed_data) <= 1.5
            assert np.min(normed_data) >= -1.5
        final_normed[nm] = normed_data

    assert np.all(min_maxs[:, 0] != np.inf)
    assert np.all(min_maxs[:, 1] != -np.inf)

    # save the result (min_maxs, normed_data)
    if save_res:
        np.savez("processed/out_%s_pn_%s.npz" % (output_dims, pattern_num),
                 min_maxs=min_maxs, **final_normed)

    return min_maxs, final_normed


def get_function_list(output_dims, pattern_num, filename="", trange=(-np.pi, np.pi)):
    """get the interpolated function"""

    function_list = []
    ptrn_file_names = pattern_file_names[:pattern_num]
    file_res = np.load(filename)

    for n_i, nm in enumerate(ptrn_file_names):

        normed_dat = file_res[nm]
        function_list.append([])
        xv = np.linspace(trange[0], trange[1], normed_dat.shape[1])

        for o_i in range(output_dims):
            function_list[-1].append(interpolate.interp1d(xv, normed_dat[o_i, :]))

    return function_list
