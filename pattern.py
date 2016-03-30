# complete test

from dmp_utils import *

import scipy.io
from scipy import interpolate
import nengo
import numpy as np
import ipdb

def d3_scale(dat, out_range=(-1, 1), in_range=None):
    """scale function for mapping from one range to another stolen from D3.js"""

    if in_range == None:
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

def make_dmp_net(functions, input_obj, output_obj, name=""):
    """create one point attractor per dimension with goals as nodes
    and one unified neuron input for inhibition

    TODO:
        - make the connections for inhibition (proxy for output of BG)
        - make the connections for scaling (proxy for output of Thal)
    """
    n = nengo.Network(label=name)
    n.pt_attractors = []
    n.conn_funcs = []
    n.f_conns = []
    with n:
        n.output = nengo.Node(size_in=len(functions))

        for d in range(len(functions)):
            goal = nengo.Node([0], label="goal_%s" %(d))
            attractor = gen_point_attractor(n, goal, n_neurons=500)
            attractor.label = "pt_attr_%s" %(d)
            nengo.Connection(attractor[0], n.output[d], synapse=None)
            n.pt_attractors.append(attractor)

    for f_i, func in enumerate(functions):
        dest = func(np.linspace(-np.pi, np.pi, ea_func_steps))
        dest = dest.reshape((-1, 1))
        force_func = gen_forcing_functions(dest, num_samples=50)[0]
        n.conn_funcs.append(lambda x, force_func=force_func: force(x, force_func))

        n.f_conns.append(nengo.Connection(input_obj, n.pt_attractors[f_i][1], synapse=None,
                         function=n.conn_funcs[f_i]))

    nengo.Connection(n.output, output_obj)
    return n

# load the patterns from matlab
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

# max is 61, but 14 is a nice leg
output_dims = 5
pattern_num = 1
pattern_file_names = pattern_file_names[:pattern_num]

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

ea_n_neurons = 300
ea_func_steps = 100
np.random.seed(3)
# maps from input value (in this case, theta) to output value

model = nengo.Network()
tau = 0.1
#model.config[nengo.Ensemble].neuron_type = nengo.Direct()
dmp_net_list = []

with model:
    osc = nengo.Ensemble(n_neurons=1, dimensions=3, neuron_type=nengo.Direct())

    def cycle(x):
        """makes a speed controlled oscillator"""
        a = 1.0
        b = 2.0 * np.pi * x[2]
        r = np.sqrt(x[0]**2.0 + x[1]**2.0)
        theta = np.arctan2(x[1], x[0])
        dr = 10.0*(-r**3.0 + a*r)
        dtheta = b

        dx = dr*np.cos(theta) - r*np.sin(theta)*dtheta
        dy = dr*np.sin(theta) + r*np.cos(theta)*dtheta

        return [x[0] + tau*dx, x[1] + tau*dy]

    nengo.Connection(osc, osc[:2], synapse=tau, function=cycle)

    rate = nengo.Node([1])
    nengo.Connection(rate, osc[2])

    bump = nengo.Node(lambda t: 1 if t < 0.5 else 0)
    nengo.Connection(bump, osc[0])

    # controllers # TODO: Make smoother transition
    #inhibit_control = nengo.Node(lambda t: [0,1] if t < 2 else [1,0])
    #inhibit_control = nengo.Node([0])
    #scale_control = nengo.Node([0]*pattern_num)

    output = nengo.networks.EnsembleArray(n_neurons=1, n_ensembles=output_dims, radius=np.pi,
                                          neuron_type=nengo.Direct(), label="output")

    # one ensemble array per output pattern
    # each ensemble array has the output dimensions
    # combine ensemble arrays to combine patterns
    # BONUS: figure out how to combine individual dimensions, once I know what they actually mean...
    # SUPER BONUS: use visual assement for the robot to be able to imitate a movement
    # ASIDE: which would be cool, because then maybe the robot could infer properties of objects from movement

    ea_n_neurons = 300

    for n_i, nm in enumerate(pattern_file_names):
        name = nm[5:]
        print(name)
        n = make_dmp_net(function_list[n_i], osc[:2], output.input, name=name)
        dmp_net_list.append(n)

    # Helper nodes for verification
    arctan = nengo.Node(size_in=1)
    nengo.Connection(osc[:2], arctan, function=lambda x: np.arctan2(x[0], x[1]), synapse=None)
    ideal = nengo.Node(lambda t, x: function_list[0][0](x), size_in=1)
    nengo.Connection(arctan, ideal)


    # probe the output
    p_out = nengo.Probe(output.output, synapse=0.15)
    p_ideal = nengo.Probe(ideal)

with nengo.Simulator(model) as sim:
    sim.run(2)

# un-normalise on export based off the original domain
# NOTE: ideal is going to break for dimensions larger than 1
tmp_out = sim.data[p_out]
tmp_ideal = sim.data[p_ideal]
reg_out = np.zeros_like(tmp_out)
reg_min = np.min(tmp_out)
reg_max = np.max(tmp_out)
out_scaled = np.zeros_like(tmp_out)
compressed = []
final_out = []
#ideal_out = []
for t_i in range(tmp_out.shape[1]):
    # additional filtering for reg_out to stop the crazy jitters
    # TODO: compare with guassian?
    tmp_scaled = d3_scale(tmp_out[:, t_i], out_range=min_maxs[t_i],
                          in_range=(reg_min, reg_max))
    print("Iinitial min: %s Initial max: %s" %(reg_min, reg_max))
    print("Scaled min: %s Scaled max: %s" %(np.min(tmp_scaled), np.max(tmp_scaled)))
    print("Intended min: %s Intended max %s" %(min_maxs[t_i][0], min_maxs[t_i][1]))

    out_scaled[:, t_i] = tmp_scaled

    good_shape = 1000
    sample_interval = int( good_shape / raw_dats[0].shape[1])
    re_size = good_shape + good_shape % sample_interval
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

    reg_out[:, t_i] = nengo.Lowpass(0.1).filt(tmp_scaled, dt=0.001)
    #ideal_out.append(d3_scale(tmp_ideal[:, t_i], out_range=min_maxs[0],
    #                       in_range=(-1, 1))[re_size:][::sample_interval])


# try running the patterns in Matlab to see if they're legit
# I GOT A LEG
#scipy.io.savemat("pattern_out.mat", {"out_scaled": out_scaled, "reg_out": reg_out, "ideal_out": ideal_out, "compressed": compressed})
scipy.io.savemat("pattern_out.mat", {"out_scaled": out_scaled, "reg_out": reg_out, "compressed": compressed})
scipy.io.savemat("final_pattern.mat", {"final_out": final_out})
