import scipy.io
from scipy import interpolate
import nengo
import numpy as np
import ipdb

def d3_scale(dat, out_range=(-1, 1), in_range=None):
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

output_dims = 61
pattern_num = 2
pattern_file_names = pattern_file_names[:2]

min_maxs = np.zeros((output_dims, 2))
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

function_list = []

for n_i, nm in enumerate(pattern_file_names):
    # make each pattern values normalised between -1, 1
    # and temporally squash them between -pi and pi too
    function_list.append([])
    raw_dat = raw_dats[n_i]
    xv = np.linspace(-np.pi, np.pi, raw_dat.shape[1])
    assert raw_dat.shape[0] == output_dims
    normed_data = np.zeros_like(raw_dat)

    for o_i in range(output_dims):
        assert min_maxs[o_i][0] <= np.min(raw_dat[o_i, :]) 
        assert min_maxs[o_i][1] >= np.max(raw_dat[o_i, :]) 
        normed_data[o_i, :] = d3_scale(raw_dat[o_i, :], in_range=min_maxs[o_i])
        assert np.max(normed_data) <= 1.5
        assert np.min(normed_data) >= -1.5
        function_list[-1].append(interpolate.interp1d(xv, normed_data[o_i, :]))


np.random.seed(3)
# maps from input value (in this case, theta) to output value

model = nengo.Network()
tau = 0.1
#model.config[nengo.Ensemble].neuron_type = nengo.Direct()
ea_list = []
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

    # what's stopping this from being a passthrough node?
    readout = nengo.Ensemble(n_neurons=1, dimensions=1, neuron_type=nengo.Direct())
    nengo.Connection(osc[:2], readout, synapse=None,
                     function=lambda x: np.arctan2(x[1], x[0]))

    # controllers # TODO: Make smoother transition
    inhibit_control = nengo.Node(lambda t: [0,1] if t < 2 else [1,0])
    #scale_control = nengo.Node([0]*pattern_num)

    output = nengo.networks.EnsembleArray(n_neurons=100, n_ensembles=output_dims,
                                          label="output")

    # one ensemble array per output pattern
    # each ensemble array has the output dimensions
    # combine ensemble arrays to combine patterns
    # BONUS: figure out how to combine individual dimensions, once I know what they actually mean...
    # SUPER BONUS: use visual assement for the robot to be able to imitate a movement
    # ASIDE: which would be cool, because then maybe the robot could infer properties of objects from movement

    for n_i, nm in enumerate(pattern_file_names):
        # make the EnsembleArray with the associated functions
        name = nm[5:]
        e = nengo.networks.EnsembleArray(n_neurons=100, n_ensembles=output_dims,
                                label=name)
        e.add_output("out_"+name, function_list[n_i])
        e.add_neuron_input()

        # make the connections for inhibition (proxy for output of BG)
        nengo.Connection(inhibit_control[n_i], e.neuron_input,
                         transform=np.ones((100*output_dims, 1)) * -3)

        # TODO: make the connections for scaling (proxy for output of Thal)

        # make input and output connections
        nengo.Connection(readout, e.input, transform=np.ones((output_dims, 1)))
        nengo.Connection(getattr(e, "out_"+name), output.input)

        ea_list.append(e)

    # probe the output
    p_out = nengo.Probe(output.output)

with nengo.Simulator(model) as sim:
    sim.run(4)

# un-normalise on export based off the original domain
tmp = sim.data[p_out]
reg_out = np.zeros_like(tmp)
for t_i in range(tmp.shape[0]):
    reg_out[t_i, :] = d3_scale(tmp[t_i, :], out_range=min_maxs[o_i],
                               in_range=(-1, 1))
ipdb.set_trace()

# try running the patterns in Matlab to see if they're legit
scipy.io.savemat("pattern_out.mat", reg_out)