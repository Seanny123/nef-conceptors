import scipy.io
from scipy import interpolate
from sklearn.preprocessing import scale
import nengo
import numpy as np

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
pattern_file_names = [
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
]

output_dims = 61
pattern_num = 2
pattern_file_names = pattern_file_names[:2]

function_list = []
min_maxs = np.zeros((output_dims, 2))
raw_dats = []

# get the actual maximum and minimums for each dimension
for nm in pattern_file_names:
    raw_dats.append(scipy.io.loadmat("section2.3_demoMotionCapture/nnData/%s.mat" %(nm))[nm])
    for o_i in range(output_dims):
        min_val = np.min(raw_dats[-1][o_i])
        if min_val < min_maxs[o_i, 0]:
            min_maxs[o_i, 0] = min_val

        max_val = np.max(raw_dats[-1][o_i])
        if max_val > min_maxs[o_i, 1]:
            min_maxs[o_i, 1] = max_val

for n_i, nm in enumerate(pattern_file_names):
    # make each pattern values normalised between -1, 1
    # and temporally squash them between pi/2 and -pi/2
    function_list.append([])
    raw_dat = raw_dats[n_i]
    xv = np.linspace(-np.pi/2, np/2, raw_dat.shape[0])
    assert raw_dat.shape[1] == output_dims
    normed_data = np.zeros_like(raw_dat)

    for o_i in range(output_dims):
        normed_data[:, o_i] = d3_scale(raw_dat[:, o_i], in_range=min_maxs[o_i])
        assert np.max(normed_data) == 1
        assert np.min(normed_data) == -1
        function_list[-1].append(interpolate.interp1d(xv, normed_dat[:, o_i]))


np.random.seed(3)
# maps from input value (in this case, theta) to output value

model = nengo.Network()
tau = 0.1
model.config[nengo.Ensemble].neuron_type = nengo.Direct()
ea_list = []
with model:
    osc = nengo.Ensemble(n_neurons=1500, dimensions=3, radius=1.4)

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
    readout = nengo.Ensemble(n_neurons=500, dimensions=1)
    nengo.Connection(osc[:2], readout,
                     function=lambda x: np.arctan2(x[1], x[0]))

    # controllers
    inhibit_control = nengo.Node([0]*pattern_num)
    #scale_control = nengo.Node([0]*pattern_num)

    output = nengo.EnsembleArray(n_neurons=100, dimensions=output_dims,
                                 label="output")

    # one ensemble array per output pattern
    # each ensemble array has the output dimensions
    # combine ensemble arrays to combine patterns
    # BONUS: figure out how to combine individual dimensions, once I know what they actually mean...
    # SUPER BONUS: use visual assement for the robot to be able to imitate a movement
    # ASIDE: which would be cool, because then maybe the robot could infer properties of objects from movement

    for n_i, nm in enumerate(pattern_file_names):
        # make the EnsembleArray with the associated functions
        name = nm[10:]
        e = nengo.EnsembleArray(n_neurons=100, n_ensembles=output_dims,
                                label=nm[10:])
        for dim in range(output_dims):
            e.add_ouput("out_"+nm[10:], function_list[n_i][dim])
        e.add_neuron_input()

        # make the connections for inhibition (proxy for output of BG)
        nengo.Connection(inhibit_control[n_i], e.neuron_input)

        # TODO: make the connections for scaling (proxy for output of Thal)

        # make input and output connections
        nengo.Connection(readout, e.input)
        nengo.Connection(e.output, output)

        ea_list.append(e)

    # probe the output
    p_out = nengo.Probe(output.output)

with nengo.Simulator(model) as sim:
    sim.run(1)

# un-normalise on export based off the original domain
tmp = sim.data["p_out"]
reg_out.zeros_like(tmp)
for t_i in range(tmp.shape[0]):
    reg_out[t_i, :] = d3_scale(tmp[t_i, :], out_range=min_maxs[o_i])


# try running the patterns in Matlab to see if they're legit
scipy.io.savemat("pattern_out.mat", reg_out)