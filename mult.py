# try boosting a pattern

from dmp_utils import *
from process import *

import scipy.io
from scipy import interpolate
import nengo
import numpy as np
import ipdb


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

def make_mult_net(n_neurons, dims, name=""):
    n = nengo.Network(label=name)
    with n:
        n.input_sig = nengo.Node(size_in=dims)
        n.input_scale = nengo.Node(size_in=1)
        n.output = nengo.Node(size_in=dims)
        # TODO: add neuron input later for inhibition

        mult_encs = nengo.dists.Choice([[1,1],[1,-1],[-1,-1],[-1,1]])
        for d in range(dims):
            mult_ens = nengo.Ensemble(n_neurons=n_neurons, dimensions=2,
                                   encoders=mult_encs, neuron_type=nengo.LIFRate(), radius=2)
            nengo.Connection(n.input_sig[d], mult_ens[0])
            nengo.Connection(n.input_scale, mult_ens[1])
            nengo.Connection(mult_ens, n.output[d], function=lambda x: x[0]*x[1])
    return n


# load the patterns from matlab
pattern_file_names = (
     "nnRawRunJog",
     "nnRawExaStride",
     "nnRawSlowWalk",
     "nnRawWalk",
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
output_dims = 15
pattern_num = 1
pattern_file_names = pattern_file_names[:pattern_num]

function_list, min_maxs = pre(output_dims, pattern_file_names)

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

    bump = nengo.Node(lambda t: 1 if t < 0.05 else 0)
    nengo.Connection(bump, osc[0])

    # controllers # TODO: Make smoother transition
    #inhibit_control = nengo.Node(lambda t: [0,1] if t < 2 else [1,0])
    #inhibit_control = nengo.Node([0])
    scale_control = nengo.Node([1]*pattern_num)

    output = nengo.networks.EnsembleArray(n_neurons=1, n_ensembles=output_dims, radius=np.pi,
                                          neuron_type=nengo.Direct(), label="output")

    # one ensemble array per output pattern
    # each ensemble array has the output dimensions
    # combine ensemble arrays to combine patterns

    ea_n_neurons = 300

    for n_i, nm in enumerate(pattern_file_names):
        name = nm[5:]
        print(name)
        # first just get inhibition working
        # then convert this normal ensemble array into a fancy multiplication array
        mult_ea = make_mult_net(500, output_dims, name="mult"+name)

        n = make_dmp_net(function_list[n_i], osc[:2], mult_ea.input_sig, name=name)

        nengo.Connection(scale_control[n_i], mult_ea.input_scale)
        nengo.Connection(mult_ea.output, output.input)
        dmp_net_list.append(n)


    # probe the output
    p_out = nengo.Probe(output.output, synapse=0.15)

with nengo.Simulator(model) as sim:
    sim.run(2)

post(sim, p_out, min_maxs)
