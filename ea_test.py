# test architecture with reduced functions

import scipy.io
from scipy import interpolate
import nengo
import numpy as np

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

output_dims = 2
pattern_num = 2
pattern_file_names = pattern_file_names[:2]

function_list = [
    [
        lambda x: x/2,
        lambda x: -x/2
    ],
    [
        lambda x: np.cos(np.pi*x)/2,
        lambda x: np.sin(np.pi*x)/2
    ]
]

np.random.seed(3)
# maps from input value (in this case, theta) to output value

model = nengo.Network()
tau = 0.1
model.config[nengo.Ensemble].neuron_type = nengo.LIFRate()
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

    # controllers
    inhibit_control = nengo.Node([0]*pattern_num)
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
        """
        for dim in range(output_dims):
            e.add_output("out_"+name+str(dim), function_list[n_i][dim])
        """
        e.add_output("out_"+name, function_list[n_i])
        e.add_neuron_input()

        # make the connections for inhibition (proxy for output of BG)
        nengo.Connection(inhibit_control[n_i], e.neuron_input,
                         transform=np.ones((100*output_dims, 1)) * -3)

        # TODO: make the connections for scaling (proxy for output of Thal)

        # make input and output connections
        nengo.Connection(readout, e.input, transform=np.ones((output_dims, 1)))
        """
        for dim in range(output_dims):
            nengo.Connection(getattr(e, "out_"+name+str(dim)), output.input[dim])
        """
        nengo.Connection(getattr(e, "out_"+name), output.input)

        ea_list.append(e)