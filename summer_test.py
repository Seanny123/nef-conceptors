import nengo
import numpy as np
from scipy import interpolate

def gen_forcing_functions(y_des, dt=.001, alpha=10, beta=10/4.):
        
        # scale our trajectory and find the center point
        y_des = y_des.T / 1e5
        goal = np.sum(y_des, axis=1) / y_des.shape[1]
                
        # interpolate our desired trajectory to smooth out the sampling
        #s this makes it more boxy by drawing lines between the sample points
        num_samples = 10
        path = np.zeros((y_des.shape[0], num_samples))
        x = np.linspace(-np.pi, np.pi, y_des.shape[1])
        for d in range(y_des.shape[0]):
            path_gen = interpolate.interp1d(x, y_des[d])
            for ii,t in enumerate(np.linspace(-np.pi, np.pi, num_samples)):
                path[d, ii] = path_gen(t)
        y_des = path
    
        # calculate velocity of y_des
        #s calculate the velocity for every point of the path
        dy_des = np.diff(y_des) / dt
        # add zero to the beginning of every row
        dy_des = np.hstack((np.zeros((y_des.shape[0], 1)), dy_des))

        # calculate acceleration of y_des
        ddy_des = np.diff(dy_des) / dt
        # add zero to the beginning of every row
        ddy_des = np.hstack((np.zeros((y_des.shape[0], 1)), ddy_des))

        forcing_functions = []
        for d in range(y_des.shape[0]):
            # find the force required to move along this trajectory
            # by subtracting out the effects of the point attractor 
            force = ddy_des[d] - alpha * \
                            (beta * (goal[d] - y_des[d]) - \
                             dy_des[d])
            # generate another interpolation function we can use 
            # to now train up our decoded oscillator output
            forcing_functions.append(lambda x, force=force:
                                    interpolate.interp1d(np.linspace(-np.pi, np.pi, num_samples), force)(x))
            
        return forcing_functions

def gen_oscillator(model, speed=.05):
    with model:
        # ------------------ Oscillator -------------------
        osc = nengo.Ensemble(n_neurons=500, dimensions=2)
        # recurrent connections
        nengo.Connection(osc, osc,
                          transform=np.eye(2) + \
                                    np.array([[1, -1],
                                              [1, 1]]) * speed)
        return osc

def gen_point_attractor(model, in_goal, n_neurons=200, alpha=10, beta=10/4.):
    # create an ensemble with point attractor dynamics
    with model: 
        yz = nengo.Ensemble(n_neurons=n_neurons, dimensions=2, radius=5)
        # set up recurrent connection for system state, which 
        # specify the rate of change of each dimension represented.
        # first row of the transform is dyz, second row is ddyz
        nengo.Connection(yz, yz, 
                         transform=np.eye(2) + \
                                      np.array([[0, 1],
                                                [-alpha*beta, -alpha]]),
                         synapse=1)
        # connect up the input signal
        nengo.Connection(in_goal, yz[1],
                         transform=[[alpha*beta]], 
                         synapse=1)
        return yz

m = nengo.Network()
with m: 
    # --------------------- Inputs --------------------------
    in_goal1 = nengo.Node([0])
    zero = nengo.Node([0])

    # ------------------- Point Attractors --------------------
    yz1 = gen_point_attractor(m, in_goal1, n_neurons=1000)
    
    # -------------------- Oscillators ----------------------
    osc = gen_oscillator(m, speed=.05)
    
    # generate our forcing function
    y_des = np.load('heart_traj.npz')['arr_0']
    y_des[:, 0] = 0
    forcing_functions = gen_forcing_functions(y_des)
    
    def force(x, function, gain=1):
        # calculate the angle theta
        theta = np.arctan2(x[0], x[1])
        # decode our function off of the theta value
        return function(theta) * gain
    # connect oscillator to point attractor
    nengo.Connection(osc, yz1[1], function=lambda x: force(x, forcing_functions[0]))

    # output for easy viewing
    output = nengo.Ensemble(n_neurons=1, dimensions=2, neuron_type=nengo.Direct())
    nengo.Connection(yz1[0], output[0], synapse=.01)
    nengo.Connection(zero, output[1], synapse=.01)