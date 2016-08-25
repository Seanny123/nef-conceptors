import nengolib
import numpy as np


class TanhWithBias(nengolib.neurons.Unit):
    """A rectified linear neuron model."""

    def __init__(self, seed=0, leak_rate=0.6):
        super(TanhWithBias, self).__init__()
        self.leak_rate = leak_rate
        self.rng = np.random.RandomState(seed=seed)

    def gain_bias(self, max_rates, intercepts):
        """Return gain and bias given maximum firing rate and x-intercept."""
        bias = 0.8 * 2 * (self.rng.randn(*max_rates.shape) - 0.5)
        return np.ones_like(max_rates), bias

    def step_math(self, dt, J, output):
        """Compute rates in Hz for input current (incl. bias)"""
        output[...] = 0.6 * np.tanh(J)
