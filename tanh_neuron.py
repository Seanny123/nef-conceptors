import nengolib
import numpy as np
from nengo.builder import Builder, Operator, Signal
from nengo.builder.neurons import SimNeurons


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

    def step_math(self, dt, J, output, last_rate):
        """Compute rates in Hz for input current (incl. bias)"""
        output[...] = (1 - self.leak_rate) * last_rate + self.leak_rate * np.tanh(J)
        last_rate = output


@Builder.register(TanhWithBias)
def build_tanhwithbias(model, tanhwithbias, neurons):
    """acquire the last_rate argument"""

    model.sig[neurons]['last_rate'] = Signal(
        np.zeros(neurons.size_in), name="%s.last_rate" % neurons)
    model.add_op(SimNeurons(neurons=tanhwithbias,
                            J=model.sig[neurons]['in'],
                            output=model.sig[neurons]['out'],
                            states=[model.sig[neurons]['last_rate']]))
