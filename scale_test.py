import numpy as np
#from sklearn.preprocessing import MinMaxScaler
import ipdb

# sklearn didn't work
"""
scale = MinMaxScaler((-1, 1))
x = np.arange(9).reshape(3,3)
res = scale.fit_transform(x)
"""

"""
def scale_it(dat, axis=0, bounds=(-1, 1)):
    d_max = np.max(dat, axis=axis)
    d_min = np.min(dat, axis=axis)
    return (np.diff(bounds) * (dat - d_min)) / (d_max - d_min) + bounds[0]
"""

def d3_scale(dat, out_range=(-1, 1)):
    domain = [np.min(dat, axis=0), np.max(dat, axis=0)]

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

print(d3_scale(np.array([-2, 0, 2], dtype=np.float)))
print(d3_scale(np.array([-3, -2, -1], dtype=np.float)))
ipdb.set_trace()