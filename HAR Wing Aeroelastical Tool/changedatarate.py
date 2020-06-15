import numpy as np
from scipy.interpolate import interp1d


def changedatarate(t, X, dt):
    tfinal = t[t.shape[0] - 1]
    tnew = np.arange(0, tfinal + dt, dt)
    Xnew = interp1d(t, X)(tnew)

    return tnew, Xnew
