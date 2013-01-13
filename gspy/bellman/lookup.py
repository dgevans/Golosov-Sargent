"""
Just a function that matches lookup.c or lookup.m from compecon
"""
import numpy as np
import pandas as pd

def lookup(tabvals, x, endadj):
    """
    Just a function that matches lookup.c or lookup.m from compecon
    """
    import numpy as np
    import pandas as pd
    n = np.prod(x.shape)
    m = tabvals.size
    if endadj >= 2:
        m = m - np.where(tabvals==tabvals[-1])[0].size

    temp_series = pd.Series(np.append(tabvals[:m], x))
    temp_series.sort()
    ind = temp_series.index.values
    temp = np.where(ind > m - 1)[0]
    j = ind[temp] - m
    ind = (temp - range(1, n+1)).reshape(x.shape)
    ind[j] = ind[:]
    if endadj == 1 or endadj == 2:
        ind[ind==0] = np.where(tabvals==tabvals[0])[0].size
    return ind
