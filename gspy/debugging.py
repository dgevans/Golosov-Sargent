"""
This file simply contains routines used to debug/test the python code
relative to the MatLab. Mostly it will just setup the workspace and
define appropriate variables.

It uses the sys.argv routine to retrieve a command line argument from
the user that specifies which routine needs to be called to set things
up

This file is called from the ipython command-line. It is done like this:

    run debugging.py 'init_coefs'

where 'init_coefs' is replaced with the phrase associated with the
name of the routine this file will be setting the workspace for.
"""
import sys

# The first argument we care about is the type of workspace that needs to be
# setup. This line gets the command line arguments and drops the first one
# becuase the first one is always file name (in this case debugging.py).
# The line below is necessary when calling this from ipython run magic b/c
# when called this way the first argument is the file name
args = sys.argv[1:]

if args[0] == 'init_coefs':
    print 'Initializing workspace for init_coefs'
    import numpy as np
    from scipy.io import loadmat
    from compeconpy import fundefn
    from set_params import DotDict

    a = np.array([-2.5,  2.2])
    b = np.array([2.5,  3.5])
    n = np.array([19, 19])
    _domain = loadmat('data/debugging/domain_them.mat')['domain_']
    V0 = loadmat('data/debugging/V0_them.mat')['V0']
    info_dict = fundefn(n, a, b)
    order = 0

    print 'Done. Initialized variables a, b, n, _domain, V0, info_dict, order, \
                DotDict'
else:
    raise ValueError('Workspace argument unknown.')
