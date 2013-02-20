"""
Created Dec 22, 2012

This file sets the parameters to be used in sovling the value function

Matches the file ./BellmanEquationCode/SetParaStruc.m
"""
from numpy import array
import os

sl = os.path.sep


class DotDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, in_dct=None):
        """
        A dictionary that supports dot notation as well as dictionary access
        notation.

        Examples
        ---------
        >>> d = DotDict() or d = DotDict({'val1':'first'})
        >>> d.val2 = 'second'  # use dot notation to set attributes
        >>> d['val2'] = 'second'  # use dict notation to set attributes
        >>> d.val2  # use dot notation get attributes
        >>> d['val2']  # use dict notation get attributes
        """
        dct = in_dct if in_dct else dict()
        for key, value in dct.items():
            if hasattr(value, 'keys'):
                value = DotDict(value)
            self[key] = value

# 1. Paramters describing the preferences
theta_1 = 2  # type of Agent 1
theta_2 = 1  # Type of Agent 2
psi = .6  # Leisure consumption substitution
sigma = 2  # risk aversion
beta = .96  # subjective time discount factor
n1 = .5  # mass of type 1 agent
n2 = .5  # mass of type 2 agent

# 2. Technology
g_l = .15  # Government expenditure in low state s_l
g_h = .17  # Government expenditure in high state s_h
g = array([g_l, g_h])
P = array([[.5, .5], [.5, .5]])  # Transition Matrix for g shocks
alpha_1 = .5  # pareto wt of agent 1
alpha_2 = 1 - alpha_1  # pareto wt of agent 2
alpha_1 = alpha_1 * n1  # population adjusted pareto wts
alpha_2 = alpha_2 * n2  # population adjusted pareto wts
alpha = array([alpha_1, alpha_2])  # Pareto Weights for agent 1 and Agent 2
sSize = 2  # Dimension of the markov state

# 3. Others
ctol = 1e-8   # stopping criteria for iner optimization
grelax = .95  # wt of the new coeff NOTE: This was set at .95 prior to meddling to try for convergence
Niter = 500  # number of value function iterations
resolve_ctr = 1  # Frequency with which the routine for unresolved points must be tried
NumSim = 10000  # Number of simulations
btild_1 = 0  # initial condition for Time0 problem
DeltaX = 1  # deviation from steadys state for the default grid
DeltaR = 0.5

# 4. Polynomial Approximation details
ApproxMethod = 'spli'  # basis poly
xGridSize = 20  # density of grid in dimension x
RGridSize = 20  # density of grid in dimension R
orderofappx_x = 19  # number of splines in dimension x
orderofappx_R = 19  # number of splines in dimension R

# 5. Path information
root_dir_temp = os.getcwd()
num_extras = len(root_dir_temp.split(sl)[-1])
if root_dir_temp[-num_extras:] == 'gspy':
    root_dir = root_dir_temp + sl
else:
    root_dir = root_dir_temp[:-len(root_dir_temp.split(sl)[-1])]
texpath = root_dir + 'tex' + sl
plotpath = root_dir + 'graphs' + sl
datapath = root_dir + 'data' + sl

new_paths = [texpath, plotpath, datapath]

for path in new_paths:
    if not os.path.exists(path):
        os.makedirs(os.makedirs)

# create parameters dictionary
params = DotDict()
params.sl = sl
params.ctol = ctol
params.theta_1 = theta_1
params.theta_2 = theta_2
params.psi = psi
params.sigma = sigma
params.beta = beta
params.g_l = g_l
params.g_h = g_h
params.g = g
params.P = P
params.alpha_1 = alpha_1
params.alpha_2 = alpha_2
params.n1 = n1
params.n2 = n2
params.Niter = Niter
params.sSize = sSize
params.xGridSize = xGridSize
params.RGridSize = RGridSize
params.rootpath = root_dir
params.texpath = texpath
params.plotpath = plotpath
params.datapath = datapath
params.ApproxMethod = ApproxMethod
params.orderofappx_R = orderofappx_R
params.orderofappx_x = orderofappx_x
params.grelax = grelax
params.resolve_ctr = resolve_ctr
params.NumSim = 10000
params.btild_1 = btild_1
params.DeltaX = DeltaX
params.DeltaR = DeltaR

# Add path info to params
params.root_dir = root_dir
params.texpath = texpath
params.plotpath = plotpath
params.datapath = datapath
