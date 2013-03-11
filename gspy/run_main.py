"""
Created Dec 22, 2012

Author(s): Spencer Lyon and Chase Coleman

Solve the G-S economy with BGP preferences of the form
 psi.c^(1-sigma)/1-sigma+(1-psi)log[1-l] with following calibrations

1. The ratio of productivities is 3.3 and the low productivity is
   normalized to 1
2. psi is choosen to get a average FE of about 0.5
3. pareto wts are such that the no-shock problem gives a tax rate of
   about 20 percent
4. Government expenditures are about 11 and 13 percent of the output
5. beta =0.9

This file mimics ./Main/RunMainWithAltSigmas.m
"""
from set_params import params
import bellman as bell
# import cyed.bellmancy as bell
# import bellman_parallel as bell
import numpy as np
from scipy.optimize import fsolve

sl = params.sl


def get_calibration_fe(x, target, th1, th2, tau, g_Y, n1, n2):
    """
    Mimics the file ./Main/GetCalibrationFrischElasticity.m
    """
    gamma = x[0]
    Y = x[1]

    trans = ((tau - g_Y) / (n1 + n2)) * Y
    l1Num = (1 - tau) * theta_1 - gamma * trans
    l2Num = (1 - tau) * theta_2 - gamma * trans
    l1Den = (1 + gamma) * (1 - tau) * theta_1
    l2Den = (1 + gamma) * (1 - tau) * theta_2
    l1 = l1Num / l1Den
    l2 = l2Num / l2Den

    AvgFE = (n1 * (1 / l1 - 1) + n2 * (1 / l2 - 1)) / (n1 + n2)
    res = np.empty(2)
    res[0] = Y - theta_1 * l1 * n1 - theta_2 * l2 * n2
    res[1] = (AvfFETarget - AvgFE)

    return res

theta_1 = 3.3  # theta high
theta_2 = 1   # theta low
g_l_y = .11  # g low
g_h_y = .13  # g high
n1 = 1
n2 = 1
tau = .2
g_Y = np.mean([g_l_y, g_h_y])
AvfFETarget = .5
z = fsolve(get_calibration_fe, (1, 1),
           args=(AvfFETarget, theta_1, theta_2, tau, g_Y, n1, n2), xtol=1e-13)
gamma = z[0]
Y = z[1]

# BASELINE GOVERNMENT EXsPENDITURE LEVELS
g = g_Y * Y

# BASELINE PSI
psi = 1 / (1 + gamma)
# BASELINE DISCOUNT FACTOR

beta = .9

# BASELINE PARETO WTS
alpha_1 = 0.69
alpha_2 = 1 - alpha_1
params.n1 = n1
params.n2 = n2
alpha_1 = alpha_1 * params.n1
alpha_2 = alpha_2 * params.n2

# BASELINE PROBABILITY MATRIX
NewPh = .5
params.P = np.array([[1 - NewPh, NewPh], [1 - NewPh, NewPh]])

# POPULATE THE PARA STRUC WITH THE BASELINE VALUES
params.beta = .9
params.alpha_1 = alpha_1
params.alpha_2 = alpha_2
params.psi = psi
params.g = np.array([g_l_y, g_h_y]) * Y
params.theta_1 = theta_1
params.theta_2 = theta_2
params.btild_1 = 0
params.alpha_1 = alpha_1
params.alpha_2 = alpha_2
params.datapath = params.root_dir + 'data' + sl + 'temp' + sl
casename = 'sigma'
params.StoreFileName = 'c' + casename + '.mat'
coeff_file_name = params.datapath + params.StoreFileName

#  --- SOLVE THE BELLMAN EQUATION --------------------------------------
params.Niter = 200  # MAXIMUM NUMBER OF ITERATION


# flagSetRGrid,flagSetxGrid
# TAKES TWO VALUES: 0 IF DEFAULT GRID OR 1 FOR USERDEFINED GRID

params.flagSetRGrid = 1
params.flagSetxGrid = 1
params.xMin = -2.5
params.xMax = 2.5

 # EXPERIMENT 1: SIGMA=1
casename = 'sigmaLow'
params.StoreFileName = 'c' + casename + '.mat'
coeff_file_name = params.datapath + params.StoreFileName
params.sigma = 1
params.RMin = 2.2
params.RMax = 3.5
bell.main(params)

# # EXPERIMENT 2: SIGMA=2
# casename = 'sigmaMed'
# params.StoreFileName = 'c' + casename + '.mat'
# coeff_file_name = params.datapath + params.StoreFileName
# params.sigma = 2
# params.RMin = 3.5
# params.RMax = 4.5
# main(params)
# # bellman.main(params)

# # EXPERIMENT 3: SIGMA=3
# casename = 'sigmaHigh'
# params.StoreFileName = 'c' + casename + '.mat'
# coeff_file_name = params.datapath + params.StoreFileName
# params.sigma = 3
# params.RMin = 4.5
# params.RMax = 5.5
# main(params)
# # bellman.main(params)

#----------------------------Simulate the Model--------------------------------#
# NumSim = 60000
# rHist0 = np.random.rand(NumSim)
# K = 3

#Need to create the three cases.  Create file names for each of the coefficients
#That are given for the policy function.  Then need to pass in appropriate inputs
#to runsimulation and save results.

## NOTE: finished through line 138 in RunMainWithAltSigmas.m file
