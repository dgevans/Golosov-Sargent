"""
This is the equivalent to MainBellman.m from the original code
This is the main file that computes the value function via
time iteration for the parameters passed in the structure parameters

Notation:
    x = u_2 btild
    R = u_2/u_1
"""
from __future__ import division
from itertools import product
import time
from mpi4py import MPI
import numpy as np
import scipy.optimize as opt
from scipy.io import savemat
# from compeconpy import fundefn, funfitxy, funeval, funeval_new
from cyed.compeconcy import fundefn, funfitxy, funeval_new
from steady_state import steady_state_res, find_steady_state
from inner_opt import uAlt, check_grad
from set_params import DotDict

# Not sure what to do with the section DEFAULT PARAMETERS

comm = MPI.COMM_WORLD
rank = comm.Get_rank()  # Which process we are on
size = comm.Get_size()  # Total number of processes


def pprint(x):
    if rank == 0:
        print(x)


def mpi_tuples(data):
    """
    This function takes the data vector or matrix and returns the send
    counts, displacement counts, and a local displacement array to be
    used by MPI.

    TODO: this docstring.
    """
    # Calculate send counts and displacements for Allgatherv
    rows = data.shape[0]

    # Calculate the send counts
    counts = np.zeros(size)
    counts[:] = rows // size
    counts[range(rows % size)] += 1
    count_tuple = tuple(counts * data.shape[1]) if data.ndim == 2 \
                                                else tuple(counts)

    # Calculate the displacements
    disps = counts.cumsum() - counts[0]
    if (rows % size) != 0:
        disps[-(size - rows % size):] += 1
    disp_tuple = tuple(disps * data.shape[1]) if data.ndim == 2 \
                                                else tuple(disps)

    # append extra element so last process has an ending range
    loc_disps = np.append(disps, [rows])

    return count_tuple, disp_tuple, loc_disps


#Build Grid
def build_grid(params):
    '''
    This is the function that executes the equivalent of BuildGrid.m.
    This function defines the grid and defines the value function.  There
    are two alternatives.  First, the user could input either the x or Rgrid.
    This should supercede any other option.  Otherwise we use the SS
    computation and use efault DeltaX,DeltaR parameters to set the deviation
    from the SS.

    params.flagSetxGrid sets the flag for either using the default grid (Value = 0)
     or using the user defined grid (Value = 1)

    params.flagsetRdgrid sets the flag for either using the default grid (Value = 0)
    or using the user defined grid (Value = 1)
    '''
    #???Do we want to name all of the elements in params like params[0] = ___
    #How are we importing it?

    #Find the SS
    [xSS, RSS, NA] = find_steady_state(0, 3, params)

    #Check the flags and define Grid Endpoints
    if 'flagSetxGrid' in params:
        flagSetxGrid = params.flagSetxGrid
    else:
        flagSetxGrid = 0

    if flagSetxGrid == 1:
        xMin = params.xMin
        xMax = params.xMax
        pprint('Msg: Using user defined grid on x')
        #^Will proceed using user defined gridpoints
    else:
        xMin = xSS - params.DeltaX
        xMax = xSS + params.DeltaX
        pprint('Msg: Using default grid around SS')
        #^Will proceed using the default gridpoints
    #Uniformly space the gridpoints
    xGrid = np.linspace(xMin, xMax, params.xGridSize, endpoint=True)

    #Update the params struct
    params.xGrid = xGrid
    params.xLL = xMin
    params.xUL = xMax

    #Check the flags and define Grid Endpoints
    if 'flagSetRGrid' in params:
        flagSetRGrid = params.flagSetRGrid
    else:
        flagSetRGrid = 0

    if flagSetRGrid == 1:
        RMin = params.RMin
        RMax = params.RMax
        pprint('Proceeding with RGrid based on user inputs')
        #^Will proceed using user defined grid
    else:
        RMin = RSS - params.DeltaR
        RMax = RSS + params.DeltaR
        pprint('Proceeding with default RGrid')
        #^Will proceed using default grid

    #Uniformly space the gridpoints
    RGrid = np.linspace(RMin, RMax, params.RGridSize, endpoint=True)
    params.RGrid = RGrid
    #Gives the total gridsize
    GridSize = params.xGridSize * params.RGridSize * params.sSize

    #Update params
    params.GridSize = GridSize
    params.xMin = xMin
    params.xMax = xMax
    params.RMax = RMax
    params.RMin = RMin

    ##Define Functional Space
    #Priorly used the CompEcon Library routine 'fundefn' to create a functional
    #space which stores key settings for the basis polynomials, domand and nodes.
    #V(s) is he functional space for the value function given the discrete shock
    #Need to look up documentation on Compecon toolbox function fundefn

    V = np.zeros(2, dtype=object)
    V[0] = fundefn(np.array([params.orderofappx_x, params.orderofappx_R]),
                   np.array([xMin, RMin]), np.array([xMax, RMax]))
    V[1] = fundefn(np.array([params.orderofappx_x, params.orderofappx_R]),
                   np.array([xMin, RMin]), np.array([xMax, RMax]))

    #It seems that we can just use our own grid that was built using linspace
    #One difference in the grid that is created by fundef is that it would appear
    #they account for degrees of freedom of some sort b/c the number of gridpoints
    #is n + 1 - k.  Don't think we need to worry about it much, but it's worth noting

    return params, V


def init_coef(params, info_dict):
    """
    Mimics ./BellmanEquationCode/InitializeCoeff.m

    INITIALIZE THE COEFF
    This section uses the stationary policies to initialze the value
    functions. The block prdoduces three outcomes
    1. domain which is the vectorized domain
    2. PolicyRulesStore : This serves as a matrix of guess for optimal
    policies at the interpolation nodes
    3. c0 : initial coeffecients
    """
    p = DotDict(params.copy())
    xGrid = p.xGrid
    RGrid = p.RGrid
    gTrue = p.g
    p.g = gTrue.mean() * np.ones((2, 1))

    #Need to initialize arrays before we fill them.
    n_size = p.xGridSize * p.RGridSize
    n_coefs = (p.xGridSize - 1) * (p.RGridSize - 1)
    c0 = np.zeros((p.sSize, n_coefs))
    V0 = np.zeros((p.sSize, n_size))
    xInit_0 = np.zeros((2, p.xGridSize * p.RGridSize, 7))

    _domain = np.array(list(product(p.xGrid, p.RGrid)))
    for _s in range(p.sSize):
        n = 0
        if _s == 0:
            for xctr in xrange(p.xGridSize):
                for Rctr in xrange(p.RGridSize):
                    _x = xGrid[xctr]
                    _R = RGrid[Rctr]

                    # NOTE: _domain is product(xGrid, RGrid). See above loop
                    # _domain[_s, n, :] = [_x, _R]

                    # Initialize the guess for Stationary Policies
                    cRat = _R ** (-1. / p.sigma)

                    c1_1 = (0.8 * (p.n1 * p.theta_1 + p.n2 * p.theta_2) - p.g[0])\
                                / (p.n1 + cRat * p.n2)

                    c1_2 = (0.8 * (p.n1 * p.theta_1 + p.n2 * p.theta_2) - p.g[1])\
                                / (p.n1 + cRat * p.n2)

                    c2_1 = cRat * c1_1

                    guess = np.array([c1_1, c1_2, c2_1]).flatten()
                    [xSS, info, exitFlag, msg] = opt.fsolve(steady_state_res,
                                                       guess, full_output=1,
                                                        args=(_x, _R, p, _s, True),
                                                        xtol=1.0e-13)

                    [res, c1_, c2_, l1_, l2_] = steady_state_res(xSS, _x, _R,
                                                                 p, _s)

                    # Present Discounted value for Stationary policies
                    V0[_s, n] = (p.alpha_1 * uAlt(c1_, l1_, p.psi, p.sigma) +
                                p.alpha_2 * uAlt(c2_, l2_, p.psi, p.sigma)).dot( \
                                p.P[_s, :].T) / (1 - p.beta)

                    xInit_0[_s, n, :] = [c1_[_s], c2_[_s], l1_[_s],
                                         l2_[_s], _x, _R, _x]

                    n += 1

            c0[_s, :], B = funfitxy(info_dict[_s], _domain, V0[_s, :])
        else:
            c0[_s, :] = c0[_s - 1, :]
            V0[_s, :] = V0[_s - 1, :]
            xInit_0[_s, :, :] = xInit_0[_s - 1, :, :]

    domain_1 = np.column_stack((_domain, np.ones((n_size, 1))))
    domain_2 = np.column_stack((_domain, np.ones((n_size, 1)) * 2))
    domain = np.concatenate((domain_1, domain_2), axis=0)

    c = c0
    data = {'c': c}
    savemat(p.datapath + 'c1.mat', data, oned_as='row')

    x_shape = xInit_0.shape
    policy_rules_store = np.empty((x_shape[0] * x_shape[1], x_shape[2] * 2)).T

    for i in range(x_shape[2]):
        cols = range(2 * i, 2 * i + 2)
        policy_rules_store[cols, :400] = xInit_0[0, :, i].squeeze()
        policy_rules_store[cols, 400:] = xInit_0[1, :, i].squeeze()

    policy_rules_store = policy_rules_store.T

    # NOTE: this is just here for comparison with the MatLab objects
    # data = {'my_c0': c0, 'my_v0': V0, 'my_xinit': xInit_0,
    #         'my_domain': domain, 'my_policyrules': policy_rules_store}
    # savemat('/Users/spencerlyon2/Documents/Research/Golosov-Sargent/' + \
    #           'gspy/data/debugging/init_coefs.mat', data)

    # NOTE: All objects in data are within 1e-13 of corresponding MatLab objs

    return [domain, c, policy_rules_store]


def main(params):
    """
    This is the main file computes the value function via time iteration for
    the parameters passsed in the structure Para.

    NOTATION:
    --------
    x = u_2 btild
    R = u_2/u_1
    """
    #BUILD GRID
    # NOTE: we are not doing this only on root process because you can't
    #       pickle params and mpi4pi uses pickling behind the scenes.
    #       Also, it doesn't really take much time.
    params, info_dict = build_grid(params)
    grid_size = int(params.GridSize)
    pprint('Msg: Completed definition of functional space')

    if rank == 0:
        #INITIALIZE THE COEFF
        pprint('Msg: Initializing the Value Function...')
        [domain, c, policy_rules_store] = init_coef(params.copy(), info_dict)
        pprint('Msg: ... Completed')

        #Iterate on the Value Function
        #This block iterates on the bellman equation
        x_slice = domain[:, 0]
        R_slice = domain[:, 1]
        s_slice = domain[:, 2]

        #Initialize The Sup Norm Error matrix and variables needed in loop
        errorinsupnorm = np.ones(params.Niter)
        # policy_rules_old = np.zeros(policy_rules_store.shape)
    else:
        policy_rules_store = None
        c = None
        x_slice = None
        R_slice = None
        s_slice = None
        errorinsupnorm = None

    comm.Barrier()  # Pause here to wait for process 1 to initialize vf
    policy_rules_store = comm.bcast(policy_rules_store, 0)
    c = comm.bcast(c, 0)
    x_slice = comm.bcast(x_slice, 0)
    R_slice = comm.bcast(R_slice, 0)
    s_slice = comm.bcast(s_slice, 0)
    errorinsupnorm = comm.bcast(errorinsupnorm, 0)
    comm.Barrier()  # Pause here to make sure all processes have the data.

    #Begin the for loops
    for it in xrange(1, params.Niter):
        # Clear index arrays to be sure they don't persist to long
        try:  # On first iteration they aren't defined yet.
            del ix_unsolved
            del ix_solved
        except:
            pass

        start_time = time.time()

        exitflag = np.zeros(int(grid_size))
        vnew = np.zeros(int(grid_size))

        vec_count_tuple, vec_disp_tuple, vec_loc_disps = mpi_tuples(vnew)
        mat_count_tuple, mat_disp_tuple, mat_loc_disps = mpi_tuples(
                                                           policy_rules_store)

        my_start = int(vec_loc_disps[rank])
        my_end = int(vec_loc_disps[rank + 1])

        local_vnew = np.zeros((vec_count_tuple[rank]))
        local_flag = np.zeros((vec_count_tuple[rank]))
        local_prs = np.zeros((vec_count_tuple[rank],
                              policy_rules_store.shape[1]))
        local_ind = 0

        for ctr in xrange(my_start, my_end):
            # Here they use a parfor loop invoking parallel type for loops
            # Will make this parallel when we speed up program
            x = x_slice[ctr]
            R = R_slice[ctr]
            s = int(s_slice[ctr] - 1)

            # Initalize Guess for the inneropt
            xInit = policy_rules_store[ctr, :]

            # Inner Optimization
            policyrules, v_new, flag = check_grad(x, R, s, c, info_dict,
                                                  xInit, params)
            local_vnew[local_ind] = v_new
            local_flag[local_ind] = flag

            if flag == 1:
                local_prs[local_ind, :] = policyrules

            local_ind += 1

        # Gather vnew
        comm.Gatherv(local_vnew, [vnew, vec_count_tuple,
                                  vec_disp_tuple, MPI.DOUBLE])

        # Gather exitflag
        comm.Gatherv(local_flag, [exitflag, vec_count_tuple,
                                  vec_disp_tuple, MPI.DOUBLE])

        # Gather policy_rules_store
        comm.Gatherv(local_prs, [policy_rules_store, mat_count_tuple,
                                 mat_disp_tuple, MPI.DOUBLE])
        #------------------Beigns HandleUnresolvedPoints.m-------------------#
        #--------------------------------------------------------------------#

        # Only do this stuff on root process
        if rank == 0:
            # Locate the unresolved points
            ix_solved = np.where(exitflag == 1)[0]
            ix_unsolved = np.where(exitflag != 1)[0]

            print 'The fraction of nodes that remain unsolved at first pass is ',\
                    ix_unsolved.size / exitflag.size

            # Resolve the FOC at the failed points
            if it % params.resolve_ctr == 0:
                num_trials = 5

                if ix_unsolved.size != 0:
                    print 'Points that failed the first round of FOC \n', \
                                                            domain[ix_unsolved, :]
                    print 'Resolving the unresolved points using alternate routine'

                #----------------Begins UnResolvedPoints.m-----------------------#
                #----------------------------------------------------------------#
                    print 'Unresolved so far ', ix_unsolved.size
                num_unsolved = ix_unsolved.size
                for i in xrange(num_unsolved):
                    ix_solved = np.where(exitflag == 1)[0]  # Reset solved points
                    uns_index = ix_unsolved[i]  # Start on first unresolved point

                    x = x_slice[uns_index]
                    R = R_slice[uns_index]
                    _s = s_slice[uns_index] - 1
                    print 'Resolving... ', [x, R, _s]

                    #Try 1
                    #----------------Begins GetInitialApproxPolicy.m-------------#
                    #------------------------------------------------------------#
                    x_store = domain[ix_solved, :]  # domain for all good points
                    x_target = np.array([x, R, _s])  # Current unsolved point
                    p_r_store = policy_rules_store[ix_solved, :]  # PRS all good

                    dist = ((x_store - x_target) ** 2).sum(1)
                    ref_id = np.argmin(dist)
                    p_r_init = p_r_store[ref_id]
                    x_ref = x_store[ref_id]
                    #----------------Ends GetInitialApproxPolicy.m---------------#
                    #------------------------------------------------------------#

                    x0 = np.zeros((2, num_trials))
                    x0[0, :] = np.linspace(x_ref[0], x, num_trials)
                    x0[1, :] = np.linspace(x_ref[1], R, num_trials)

                    for tr_indx in xrange(num_trials):
                        policyrules, v_new, flag = check_grad(x0[0, tr_indx],
                                                            x0[1, tr_indx], _s, c,
                                                            info_dict, p_r_init,
                                                            params)
                        p_r_init = policyrules

                    if flag != 1:
                        #Try method 2 if 1 doesn't work
                        xguess2 = np.array([x, R, _s]) * (1 +
                                np.sign(np.array([x, R, _s]) - x_ref) * .05)
                        #--------------Begins GetInitialApproxPolicy.m-----------#
                        #--------------------------------------------------------#
                        x_target = xguess2  # Current unsolved point
                        dist = ((x_store - x_target) ** 2).sum(1)
                        ref_id = np.argmin(dist)
                        p_r_init = p_r_store[ref_id]
                        x_ref = x_store[ref_id]
                        #--------------Ends GetInitialApproxPolicy.m-------------#
                        #--------------------------------------------------------#
                        x0 = np.zeros((2, num_trials))
                        x0[0, :] = np.linspace(x_ref[0], x, num_trials)
                        x0[1, :] = np.linspace(x_ref[1], R, num_trials)

                        for tr_indx in xrange(num_trials):
                            policyrules, v_new, flag = check_grad(x0[0, tr_indx],
                                                            x0[1, tr_indx], _s, c,
                                                            info_dict, p_r_init,
                                                            params)
                            p_r_init = policyrules

                    exitflag[uns_index] = flag
                    vnew[uns_index] = v_new
                    if flag == 1:
                        policy_rules_store[uns_index, :] = policyrules

                    # NOTE: Why do we this at top and bottom of the for loop?
                    ix_solved = np.where(exitflag == 1)[0]

                ix_solved = np.where(exitflag == 1)[0]
                ix_unsolved = np.where(exitflag != 1)[0]
                numresolved = num_unsolved - ix_unsolved.size

                if ix_unsolved.size != 0 or numresolved != 0:
                    print 'Number of points solved with alternative guess ', \
                                                    numresolved

                #-------------------Ends UnResolvedPoints.m----------------------#
                #----------------------------------------------------------------#

            if numresolved - num_unsolved != 0:
                NumTrials = 10
                print 'Resolving the unresolved points using alternate routine'

            ix_unsolved = np.where(exitflag != 1)[0]
            ix_solved = np.where(exitflag == 1)[0]

            # TODO: Verify index in the two lines below (should I do -1 or not?)
            ix_solved_1 = ix_solved[ix_solved <= (grid_size // params.sSize - 1)]
            ix_solved_2 = ix_solved[ix_solved > (grid_size // params.sSize - 1)]

            #-------------------Ends HandleUnresolvedPoints.m--------------------#
            #--------------------------------------------------------------------#
            #-------------------Begins UpdateCoefficients.m----------------------#
            #--------------------------------------------------------------------#
            c_temp, junk = funfitxy(info_dict[0],
                                    domain[ix_solved_1, :2],
                                    vnew[ix_solved_1])

            c_new = np.tile(c_temp, (2, 1))
            if it == 1:  # Create c_diff on first iteration
                c_diff = np.abs(c - c_new).sum(0)
            else:  # Append to it on all other iterations
                c_diff = np.row_stack([c_diff, np.abs(c - c_new).sum(0)])

            c_old = c

            # Take convex combination to update coefficients
            c = c_new * params.grelax + (1 - params.grelax) * c_old

            # Calculate error in supremum norm
            errorinsupnorm[it - 1] = np.max(np.abs(vnew[ix_solved_1] -
                                                funeval_new(c_old[0, :], info_dict[0],
                                                domain[ix_solved_1, :2],
                                                np.array([0, 0]))))

            end_time = time.time()
            elapsed = end_time - start_time
            print 'Completed iteration %i. Time required %.3f. Norm %.4e: ' % \
                                            (it, elapsed, errorinsupnorm[it - 1])

            save_path = params.datapath + params.StoreFileName[1:-4] + params.sl
            save_name = save_path + 'c_' + str(it) + '.mat'

            data = {'c': c,
                    'errorinsupnorm': errorinsupnorm,
                    'c_diff': c_diff,
                    'ix_solved': ix_solved,
                    'ix_unsolved': ix_unsolved,
                    'policy_rules_store': policy_rules_store,
                    'vnew': vnew,
                    'domain': domain,
                    'params': params,
                    'info_dict': info_dict}

            savemat(save_name, data, oned_as='row')
            #-------------------Ends UpdateCoefficients.m------------------------#
            #--------------------------------------------------------------------#

            if ix_unsolved.size / grid_size > 0.02:
                print 'Exiting for a new grid'
                break

            if errorinsupnorm[it - 1] < params.ctol:
                print 'Convergence criterion met. Exiting bellman.main'
                break

        comm.Barrier()  # Pause all processes here to wait for root
        c = comm.bcast(c, 0)
        policy_rules_store = comm.bcast(policy_rules_store, 0)
        comm.Barrier()  # pause to make sure everyone has the data

    if rank == 0:  # Only save data on root process
        file_name = params.datapath + params.StoreFileName
        data = {'c': c,
                    'errorinsupnorm': errorinsupnorm,
                    'c_diff': c_diff,
                    'ix_solved': ix_solved,
                    'ix_unsolved': ix_unsolved,
                    'policy_rules_store': policy_rules_store,
                    'vnew': vnew,
                    'domain': domain,
                    'params': params,
                    'info_dict': info_dict}

        savemat(file_name, data)

    return
