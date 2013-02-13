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
import numpy as np
import scipy.linalg as la
import scipy.optimize as opt
from scipy.io import savemat
# from compeconpy import fundefn, funfitxy, funeval
from compeconcy import fundefn, funfitxy, funeval_new
from steady_state import steady_state_res, find_steady_state
from inner_opt import uAlt, check_grad
from set_params import DotDict

cimport cython
cimport numpy as np

DTYPE = np.float

ctypedef np.float_t DTYPE_t

@cython.boundscheck(False)
def main(object params):
    """
    This is the main file computes the value function via time iteration for
    the parameters passsed in the structure Para.

    NOTATION:
    --------
    x = u_2 btild
    R = u_2/u_1
    """
    #BUILD GRID
    cdef object info_dict
    
    params, info_dict = build_grid(params)
    print('Msg: Completed definition of functional space')

    #INITIALIZE THE COEFF
    cdef np.ndarray domain, c, policy_rules_store
    
    print('Msg: Initializing the Value Function...')
    [domain, c, policy_rules_store] = init_coef(params.copy(), info_dict)
    print('Msg: ... Completed')

    #Iterate on the Value Function
    #This block iterates on the bellman equation
    cdef np.ndarray x_slice, R_slice, s_slice
    cdef int grid_size
    
    x_slice = domain[:, 0]
    R_slice = domain[:, 1]
    s_slice = domain[:, 2]
    grid_size = params.GridSize

    #Initialize The Sup Norm Error matrix and variables needed in loop
    cdef np.ndarray errorinsupnorm, policy_rules_old
    
    errorinsupnorm = np.ones(params.Niter)
    policy_rules_old = np.zeros(policy_rules_store.shape)

    #Begin the for loops
    # for it in xrange(1, params.Niter):
    cdef int it, ctr
    for it in xrange(1, 11):
        #Record Start Time.  Total time will be starttime-endtime
    cdef np.ndarray ix_unsolved, ix_solved
        # Clear index arrays to be sure they don't persist to long
        try:  # On first iteration they aren't defined yet.
            del ix_unsolved
            del ix_solved
        except:
            pass

        start_time = time.time()
        
        cdef np.ndarray exitflag, vnew
        
        exitflag = np.zeros(int(grid_size))
        vnew = np.zeros(int(grid_size))

        #Clear Records of arrays that store index of failure
        #^We can add if we need to

        #Initialize the initial guess for the policy rules that the inneropt
        #will solve
        policy_rules_old = policy_rules_store
        
        # Do I need to redefine ctr as an int?
        cdef int ctr
        
        for ctr in xrange(grid_size // 2):
            # Here they use a parfor loop invoking parallel type for loops
            # Will make this parallel when we speed up program
            cdef double x, R, v_new
            cdef int s, flag
            
            x = x_slice[ctr]
            R = R_slice[ctr]
            s = s_slice[ctr] - 1

            # Initalize Guess for the inneropt
            cdef np.ndarray xInit, policyrules

            xInit = policy_rules_store[ctr, :]

            # Inner Optimization
            policyrules, v_new, flag = check_grad(x, R, s, c, info_dict,
                                                  xInit, params)
            vnew[ctr] = v_new
            exitflag[ctr] = flag

            if flag == 1:
                policy_rules_store[ctr, :] = policyrules

        #---IID Case-----#
        # In the IID case we solve it for s = 1 and use the solution
        # to populate s = 2
        exitflag[grid_size // 2:] = exitflag[:grid_size // 2]
        vnew[grid_size // 2:] = vnew[:grid_size // 2]
        policy_rules_store[grid_size // 2:] = \
                                        policy_rules_store[:grid_size // 2]
        #----------------#

        #------------------Beigns HandleUnresolvedPoints.m-------------------#
        #--------------------------------------------------------------------#

        # Locate the unresolved points
        # Defined earlier at the try statement
        ix_solved = np.where(exitflag == 1)[0]
        ix_unsolved = np.where(exitflag != 1)[0]

        print 'The fraction of nodes that remain unsolved at first pass is ',\
                ix_unsolved.size / exitflag.size

        # Resolve the FOC at the failed points
        # Do I need to redefine it?
        cdef int it, num_trials
        
        if it % params.resolve_ctr == 0:
            num_trials = 5

            if ix_unsolved.size != 0:
                print 'Points that failed the first round of FOC \n', \
                                                        domain[ix_unsolved, :]
                print 'Resolving the unresolved points using alternate routine'

            #----------------Begins UnResolvedPoints.m-----------------------#
            #----------------------------------------------------------------#
                print 'Unresolved so far ', ix_unsolved.size
            # TODO: Need to doulbe check a couple of these
            # ref_id and p_r_init
            cdef int num_unsolved, i, uns_index, _s
            cdef double x, R, dist, x_ref, ref_id
            cdef np.ndarray x_store, x_target, p_r_store, x0, p_r_init
            
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

                cdef int tr_indx, numresolved

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
        
        cdef np.ndarray ix_solved_1, ix_solved_2
        
        ix_solved_1 = ix_solved[ix_solved <= (grid_size // params.sSize - 1)]
        ix_solved_2 = ix_solved[ix_solved > (grid_size // params.sSize - 1)]

        #-------------------Ends HandleUnresolvedPoints.m--------------------#
        #--------------------------------------------------------------------#
        #-------------------Begins UpdateCoefficients.m----------------------#
        #--------------------------------------------------------------------#
        cdef np.ndarray c_temp, c_new, c_diff, c_old
        cdef unsigned junk
        
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

        save_name = params.datapath + params.StoreFileName[1:-4] + params.sl +\
                    'c_' + str(it) + '.mat'
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

        savemat(save_name, data)
        #-------------------Ends UpdateCoefficients.m------------------------#
        #--------------------------------------------------------------------#

        if ix_unsolved.size / grid_size > 0.02:
            print 'Exiting for a new grid'
            break

        if errorinsupnorm[it - 1] < params.ctol:
            print 'Convergence criterion met. Exiting bellman.main'
            break

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

    return it
