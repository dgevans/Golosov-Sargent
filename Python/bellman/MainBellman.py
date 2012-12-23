#This is the equivalent to MainBellman.m from the original code
#This is the main file that computes the value function via
#time iteration for the parameters passed in the structure parameters

#Notation
# x = u_2 btild
# R = u_2/u_1

#Imports would go here, I think the first file will make all necessary imports
#for the process.

##Not sure what to do with the section DEFAULT PARAMETERS

#Build Grid
def Buildgrid(params):
    '''This is the function that executes the equivalent of BuildGrid.m.
    This function defines the grid and defines the value function.  There 
    are two alternatives.  First, the user could input either the x or Rgrid.
    This should supercede any other option.  Otherwise we use the SS
    computation and use efault DeltaX,DeltaR parameters to set the deviation
    from the SS.

    params.flagSetxGrid sets the flag for either using the default grid (Value = 0)
     or using the user defined grid (Value = 1)

    params.flagsetRdgrid sets the flag for either using the default grid (Value = 0)
    or using the user defined grid (Value = 1)'''
    #???Do we want to name all of the elements in params like params[0] = ___
    #How are we importing it?
    
    #Find the SS
    [xSS, RSS, NA] = findSteadyState(0,3,params)

    #Check the flags and define Grid Endpoints
    if params.flagSetxGrid == 1:
        xMin = params.xMin
        xMax = params.xMax
        disp('Msg: Using user defined grid on x')
        #^Will proceed using user defined gridpoints
    else:
        xMin = xSS - params.DeltaX
        xMax = xSS + params.DeltaX
        disp('Msg: Using default grid around SS')
        #^Will proceed using the default gridpoints
    #Uniformly space the gridpoints
    xGrid = np.linspace(xMin, xMax, params.xGridSize, endpoint = True)
    
    #Update the params struct
    params.xGrid = xGrid
    params.xLL = xMin
    params.xUL = xMax

    #Check the flags and define Grid Endpoints
    if params.flagSetRGrid == 1:
        RMin = params.RMin
        Rmax = params.Rmax
        disp('Proceeding with RGrid based on user inputs')
        #^Will proceed using user defined grid
    else:
        RMin = RSS - params.DeltaR
        RMax = RSS + params.DeltaR
        disp('Proceeding with default RGrid')
        #^Will proceed using default grid

    #Uniformly space the gridpoints
    RGrid = np.linspace(RMin, RMax, params.RGridSize, endpoint = True)
    params.RGrid = RGrid
    #Gives the total gridsize
    GridSize = params.xGridSize * params.RGridSize * params.sSize

    #Update params struct
    param.Gridsize = GridSize
    param.xMin = xMin
    params.xMax = xMax
    params.RMax = RMax
    params.RMin = RMin

    ##Define Functional Space
    #Priorly used the CompEcon Library routine 'fundefn' to create a functional
    #space which stores key settings for the basis polynomials, domand and nodes.
    #V(s) is he functional space for the value function given the discrete shock
    #Need to look up documentation on Compecon toolbox function fundefn

    V = np.zeros(2)
    V[0] = fundefn(___)
    V[1] = V[0]
    #We return the updated params and the funcitonal space
    return params, V

def InitializeCoeff(params, V):
    #This function is the equivalent of InitializeCoeff.m
    ''' INITIALIZE THE COEFF
     This section uses the stationary policies to initialze the value
     functions. The block prdoduces three outcomes 
     1. domain which is the vectorized domain
     2. PolicyRulesStore : This serves as a matrix of guess for optimal
     policies at the interpolation nodes
    3. c0 : initial coeffecients'''
    xGrid = params.xGrid
    RGrid = params.RGrid
    for s_ in range(params.sSize):
        n = 1
        if s==1:
            for xctr in range(params.RGridSize):
                x_ = XGrid[xctr]
                R_ = RGrid[Rctr]
                domain_[s_, n, :] = [x_, R_]
                #Initialize the guess for Stationary Policies
                cRat = R_**(-1. / params.sigma)
                c1_1 = (0.8*(params.n1*params.theta_1+params.n2*params.theta_2)-params.g[1])/(params.n1+cRat*params.n2)
                c1_2 = (0.8*(params.n1*params.theta_1+params.n2*params.theta_2)-params.g[2])/(params.n1+cRat*params.n2)
                c2_1 = cRat*c1_1
                ###Need to finish later.  Finished driving and had to stop.
                #Compute the Stationary Policies using the
                #SteadyStateResiduals Routine
                [xSS, , exitFlag] = opt.fsolve()
                


#BUILD GRID
[params, V] = BuildGrid(params)
disp('Msg: Completed definition of functional space')

#INITIALIZE THE COEFF
disp('Msg: Initializing the Value Function...')
[domain, c, PolicyRulesStore] = InitializeCoeff(params, V)
disp('Msg: ... Completed')
