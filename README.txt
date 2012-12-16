This files dicusses the stucture of the code for solving the Golosov-Sargent Economy

The code is subdivided into folders that contain the files for the specific tasks


1. MAIN : This folder contains the main file that currently is setup to compute and simualte 3 economies with sigma =1,2,3. 
	+ RunMainWithAltSigma.m : This file solves the G-S economy with BGP preferences of the form
		psi.c^(1-sigma)/1-sigma+(1-psi)log[1-l] and then simualtes the economy
		
 2. BELLMANEQUATIONCODE : This folder contains the file that organizes the task of running the iterations on the bellman equation from t=1. It also has the code that helps in solving the T=0 problem. Below we discuss the key parts of the code
	+ SetParaStruc : Sets the default values for Para struc, updates the path for libraries and storing data
	
	+ MainBellman.m :  This is the main file computes the value function via time iteration for the parameters passsed in the structure Para. 
	
	+ BuildGrid.m : This function defines the grid  and defines the value function  We have two alternatives First the user could input either the x or Rgrid. this should supersede any other option. Otherwise we use the steady state computation and use the default DeltaX,DeltaR paramters to set the deviation from the steady state. 
	
	+ InitializeCoeff.m : This file uses the stationary policies to initialze the value functions. 
	
	+ CheckGradNAG.m : This performs the inner optimization. Refer below for dependecies and details
	
	+ HandleUnresovledPoints.m :  This file takes care of the points that the inner optimization failed using a homotopy - kind of approach
	
	+ UpdateCoeff.m : This file updates coeff of the value function using the CompEcon routine t fit multidimendional splines
	
3.INNEROPTIMIZATIONCODE : This folder contains the files for doing the inner optimization for a specific point in the state space using a root solver from the NAG library  The key files are 
	
	+ CheckGradNAG.m : The inner optimization is solved in two steps First ggnoring  the bounds on $x$, we use the  FOCs,  solve the unconstrained problem If the implied policies for $x(s)$ in step 1 violate the bounds we use the KKT conditions to solve the optimization

	+ BelObjectiveUncondGradNAGBGP.m : Computes the gradient of the bellman equation objective with respect to c_1(1), c_1(2) and c_2(1).  Substitues out for the rest of the variables using their respective gradients. 
	
	+ resFOCBGP_alt.m  :  Computes the gradient of the bellamn equation objective under the constraints that xLL <= xprime <= xUL. 

4. SIMULATIONCODE : This folder containts the files for computing the simulations. 
	
	+ RunSimulations : Computes the simulations using a user given random seed and initial conditions b_{-1},s0. It first solves the T-0 problem and then uses stored coeff for the value function to solve for the policies and state transitions

5. PLOTTINGCODE : This contains the code for plotting the diagnostics and other graphs 

6.LINEARAPPROXIMATIONCODE : This contains the code for computing the linear approximation around the steady state using only the primitives

7. STEADYSTATECODE: This contains the code for computing the steady state 

	+findSteadyState .m : Computes the steady state for a given Para structure