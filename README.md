# coastal-management
Identifies adaptation strategies consisting of beach nourishment and managed retreat using discrete dynamic programming

Reference for accompnanying paper:
Cutler, E.M., Albert, M.R., White, K.D., 2019. "A Dynamic Programming Approach to Coastal Adaptation."

The model can be run from the file main.m which takes a parameter struct as an optional input. An appropriate struct can be generated using parameters.m. Outputs include an optimal action vector and corresponding optimal states and value function. The script tests.m can be used to run sensitivity analysis on selected input parameters

The model is solved as a deterministic system, but optimal strategies can be simulated with stochastic storm events, similar to the model in the repository emcutler/ShorelineDynamics, using coupledSystem.m

Requires the MDPSolve toolbox: https://github.com/PaulFackler/MDPSolve
