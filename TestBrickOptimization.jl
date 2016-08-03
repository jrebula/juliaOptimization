
# might have to run these two lines if you haven't installed NLOpt yet:
#Pkg.add("NLopt")
#Pkg.update()

using Base.Test
using NLopt
using BrickModel
using Base.Test

import Base.convert

#include("C:\\Users\\john\\.julia\\NLopt\\src\\NLopt.jl")
#"C:\Users\john\.julia\NLopt\\src"

num_states = 10;
dt = 0.01;

include("BrickTrajectoryOptimization.jl")

state_traj = BrickModel.BrickStateTrajectory(num_states) + 4;
input_traj = BrickModel.BrickInputTrajectory(num_states);
input_traj.inputs[1].fX = 10;

initial_state = BrickTrajectoryOptimization.OptimizationState(state_traj, input_traj, dt)

BrickTrajectoryOptimization.toVector(initial_state);

state_mask_describing_changable_states = BrickTrajectoryOptimization.OptimizationState(1+(state_traj * 0), input_traj * 0, 0)

opt = BrickTrajectoryOptimization.Optimization(state_mask_describing_changable_states);

@test typeof(opt) == BrickTrajectoryOptimization.Optimization

(optimal_state, optimal_cost, return_flag) =
  BrickTrajectoryOptimization.optimize(opt, initial_state);



#    BrickTrajectoryOptimization.
#    addEqualityConstraint(opt,
#                          (returnVal, state, gradient) -> )
