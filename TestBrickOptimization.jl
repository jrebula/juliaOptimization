
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

custom_test_handler(r::Test.Success) = nothing;
custom_test_handler(r::Test.Failure) = error("test failed: $(r.expr)");
custom_test_handler(r::Test.Error)   = rethrow(r);

Test.with_handler(custom_test_handler) do
    include("BrickTrajectoryOptimization.jl")
    #
    state_traj = BrickModel.BrickStateTrajectory(num_states) + 4;
    input_traj = BrickModel.BrickInputTrajectory(num_states);
    input_traj.inputs[1].fX = 10;
    initialState = BrickTrajectoryOptimization.OptimizationState(state_traj, input_traj, dt)
    #
    state_mask_describing_changable_states = 1 + (initialState * 0)
    stateMaskDescribingChangableStates.inputTraj = state_mask_describing_changable_states.inputTraj * 0;
    state_mask_describing_changable_states.dt = 0
end

    @test opt = BrickTrajectoryOptimization.Optimization(state_mask_describing_changable_states)

    (optimal_state, optimal_cost, return_flag) =
        BrickTrajectoryOptimization.optimize(opt, initial_state);



#    BrickTrajectoryOptimization.
#    addEqualityConstraint(opt,
#                          (returnVal, state, gradient) -> )
