


push!(LOAD_PATH, "C:\\Users\\john\\workspace\\juliaOptimization\\juliaOptimization\\TrajectoryOptimization");


#include("TrajectoryOptimization.jl");

using TrajectoryOptimization

include("TestTrajectoryOptimization.jl");


a = TestTrajectoryOptimization.TestSystemState();
typeof(a) <: TrajectoryOptimization.SystemState



numPoints = 4;
dt = 0.1;
initialTrajectory = TrajectoryOptimization.Trajectory{
  TestTrajectoryOptimization.TestSystemState,
  TestTrajectoryOptimization.TestSystemInput}(dt, numPoints)

show(initialTrajectory)

for i in 1:numPoints
  initialTrajectory.states[i].x = (i / numPoints);
end

function differentialEquation!(
    state::TestTrajectoryOptimization.TestSystemState,
    input::TestTrajectoryOptimization.TestSystemInput,
    stateDot::TestTrajectoryOptimization.TestSystemState)
  stateDot.x = state.xDot;
  stateDot.y = state.yDot;
  stateDot.xDot = input.fX;
  stateDot.yDot = input.fY - 1;
end

integratedTrajectory = deepcopy(initialTrajectory);
TrajectoryOptimization.integrateEachStateInTrajectoryForward!(integratedTrajectory);
errTraj = TrajectoryOptimization.calculateDynamicsErrorTrajectory(initialTrajectory);
TrajectoryOptimization.calculateDynamicsErrorTrajectory!(initialTrajectory, errTraj);


optim = TrajectoryOptimization.Optimization(initialTrajectory);

TrajectoryOptimization.setValueFunction!(optim, (t::TrajectoryOptimization.Trajectory) -> sum([(input.fX^2 + input.fY^2) for input in t.inputs])

TrajectoryOptimization.addEqualityConstraint!(optim, (t::TrajectoryOptimization.Trajectory) -> [t.states[1].x, t.states[1].y])



