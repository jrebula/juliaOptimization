


#include("TrajectoryOptimization.jl");

include("TestTrajectoryOptimization.jl");


a = TestTrajectoryOptimization.TestSystemState();

t = typeof(a)

show(t)
parent(t)


typeof(a) <: TrajectoryOptimization.SystemState


numPoints = 10;
initialTrajectory = TrajectoryOptimization.Trajectory{
  TestTrajectoryOptimization.TestSystemState,
  TestTrajectoryOptimization.TestSystemInput}(numPoints)

for i in 1:numPoints
  initialTrajectory.states[i].x = (i / numPoints);
end

function differentialEquation(state::TestSystemState, input::TestSystemInput)

end







