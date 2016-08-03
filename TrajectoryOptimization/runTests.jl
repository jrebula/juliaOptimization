
push!(LOAD_PATH, "C:\\Users\\john\\workspace\\juliaOptimization\\juliaOptimization\\TrajectoryOptimization");

using Base.Test

#include("TrajectoryOptimization.jl");

using TrajectoryOptimization

include("TestTrajectoryOptimization.jl");


a = TestTrajectoryOptimization.TestSystemState();
typeof(a) <: TrajectoryOptimization.SystemState



numPoints = 4;
dt = 0.1;
gravity = -1;
initialTrajectory = TrajectoryOptimization.Trajectory{
  TestTrajectoryOptimization.TestSystemState,
  TestTrajectoryOptimization.TestSystemInput}(dt, numPoints)

#show(initialTrajectory)

for i in 1:numPoints
  initialTrajectory.states[i].x = (i / numPoints);
end

function TrajectoryOptimization.differentialEquation!(
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
for i = 1:length(integratedTrajectory.states)
  @test integratedTrajectory.states[i].yDot == gravity * dt
end


errTraj1 = TrajectoryOptimization.calculateDynamicsErrorTrajectory(initialTrajectory);
errTraj2 = deepcopy(errTraj1);
for i = 1:length(errTraj1.states) errTraj2.states[i].x = i*30.0; end
TrajectoryOptimization.calculateDynamicsErrorTrajectory!(initialTrajectory, errTraj2);


for i = 1:length(errTraj1.states)
  #show(TrajectoryOptimization.toVector(errTraj1.states[i]))
  @test TrajectoryOptimization.toVector(errTraj1.states[i]) == TrajectoryOptimization.toVector(errTraj2.states[i])
end


println();
show(errTraj1);
println();


@test initialTrajectory.numElements == numPoints;
@test initialTrajectory.sizeOfState == 4;
@test initialTrajectory.sizeOfInput == 2;

TrajectoryOptimization.checkRep(initialTrajectory);


vec = zeros(1 + 4*numPoints + 2 * numPoints);
TrajectoryOptimization.toVector!(initialTrajectory, vec);
vec2 = TrajectoryOptimization.toVector(initialTrajectory);
@test vec == vec2


initialTrajectory2 = TrajectoryOptimization.Trajectory{
  TestTrajectoryOptimization.TestSystemState,
  TestTrajectoryOptimization.TestSystemInput}(dt, numPoints)




TestTrajectoryOptimization.fromVector!(initialTrajectory2, vec);

@test initialTrajectory2.dt == initialTrajectory.dt
for i = 1 : numPoints
  @test initialTrajectory2.states[i].x == initialTrajectory.states[i].x
  @test initialTrajectory2.states[i].y == initialTrajectory.states[i].y
  @test initialTrajectory2.states[i].xDot == initialTrajectory.states[i].xDot
  @test initialTrajectory2.states[i].yDot == initialTrajectory.states[i].yDot

  @test initialTrajectory2.inputs[i].fX == initialTrajectory.inputs[i].fX
  @test initialTrajectory2.inputs[i].fY == initialTrajectory.inputs[i].fY
end


optim = TrajectoryOptimization.Optimization(initialTrajectory);



#TrajectoryOptimization.setValueFunction!(optim, (t::TrajectoryOptimization.Trajectory) -> sum([(input.fX^2 + input.fY^2) for input in t.inputs])
#TrajectoryOptimization.addEqualityConstraint!(optim, (t::TrajectoryOptimization.Trajectory) -> [t.states[1].x, t.states[1].y])



