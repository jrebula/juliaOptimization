

include("brickModel.jl")
using Base.Test



a = BrickModel.BrickState(1.0);
b = BrickModel.BrickState(2);
c = a+b;
@test c.x == 3

a = BrickModel.BrickInputTrajectory(5);
b = a + 1.0;
c = (b+1) * (b+2);
for i = 1 : 5
  for field in {:fX, :fY}
    @test getfield(c.inputs[i], field) == 6
  end
end

@test typeof(BrickModel.toVector(a)) == Array{Float64, 2}

num_states = 10;
dt = 0.01;

state_traj = BrickModel.BrickStateTrajectory(num_states) + 4;
input_traj = BrickModel.BrickInputTrajectory(num_states);

input_traj = input_traj + 1;
input_traj = 0 * input_traj;

input_traj.inputs[1].fX = 10;

dXdX = rand(4,4)
dXdU = rand(4,2)

BrickModel.simulate!(state_traj.states[1], input_traj.inputs[1], dt,
                     derivativesWithRespectToState = dXdX,
                     derivativesWithRespectToInput = dXdU)

#println("dXdU"), show(dXdU), println("")
#println("dXdX"), show(dXdX), println("")

dXdX = zeros(4*num_states, 4*num_states) * 4
dXdU = zeros(4*num_states, 2*num_states) * 4

BrickModel.simulate!(state_traj, input_traj, dt,
                     derivativesWithRespectToState = dXdX,
                     derivativesWithRespectToInput = dXdU)

#println("dXdU"), show(dXdU), println("")
#println("dXdX"), show(dXdX), println("")

show(state_traj)

v = BrickModel.toVector(state_traj)
BrickModel.fromVector!(state_traj, v)

show(state_traj)
show(v)


