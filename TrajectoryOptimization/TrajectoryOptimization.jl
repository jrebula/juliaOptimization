

module TrajectoryOptimization


abstract SystemState
abstract SystemInput

type Trajectory{S <: SystemState, I <: SystemInput}
  dt::Float64
  states::Array{S}
  inputs::Array{I}

  Trajectory(dt::Float64, numElements::Int) =
    new(dt, [S() for i = 1:numElements], [I() for i = 1:numElements])
end





end




