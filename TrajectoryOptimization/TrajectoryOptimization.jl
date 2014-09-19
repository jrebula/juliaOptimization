

module TrajectoryOptimization


abstract SystemState

abstract SystemInput

type Trajectory{S <: SystemState, I <: SystemInput}
  states::Array{S}
  inputs::Array{I}

  #  Trajectory(numElements::Int) = new(Array(S, numElements), Array(I, numElements))
  Trajectory(numElements::Int) = new([S() for i = 1:numElements], [I() for i = 1:numElements])
end




#=
type TestState <: SystemState

end

type TestInput <: SystemInput

end

t = Trajectory{TestState, TestInput}(3)
=#



end




