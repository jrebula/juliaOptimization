

module TestTrajectoryOptimization

using TrajectoryOptimization

type TestSystemState <: TrajectoryOptimization.SystemState
  x::Float64
  y::Float64
  xDot::Float64
  yDot::Float64
end

TestSystemState() = TestSystemState(0, 0, 0, 0)

type TestSystemInput <: TrajectoryOptimization.SystemInput
  fX::Float64
  fY::Float64
end
TestSystemInput() = TestSystemInput(0, 0)



end




