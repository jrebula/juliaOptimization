

abstract OptimizationState


function getVector(s::OptimizationState)
  error("must implement OptimizationState.getVector(s::$(typeof(s)))")
end

function fromVector!(s::OptimizationState, v::Array{Float64, 2})
  error("must implement OptimizationState.fromVector!(s::$(typeof(s))), v::Array{Float64, 2}")
end


#function fillVector!(s::OptimizationState, v::Vector)
#  error("must implement OptimizationState.fillVector(s::$(typeof(s)), v::Vector)!")
#end


#=
type TestOptimizationState <: OptimizationState
  n
end


t = TestOptimizationState(5)

getVector(t)

=#
