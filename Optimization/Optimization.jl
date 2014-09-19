


module Optimization

include("OptimizationState.jl")

type OptimizationProblem
  s::OptimizationState
end






end


module TestOptimization

using Optimization

type TestState <: Optimization.OptimizationState
  a
  b
  c
  d
end


s = TestState(1,2,3,4)

opt = Optimization.OptimizationProblem(s)

Optimization.setValueFunction(opt, (s::OptimizationState) -> ((s.a^2 + s.b^2) - (s.c^2 + s.d^2)))
Optimization.addEqualityConstraint(opt, (s::OptimizationState) -> s.a - 3)
Optimization.addInequalityConstraint(opt, (s::OptimizationState) -> s.b)

optimalState = opt.optimize(s)

show(optimalState)

end

