


# might have to run these two lines if you haven't installed NLOpt yet:
#Pkg.add("NLopt")
#Pkg.update()

module BrickTrajectoryOptimization

performSafeChecking = true;

using Base.Test
using NLopt

import Base.convert

#include("BrickModel.jl")
using BrickModel

#include("C:\\Users\\john\\.julia\\NLopt\\src\\NLopt.jl")
#"C:\Users\john\.julia\NLopt\\src"




type OptimizationState
    state_traj::BrickModel.BrickStateTrajectory
    input_traj::BrickModel.BrickInputTrajectory
    dt::Float64
end

OptimizationState() = OptimizationState(1.0);
OptimizationState(dt::Float64) = OptimizationState(dt, 0);
OptimizationState(dt::Float64, numStates::Int) = OptimizationState(BrickModel.BrickStateTrajectory(numStates),
                                                                   BrickModel.BrickInputTrajectory(numStates),
                                                                   dt);

function toVector(o::OptimizationState)
  [BrickModel.toVector(o.state_traj) BrickModel.toVector(o.input_traj) dt]
end

function fromVector!(o::OptimizationState, vec::Vector)
  n = length(o.stateTraj)
  o.state_traj.fromVector!(vec[1:(n * 4)]);
  o.input_traj.fromVector!(vec[(n * 4) +  (1:2*n)]);
  o.dt = vec[n*6 + 1];

  [BrickModel.toVector(o.state_traj) BrickModel.toVector(o.input_traj) dt]
end


a = [1, 2, 3];

a[1:2]




type Optimization
  opt

  stateMaskDescribingChangableStates::OptimizationState
  vectorMaskDescribingChangableStates::Vector{Bool}

  lowerBounds::Vector
  upperBounds::Vector
end


function Optimization(stateMaskDescribingChangableStates::OptimizationState)
  opt = [];

  vectorMaskDescribingChangableStates = vec(toVector(stateMaskDescribingChangableStates) .> 0);

  lowerBounds = vec(ones(size(vectorMaskDescribingChangableStates)) * -Inf);
  upperBounds = vec(ones(size(vectorMaskDescribingChangableStates)) * Inf);


  return Optimization(opt, stateMaskDescribingChangableStates, vectorMaskDescribingChangableStates,
                      lowerBounds, upperBounds)
end

numStates = 5;
dt = 0.1;
stateMask = OptimizationState(BrickModel.BrickStateTrajectory(numStates), BrickModel.BrickInputTrajectory(numStates), dt)

opt = Optimization(stateMask)






function initializeNLOptProblem(o::Optimization, initialState::OptimizationState)
  opt = NLopt.Opt(NLopt.LN_COBYLA, length(v))

  NLopt.lower_bounds!(o.opt, vec(BrickModel.toVector(stateTraj * -Inf)))
  NLopt.xtol_rel!(o.opt, 1e-4)

  workingStateObjective = deepcopy(initialState);
  NLopt.min_objective!(o.opt,
                             (x, g) -> begin
                               workingStateObjective.fromVector!(x),
                               valueFunction(workingStateObjective, g)
                             end)

  numConstraints = length(BrickModel.toVector(stateTraj));
  #numConstraints += 4

  tolerances = ones(numConstraints) * 1e-4


  workingStateConstraint = deepcopy(initialState);
  NLopt.equality_constraint!(o.opt,
                             (r, x, g) -> begin
                               workingStateConstraint.fromVector!(x),
                               dynamics_constraint_function(workingStateConstraint, r, g)
                             end,
                             tolerances, inputTraj)

end



function optimize(o::Optimization, initialState::OptimizationState)

  if o.opt == None
    initializeNLOptProblem(o, initialState)
  end

  show(vec(BrickModel.toVector(stateTraj)))

  minx = minf = ret = 0
  @time for i in 1:2
    (minf,minx,ret) = NLopt.optimize(opt, vec(BrickModel.toVector(stateTraj)))
  end

  BrickModel.fromVector!(stateTraj, minx)

  println("got $minf at after $numberOfFunctionEvaluations iterations (returned $ret)")

  show(stateTraj)
end






function checkNotNaN(x)
  if any(isnan(x))
    error("some of the given x vector was NaN!")
  end
end


function valueFunction(s::OptimizationState, grad::Vector)
  checkNotNaN(x)
  if length(grad) > 0
      grad[:] = 0;
  end
  global numberOfFunctionEvaluations
  numberOfFunctionEvaluations::Int += 1

  s.


  val = norm(x[1:4])
  return val
end



function dynamics_constraint_function(stateTraj::OptimizationState,
                                      constraints::Vector,
                                      gradient::Matrix)
    if (performSafeChecking)
        checkNotNaN(BrickModel.toVector(stateTraj))
        checkNotNaN(BrickModel.toVector(inputTraj))
    end

    if length(gradient) > 0
        gradient[:,:] = 0;
    end
    xStep = 0.1;

    BrickModel.fromVector!(stateTraj, x)
    BrickModel.calculateSimulatedTrajectoryError!(stateTraj, inputTraj, dt, errorTrajectory)

    v = vec(BrickModel.toVector(errorTrajectory));
    i = 1
    for i in 1:length(result)
        result[i] = v[i];
    end
    #i++;

    # show(typeof(BrickModel.toVector(stateTraj.states[end] - stateTraj.states[1])))
    # show(typeof([0 0 0 0.0]))

    #BrickModel.toVector(stateTraj.states[end] - stateTraj.states[1]) - [0 0 0 0.0]

    #result[i + (1:4)] = vec(BrickModel.toVector(stateTraj.states[end] - stateTraj.states[1])) - [xStep 0 0 0.0]

    #println("finished constraint")
    #x[1]
end



end



