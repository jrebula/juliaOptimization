

# might have to run these two lines if you haven't installed NLOpt yet:
#Pkg.add("NLopt")
#Pkg.update()

module Optimizer

performSafeChecking = true;

using Base.Test
using NLopt


type Optimization
  opt

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



function optimize(o::Optimization, initialState::Trajectory)

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




function dynamics_constraint_function(stateTraj::Trajectory,
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

end



end



