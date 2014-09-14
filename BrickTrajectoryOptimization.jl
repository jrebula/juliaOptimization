


# might have to run these two lines if you haven't installed NLOpt yet:
#Pkg.add("NLopt")
#Pkg.update()

module BrickTrajectoryOptimization

performSafeChecking = true;

using Base.Test
using NLopt

#include("BrickModel.jl")
using BrickModel

#include("C:\\Users\\john\\.julia\\NLopt\\src\\NLopt.jl")
#"C:\Users\john\.julia\NLopt\\src"


function checkNotNaN(x)
  if any(isnan(x))
    error("some of the given x vector was NaN!")
  end
end


function valueFunction(x::Vector, grad::Vector)
  checkNotNaN(x)
  if length(grad) > 0
      grad[:] = 0;
  end
  global numberOfFunctionEvaluations
  numberOfFunctionEvaluations::Int += 1

  val = norm(x[1:4])
  return val
end

function constraintFunction(result::Vector, x::Vector, grad::Matrix)
  checkNotNaN(x)
  if length(grad) > 0
      grad[:,:] = 0;
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


function constraintFunction(stateTraj::BrickModel.BrickStateTrajectory, inputTraj::BrickModel.BrickInputTrajectory,
                            constraints::Vector, gradient::Matrix)
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



function performAnOptimization(stateTraj::BrickModel.BrickStateTrajectory, inputTraj::BrickModel.BrickInputTrajectory)

  opt = NLopt.Opt(NLopt.LN_COBYLA, length(BrickModel.toVector(stateTraj)))
  #opt = NLopt.Opt(NLopt.LD_SLSQP, length(BrickModel.toVector(trajectory)))

  NLopt.lower_bounds!(opt, vec(BrickModel.toVector(stateTraj * -Inf)))
  NLopt.xtol_rel!(opt,1e-4)

  NLopt.min_objective!(opt, myfunc)

  numConstraints = length(BrickModel.toVector(stateTraj));
  #numConstraints += 4

  tolerances = ones(numConstraints) * 1e-4
  show(tolerances)

  NLopt.equality_constraint!(opt, (r, x, g) -> myconstraint(r, x, g), tolerances, inputTraj)


  show(vec(BrickModel.toVector(stateTraj)))

  minx = minf = ret = 0
  @time for i in 1:2
    (minf,minx,ret) = NLopt.optimize(opt, vec(BrickModel.toVector(stateTraj)))
  end

  BrickModel.fromVector!(stateTraj, minx)

  println("got $minf at after $numberOfFunctionEvaluations iterations (returned $ret)")

  show(stateTraj)

  end

end



