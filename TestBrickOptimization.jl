
# might have to run these two lines if you haven't installed NLOpt yet:
#Pkg.add("NLopt")
#Pkg.update()


using Base.Test
using NLopt

include("BrickModel.jl")

#include("C:\\Users\\john\\.julia\\NLopt\\src\\NLopt.jl")
#"C:\Users\john\.julia\NLopt\\src"

numStates = 10;
dt = 0.01;

stateTraj = BrickModel.BrickStateTrajectory(numStates);
inputTraj = BrickModel.BrickInputTrajectory(numStates);

inputTraj.inputs[1].fX = 10;

#show(stateTraj)
#show(inputTraj)

stateTraj = BrickModel.integrateStateTrajectory(stateTraj, inputTraj, dt)


e = BrickModel.calculateSimulatedTrajectoryError(stateTraj, inputTraj, dt)

#show(stateTraj)
#show(e)

A = rand(8,8)
B = rand(8,1)

dXdX = rand(4,4)
dXdU = rand(4,2)

BrickModel.simulate!(stateTraj.states[1], inputTraj.inputs[1], dt,
                     derivativesWithRespectToState = dXdX,
                     derivativesWithRespectToInput = dXdU)

#println("dXdU"), show(dXdU), println("")
#println("dXdX"), show(dXdX), println("")

dXdX = zeros(4*numStates, 4*numStates) * 4
dXdU = zeros(4*numStates, 2*numStates) * 4

BrickModel.simulate!(stateTraj, inputTraj, dt,
                     derivativesWithRespectToState = dXdX,
                     derivativesWithRespectToInput = dXdU)

#println("dXdU"), show(dXdU), println("")
#println("dXdX"), show(dXdX), println("")

v = BrickModel.toVector(stateTraj)

#show(stateTraj)
#show(v)

BrickModel.fromVector!(stateTraj, v)



#if (true) #false) #

println("testing out NLopt!!")

numberOfFunctionEvaluations = 0

trajectory = BrickModel.BrickStateTrajectory(numStates) + 4
inputTraj = BrickModel.BrickInputTrajectory(numStates)

inputTraj.inputs[1].fX = 1

function checkNotNaN(x)
  if any(isnan(x))
    error("some of the given x vector was NaN!")
  end
end


function myfunc(x::Vector, grad::Vector)
  #println("started value function")
  checkNotNaN(x)
  if length(grad) > 0
      grad[:] = 0;
  end
  global numberOfFunctionEvaluations
  numberOfFunctionEvaluations::Int += 1
  #println("f_$numberOfFunctionEvaluations($x)")

  val = norm(x[1:4])
  #println("finished value function $(val)")
  val
end

function myconstraint(result::Vector, x::Vector, grad::Matrix)
  checkNotNaN(x)
  if length(grad) > 0
      grad[:,:] = 0;
  end

  BrickModel.fromVector!(trajectory, x)
  e = BrickModel.calculateSimulatedTrajectoryError(trajectory, inputTraj, dt)

  v = vec(BrickModel.toVector(e));
  for i in 1:length(result)
    result[i] = v[i]
  end

  #println("finished constraint")
  #x[1]
end


opt = NLopt.Opt(NLopt.LN_COBYLA, length(BrickModel.toVector(trajectory)))
#opt = NLopt.Opt(NLopt.LD_SLSQP, length(BrickModel.toVector(trajectory)))

NLopt.lower_bounds!(opt, vec(BrickModel.toVector(trajectory * -Inf)))
NLopt.xtol_rel!(opt,1e-4)

NLopt.min_objective!(opt, myfunc)


tolerances = ones(length(BrickModel.toVector(trajectory))) * 1e-4

NLopt.equality_constraint!(opt, (r, x, g) -> myconstraint(r, x, g), tolerances)


@time for i in 1:10
  (minf,minx,ret) = NLopt.optimize(opt, vec(BrickModel.toVector(trajectory)))
end

BrickModel.fromVector!(trajectory, minx)

println("got $minf at $minx after $numberOfFunctionEvaluations iterations (returned $ret)")

show(trajectory)



