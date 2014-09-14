


# might have to run these two lines if you haven't installed NLOpt yet:
#Pkg.add("NLopt")
#Pkg.update()

module BrickTrajectortyOptimization

using Base.Test
using NLopt

#include("BrickModel.jl")
using BrickModel

#include("C:\\Users\\john\\.julia\\NLopt\\src\\NLopt.jl")
#"C:\Users\john\.julia\NLopt\\src"

numStates = 10;
dt = 0.01;

stateTraj = BrickModel.BrickStateTrajectory(numStates) + 4;
inputTraj = BrickModel.BrickInputTrajectory(numStates);
errorTrajectory = deepcopy(stateTraj);

inputTraj.inputs[1].fX = 10;





BrickTrajectoryOptimization.performAnOptimization()


##################
end
