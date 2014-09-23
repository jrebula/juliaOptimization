
#Pkg.add("Debug")
#Pkg.update()


using NLopt
using Base.Test
using Debug

module TrajectoryOptimization

abstract SystemState
abstract SystemInput

type Trajectory{S <: SystemState, I <: SystemInput}
  dt::Float64
  states::Array{S}
  inputs::Array{I}

  numElements::Int64;
  sizeOfState::Int64;
  sizeOfInput::Int64;
  totalTrajectorySize::Int64;

  Trajectory(dt::Float64, numElements::Int) = begin
    sizeOfState = length(toVector(S()));
    sizeOfInput = length(toVector(I()));
    totalTrajectorySize = 1 + sizeOfState * numElements + sizeOfInput * numElements;

    stateIndeces = [];
    inputIndeces = [];
    lastStateIndex = 1 + traj.numElements * traj.sizeOfState;
    for stateNumber = 1 : traj.numElements
      push!(stateIndeces, 1 + (stateNumber-1) * traj.sizeOfState + (1:traj.sizeOfState));
      push!(inputIndeces, lastStateIndex + (inputNumber-1) * traj.sizeOfInput + (1:traj.sizeOfInput));
    end

    new(dt,
          [S() for i = 1:numElements],
          [I() for i = 1:numElements],
          numElements,
          sizeOfState,
          sizeOfInput,
          totalTrajectorySize);
    end
end


type Optimization
  initialCondition::TrajectoryOptimization.Trajectory
  opt;

  Optimization(initialCondition::TrajectoryOptimization.Trajectory) = begin
    ret = new(initialCondition, NLopt.Opt(NLopt.LN_COBYLA, length(toVector(initialCondition))));
    addEqualityConstraint!(optim::Optimization, dynamics_constraint_function(workingStateConstraint, r, g));
    ret;
  end
end


end





function differentialEquation!(state::TrajectoryOptimization.SystemState,
                               input::TrajectoryOptimization.SystemInput,
                               stateDot::TrajectoryOptimization.SystemState)
  if typeof(state) != typeof(stateDot)
    error("called differential equation with state type of $(typeof(state)), and stateDot type of $(typeof(stateDot))")
  end
  error("must define differentialEquation for types of $(typeof(state)), $(typeof(input)), $(typeof(stateDot))");
end


defineVectorConversionsForTypeWithFloat64Fields = (typeName) -> begin
  @eval toVector(s::$typeName) = begin
    fieldNames = names(s);
    vec = zeros(length(fieldNames));
    toVector!(s, vec)
    vec;
  end

  @eval toVector!(s::$typeName, vec::AbstractArray) = begin
    fieldNames = names(s);
    i = 1;
    for fieldName in fieldNames
      vec[i] = getfield(s, fieldName);
      i += 1;
    end
  end

  @eval fromVector!(s::$typeName, vec::AbstractArray) = begin
    fieldNames = names(s);
    i = 1;
    for fieldName in fieldNames
      setfield!(s, fieldName, vec[i]);
      i += 1;
    end
  end
end

defineFourFunctionMathForTypeWithNumericFields = (typeName) -> begin
  for op = (:+, :*, :-, :/)
    @eval ($op)(a::$typeName, b::$typeName) = begin
      c = deepcopy(a);
      for field in names(a)
        setfield!(c, field, ($op)(getfield(a, field), getfield(b, field)))
      end
      c
    end
  end
end



defineFourFunctionMathForTypeWithNumericFields(:(SystemState));
defineFourFunctionMathForTypeWithNumericFields(:(TrajectoryOptimization.SystemInput));

defineVectorConversionsForTypeWithFloat64Fields(:(TrajectoryOptimization.SystemState));
defineVectorConversionsForTypeWithFloat64Fields(:(TrajectoryOptimization.SystemInput));



function checkRep(traj::TrajectoryOptimization.Trajectory)
  @test traj.numElements == length(traj.states);
  @test traj.numElements == length(traj.inputs);
  if (length(traj.states) > 0)
    @test traj.sizeOfState == length(toVector(traj.states[1]));
    @test traj.sizeOfInput == length(toVector(traj.inputs[1]));
  end
  @test traj.totalTrajectorySize == 1 + (traj.sizeOfState+ traj.sizeOfInput) * traj.numElements
end


function toVector(traj::Trajectory)
  vec = zeros(1 + traj.sizeOfState * traj.numElements + traj.sizeOfInput * traj.numElements);
  toVector!(traj, vec);
end

function toVector!(traj::Trajectory, vec::AbstractArray)
  @test length(vec) == traj.totalTrajectorySize

  vec[1] = traj.dt;
  for stateNumber = 1 : traj.numElements
    stateIndeces = 1 + (stateNumber-1) * traj.sizeOfState + (1:traj.sizeOfState);
    vecStates = sub(vec, stateIndeces);
    toVector!(traj.states[stateNumber], vecStates);
  end
  lastStateIndex = 1 + traj.numElements * traj.sizeOfState;
  for inputNumber = 1 : traj.numElements
    inputIndeces = lastStateIndex + (inputNumber-1) * traj.sizeOfInput + (1:traj.sizeOfInput);
    vecInputs = sub(vec, inputIndeces);
    toVector!(traj.inputs[inputNumber], vecInputs);
  end
  return vec;
end


function fromVector!(traj::Trajectory, vec::AbstractArray)
  @test length(vec) == traj.totalTrajectorySize

  traj.dt = vec[1];
  for stateNumber = 1 : traj.numElements
    stateIndeces = 1 + (stateNumber-1) * traj.sizeOfState + (1:traj.sizeOfState);
    vecStates = sub(vec, stateIndeces);
    fromVector!(traj.states[stateNumber], vecStates);
  end
  lastStateIndex = 1 + traj.numElements * traj.sizeOfState;
  for inputNumber = 1 : traj.numElements
    inputIndeces = lastStateIndex + (inputNumber-1) * traj.sizeOfInput + (1:traj.sizeOfInput);
    vecInputs = sub(vec, inputIndeces);
    fromVector!(traj.inputs[inputNumber], vecInputs);
  end
  return traj;
end






function integrateEachStateInTrajectoryForward!(trajectory::Trajectory)
  stateDot = deepcopy(trajectory.states[1]);
  for i = 1:length(trajectory.states)
     differentialEquation!(trajectory.states[i],
                           trajectory.inputs[i],
                           stateDot);
    integrateStateForward!(trajectory.states[i], stateDot, trajectory.dt)
  end
end

function integrateStateForward!(state::SystemState, stateDot::SystemState, dt)
  fieldNames = names(state);
  for fieldName in fieldNames
    setfield!(state, fieldName,
              getfield(state, fieldName)::Float64 +
              dt * getfield(stateDot, fieldName)::Float64)
  end
end


function calculateDynamicsErrorTrajectory(trajectory::Trajectory)
  ret = deepcopy(trajectory);
  splice!(ret.states, length(ret.states));
  splice!(ret.inputs, length(ret.inputs));
  calculateDynamicsErrorTrajectory!(trajectory, ret);
  return ret;
end

function calculateDynamicsErrorTrajectory!(trajectory::Trajectory, errTraj::Trajectory)
  integratedTrajectory = deepcopy(trajectory);
  integrateEachStateInTrajectoryForward!(integratedTrajectory);
  for i = 1 : (length(trajectory.states)-1)
    errTraj.states[i] = trajectory.states[i + 1] - integratedTrajectory.states[i];
    errTraj.inputs[i] = trajectory.inputs[i];
  end
end




performSafeChecking = true;

function dynamics_constraint_function(stateTraj::Trajectory,
                                      constraints::Vector,
                                      gradient::Matrix)
  if (performSafeChecking)
    checkNotNaN(BrickModel.toVector(stateTraj));
    checkNotNaN(BrickModel.toVector(inputTraj));
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









function setCostFunction!(optim::Optimization, valFunction)
  workingStateObjective = deepcopy(optimization.initialCondition)
  NLopt.min_objective!(optim.opt,
                       (x, g) -> begin
                         fromVector!(workingStateObjective, x),
                         valueFunction(workingStateObjective, g)
                       end);
end

function addEqualityConstraint!(optim::Optimization, constraintFunction)
  # (t::TrajectoryOptimization.Trajectory) -> [t.states[1].x, t.states[1].y])
  workingStateConstraint = constraintFunction(optim.initialTrajectory);
  tolerances = ones(length(toVector(workingStateConstraint))) * 1e-4
  NLopt.equality_constraint!(o.opt,
                             (r, x, g) -> begin
                               fromVector!(workingStateConstraint, x),
                               constraintFunction(workingStateConstraint, r, g)
                             end,
                             tolerances)
end


function optimize(optim::Optimization, initialCondition::Trajectory)
  NLopt.lower_bounds!(o.opt, vec(BrickModel.toVector(stateTraj * -Inf)))
  NLopt.xtol_rel!(o.opt, 1e-4)

  workingStateObjective = deepcopy(initialState);

  numConstraints = length(BrickModel.toVector(stateTraj));
  #numConstraints += 4
end






