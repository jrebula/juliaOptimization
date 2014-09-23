
module BrickModel

import Base.show
import Base.convert
using Base.Test

#export getindex


defineFourFunctionMathForTypeWithNumericFields = (typeName) -> begin
  for op = (:+, :*, :-, :/)
    @eval ($op)(a::$typeName, b::$typeName) = begin
      c = deepcopy(a);
      for field in names(a)
        setfield!(c, field, ($op)(getfield(a, field), getfield(b, field)))
      end
      c
    end
    @eval ($op)(a::$typeName, b::Number) = ($op)(a, convert($typeName, b))
    @eval ($op)(a::Number, b::$typeName) = ($op)(convert($typeName, a), b)
  end
end

#  @eval $typeName(a::Number) = begin
#    c = $typeName();
#    for field in names(a)
#      setfield!(c, field, a)
#    end
#    c
#  end



type BrickState
  x::Real
  y::Real
  xDot::Real
  yDot::Real
end

BrickState() = BrickState(0.0, 0.0, 0.0, 0.0);
BrickState(x::Number) = BrickState(x,x,x,x)

function show(io::IO, b::BrickState)
  print(io, "pos = ($(b.x), $(b.y)), vel = ($(b.xDot), $(b.yDot))")
end

function toVector(state::BrickState)
  [state.x state.y state.xDot state.yDot]
end

function fromVector!(state::BrickState, vec::Array)
  (state.x, state.y, state.xDot, state.yDot) = vec
end

Base.convert(::Type{BrickState}, x::Number) = BrickState(x)


#=
for op = (:+, :*, :-, :/)
  @eval ($op)(a::BrickState, b::BrickState) = begin
    c = deepcopy(a);
    for field in names(a)
      setfield!(c, field, ($op)(getfield(a, field), getfield(b, field)))
    end
    c
  end
  @eval ($op)(a::BrickState, b::Number) = ($op)(a, convert(BrickState, b))
  @eval ($op)(a::Number, b::BrickState) = ($op)(convert(BrickState, a), b)
end
=#

defineFourFunctionMathForTypeWithNumericFields(:BrickState)


#a = BrickState(1, 2, 3, 4)
#b = a*3


#promote_rule(::Type{BrickState}, ::Type{Float64} ) = BrickState
#promote_type(Number, BrickState)
#promote(4.0, a)
#convert(BrickState, 4)



type BrickStateTrajectory
  states::Array{BrickState, 1}
end

BrickStateTrajectory(numberOfStates::Int) = BrickStateTrajectory(numberOfStates, BrickState())

BrickStateTrajectory(numberOfStates::Int, state::BrickState) = begin
  ret = BrickState[]
  for i = 1:numberOfStates
    push!(ret, deepcopy(state))
  end
  BrickStateTrajectory(ret)
end

function show(io::IO, b::BrickStateTrajectory)
  println(io, "BrickStateTrajectory ($(length(b.states))): ")
  for el in b.states
    print(io, "  ");
    show(io, el)
    println(io, "");
  end
end

function toVector(traj::BrickStateTrajectory)
  ret = zeros(1,0) #[];
  for i in 1:length(traj.states)
    ret = [ret toVector(traj.states[i])]
  end
  ret
end

function fromVector!(traj::BrickStateTrajectory, vec)
  inds = 1:4
  for i in 1:length(traj.states)
    fromVector!(traj.states[i], vec[inds + 4*(i-1)])
  end
  traj
end



for op = (:+, :*, :-, :/)
  @eval ($op)(a::BrickStateTrajectory, b::BrickStateTrajectory) = begin
    c = BrickStateTrajectory(0);
    for i in 1:length(a.states)
      push!(c.states, ($op)(a.states[i], b.states[i]))
    end
    c
  end

  @eval ($op)(a::BrickStateTrajectory, b::Number) = ($op)(a, convert(BrickStateTrajectory, b, length(a.states)))
  @eval ($op)(a::Number, b::BrickStateTrajectory) = ($op)(convert(BrickStateTrajectory, a, length(b.states)), b)

  @eval ($op)(a::BrickStateTrajectory, b::BrickState) = ($op)(a, convert(BrickStateTrajectory, b, length(a.states)))
  @eval ($op)(a::BrickState, b::BrickStateTrajectory) = ($op)(convert(BrickStateTrajectory, a, length(b.states)), b)
end

Base.convert(::Type{BrickStateTrajectory}, x::Number, numElements::Int) = convert(BrickStateTrajectory, convert(BrickState, x), numElements)
Base.convert(::Type{BrickStateTrajectory}, x::BrickState, numElements::Int) = BrickStateTrajectory(numElements, x)


macro forField(fieldName, a, expressions)
  quote
    for $fieldName in names(a)
      $expressions
    end
  end
end




type BrickInput
  fX::Real
  fY::Real
end

BrickInput(n::Number) = BrickInput(n, n);
BrickInput() = BrickInput(0.0);
Base.convert(::Type{BrickInput}, x::Number) = BrickInput(x)

defineFourFunctionMathForTypeWithNumericFields(:BrickInput)

function toVector(input::BrickInput)
  [input.fX input.fY]
end

function fromVector!(input::BrickInput, vec::Array)
  (input.fX, input.fY) = vec
end

function show(io::IO, b::BrickInput)
  print(io, "fX = $(b.fX), fY = $(b.fY)")
end






type BrickInputTrajectory
  inputs::Array{BrickInput, 1}
end

BrickInputTrajectory(numberOfInputs::Int) = BrickInputTrajectory([BrickInput() for i=1:numberOfInputs])
BrickInputTrajectory() = BrickInputTrajectory(0)

BrickInputTrajectory(numberOfInputs::Int, state::BrickInput) = begin
  ret = BrickInput[]
  for i = 1:numberOfInputs
    push!(ret, deepcopy(state))
  end
  BrickInputTrajectory(ret)
end

Base.convert(::Type{BrickInputTrajectory}, x::Number, numElements::Int) = convert(BrickInputTrajectory, convert(BrickInput, x), numElements)
Base.convert(::Type{BrickInputTrajectory}, x::BrickInput, numElements::Int) = BrickInputTrajectory(numElements, x);


for op = (:+, :*, :-, :/)
  @eval ($op)(a::BrickInputTrajectory, b::BrickInputTrajectory) = begin
    c = BrickInputTrajectory(0);
    for i in 1:length(a.inputs)
      push!(c.inputs, ($op)(a.inputs[i], b.inputs[i]))
    end
    c
  end

  @eval ($op)(a::BrickInputTrajectory, b::Number) = ($op)(a, convert(BrickInputTrajectory, b, length(a.inputs)))
  @eval ($op)(a::Number, b::BrickInputTrajectory) = ($op)(convert(BrickInputTrajectory, a, length(b.inputs)), b)

  @eval ($op)(a::BrickInputTrajectory, b::BrickInput) = ($op)(a, convert(BrickInputTrajectory, b, length(a.inputs)))
  @eval ($op)(a::BrickInput, b::BrickInputTrajectory) = ($op)(convert(BrickInputTrajectory, a, length(b.inputs)), b)
end

Base.convert(::Type{BrickInputTrajectory}, x::Number, numElements::Int) = convert(BrickInputTrajectory, convert(BrickInput, x), numElements)
Base.convert(::Type{BrickInputTrajectory}, x::BrickInput, numElements::Int) = BrickInputTrajectory(numElements, x)


function toVector(traj::BrickInputTrajectory)
  ret = zeros(1,0)
  for i in 1:length(traj.inputs)
    ret = [ret toVector(traj.inputs[i])]
  end
  ret
end

function fromVector!(traj::BrickInputTrajectory, vec)
  inds = 1:2
  for i in 1:length(traj.inputs)
    fromVector!(traj.inputs[i], vec[inds + 2*(i-1)])
  end
  traj
end








function show(io::IO, b::BrickInputTrajectory)
  #println(io, "BrickInputTrajectory: $( map(x -> string(x, " "), b.inputs) )")
  #println(io, "BrickInputTrajectory: $( b.inputs )")
  #println(io, "BrickInputTrajectory: $( b.inputs[:] )")
  #println(io, "BrickInputTrajectory: $( b.inputs... )")
  #println(io, string("BrickInputTrajectory: [", map(x -> string(x, ", "), b.inputs)... , "]"))
  #println(io, "BrickInputTrajectory: $( b.inputs... )")
  println(io, "BrickInputTrajectory ($(length(b.inputs))): ")
  for el in b.inputs
    print(io, "  ");
    show(io, el)
    println(io, "");
  end
end







function simulate!(
    state::BrickState,
    input::BrickInput,
    dt::Real;
    derivativesWithRespectToState=None,
    derivativesWithRespectToInput=None)

  state.x += state.xDot * dt;
  state.y += state.yDot * dt;
  state.xDot += input.fX * dt;
  state.yDot += input.fY * dt;

  conformDerivativeMatrixArguments(4, 2, derivativesWithRespectToState, derivativesWithRespectToInput)

  if (derivativesWithRespectToState != None)
    for i in 1:4
      derivativesWithRespectToState[i,i] = 1
    end
    derivativesWithRespectToState[1,3] = 1
    derivativesWithRespectToState[2,4] = 1
  end

  if (derivativesWithRespectToInput != None)
    derivativesWithRespectToInput[3,1] = 1
    derivativesWithRespectToInput[4,2] = 1
  end

  state
end


function simulate!(
    stateTraj::BrickStateTrajectory,
    inputTraj::BrickInputTrajectory,
    dt::Real;
    derivativesWithRespectToState=None,
    derivativesWithRespectToInput=None
    )
  if (length(stateTraj.states) != length(inputTraj.inputs))
    error("state trajectory and input trajectory must be of the same length")
  end

  numStates = length(stateTraj.states)
  n = 4 * numStates
  m = 2 * numStates

  conformDerivativeMatrixArguments(n, m, derivativesWithRespectToState, derivativesWithRespectToInput)

  derivativesWithRespectToInputSlice = None
  derivativesWithRespectToStateSlice = None

  for i in 1 : numStates
    stateInds = 4*(i-1) + (1:4)
    if (derivativesWithRespectToState != None)
      derivativesWithRespectToStateSlice = sub(derivativesWithRespectToState, stateInds, stateInds);
    end
    if (derivativesWithRespectToInput != None)
      inputInds = 2*(i-1) + (1:2)
      derivativesWithRespectToInputSlice = sub(derivativesWithRespectToInput, stateInds, inputInds);
    end

    BrickModel.simulate!(stateTraj.states[i], inputTraj.inputs[i], dt,
                         derivativesWithRespectToState = derivativesWithRespectToStateSlice,
                         derivativesWithRespectToInput = derivativesWithRespectToInputSlice)
  end
end



function conformDerivativeMatrixArguments(n, m, derivativesWithRespectToState, derivativesWithRespectToInput)
  if (derivativesWithRespectToState != None)
    @assert size(derivativesWithRespectToState) == (n,n) "derivativesWithRespectToState must be either None or $(n)x$(n)!"
    derivativesWithRespectToState[:,:] = 0
  end
  if (derivativesWithRespectToInput != None)
    @assert size(derivativesWithRespectToInput) == (n,m) "derivativesWithRespectToInput must be either None or $(n)x$(m)!"
    derivativesWithRespectToInput[:,:] = 0
  end
end





function integrateStateTrajectory(
    stateTraj::BrickStateTrajectory,
    inputTraj::BrickInputTrajectory,
    dt::Real;
    derivatives=None)

  integratedStates = deepcopy(stateTraj)
  simulate!(integratedStates, inputTraj, dt);

  deleteat!(integratedStates.states, length(integratedStates.states))
  integratedStates.states = [deepcopy(stateTraj.states[1]), integratedStates.states]

  #prepend!(integratedStates.states, deepcopy(stateTraj.states[1]))
  #pop!(integratedStates.states)
  integratedStates
end


function calculateSimulatedTrajectoryError(
    stateTraj::BrickStateTrajectory,
    inputTraj::BrickInputTrajectory,
    dt::Real)
  integratedStates = integrateStateTrajectory(stateTraj, inputTraj, dt);
  errorTrajectory = integratedStates - stateTraj
  errorTrajectory
end

function calculateSimulatedTrajectoryError!(
    stateTraj::BrickStateTrajectory,
    inputTraj::BrickInputTrajectory,
    dt::Real,
    errorTraj::BrickStateTrajectory)
  integratedStates = integrateStateTrajectory(stateTraj, inputTraj, dt);
  for i = 1:length(stateTraj.states), fieldName = [:x, :y, :xDot, :yDot]
    setfield!(errorTraj.states[i],
              fieldName,
              getfield(integratedStates.states[i], fieldName) -
                getfield(stateTraj.states[i], fieldName))
  end
end

end



