


workspace()

immutable TestSystemState
    x::Float64
    y::Float64
    xDot::Float64
    yDot::Float64
end

*(s::TestSystemState, n) = TestSystemState(s.x*n, s.y*n, s.xDot*n, s.yDot*n)

+(a::TestSystemState, b::TestSystemState) =
    TestSystemState(a.x + b.x, a.y + b.y, a.xDot + b.xDot, a.yDot + b.yDot)


immutable TestSystemInput
    fX::Float64
    fY::Float64
end

function getStateDots!(states::Array{TestSystemState},
                      inputs::Array{TestSystemInput},
                      stateDots::Array{TestSystemState})
    for i = 1:length(states)
        state = states[i]
        input = inputs[i]
        stateDots[i] = TestSystemState(state.xDot, state.yDot, input.fX, input.fY)
    end
end




immutable TestSystemStateTrajectory
    states::Array{TestSystemState, 1}
end
TestSystemStateTrajectory(numStates) =
    TestSystemStateTrajectory(Array(TestSystemState,
                                    numStates))


immutable TestSystemInputTrajectory
    inputs::Array{TestSystemInput, 1}
end
TestSystemInputTrajectory(numStates) =
    TestSystemInputTrajectory(Array(TestSystemInput,
                                    numStates))


immutable TestSystemTrajectory
    stateTrajectory::TestSystemStateTrajectory
    inputTrajectory::TestSystemInputTrajectory
end
TestSystemTrajectory(numStates) = 
    TestSystemTrajectory(TestSystemStateTrajectory(numStates),
                         TestSystemInputTrajectory(numStates))


function integrateMultipleShooting!(states::Array{TestSystemState},
                                           inputs::Array{TestSystemInput},
                                           finalStates::Array{TestSystemState},
                                           dt::Float64)
    getStateDots!(states, inputs, finalStates)
    for i = 1:length(states)
        currentState = states[i]
        stateDot = finalStates[i]
        finalStates[i] = currentState + stateDot * dt
    end
    return
end


immutable Integrator
    state::Array{TestSystemState, 1}
    input::Array{TestSystemInput, 1}
    stateDot::Array{TestSystemState, 1}

    Integrator() = new(Array(TestSystemState, 1),
                       Array(TestSystemInput, 1),
                       Array(TestSystemState, 1))
end


function integrateFromInitialConditions!(integrator::Integrator,
                                         initialState::TestSystemState, 
                                         states::Array{TestSystemState},
                                         inputs::Array{TestSystemInput},
                                         dt::Float64)
    states[1] = initialState
    for i = 2:length(states)
        currentInput = inputs[i-1]
        currentState = states[i-1]
        
        integrator.state[1] = currentState
        integrator.input[1] = currentInput
        
        getStateDots!(integrator.state, integrator.input, integrator.stateDot)

        states[i] = currentState + integrator.stateDot[1] * dt
    end
    return
end



#=
    function integrateStatesWithInputsSpecific!(states::Array{TestSystemState},
                                    inputs::Array{TestSystemInput},
                                    finalStates::Array{TestSystemState},
                                    dt::Float64)
    for i = 1:length(states)
        state = states[i]
        input = inputs[i]
        finalStates[i] = TestSystemState(state.x + dt * state.xDot,
                                         state.y + dt * state.yDot,
                                         state.xDot + dt * input.fX,
                                         state.yDot + dt * input.fY)
    end
    return
end
=#



macro timeMany(check, func, numIters)
    :(if $check
        $func();
        @time for i = 1:$numIters
            $func();
        end        
    end)
end








numStates = 1000
dt = 0.001

timeMultipleShootingIntegration = true #false #
timeIntegrationFromInitialConditions = true #false #

initialTrajectory = TestSystemTrajectory(numStates);
for i = 1 : numStates
    initialTrajectory.stateTrajectory.states[i] = TestSystemState(0.0, 0.0, 0.0, 0.0)
    initialTrajectory.inputTrajectory.inputs[i] = TestSystemInput(0.1, 0.2)
end

stateDots = TestSystemStateTrajectory(numStates);
getStateDots!(initialTrajectory.stateTrajectory.states,
              initialTrajectory.inputTrajectory.inputs,
              stateDots.states)

integratedStates = TestSystemStateTrajectory(numStates);

integrateMultipleShooting!() = 
    integrateMultipleShooting!(initialTrajectory.stateTrajectory.states,
                               initialTrajectory.inputTrajectory.inputs,
                               integratedStates.states, dt)

@timeMany(timeMultipleShootingIntegration, integrateMultipleShooting!, 1000)


integratedTrajectory = TestSystemTrajectory(numStates);
integratedTrajectory.inputTrajectory.inputs[:] = initialTrajectory.inputTrajectory.inputs;

integrator = Integrator();
initialCondition = TestSystemState(0.0, 0.0, 0.0, 0.0)
integrateFromInitialConditions!() =
    integrateFromInitialConditions!(integrator,
                                    initialCondition,
                                    integratedTrajectory.stateTrajectory.states,
                                    integratedTrajectory.inputTrajectory.inputs,
                                    dt)

integrateFromInitialConditions!()

#integratedTrajectory.stateTrajectory.states
@timeMany(timeIntegrationFromInitialConditions, integrateFromInitialConditions!, 1000)


#Pkg.add("Winston")
#Pkg.build("Winston")
#using Winston
using PyPlot

allStates = integratedTrajectory.stateTrajectory.states

times = linspace(0, dt * numStates, numStates)

for name = names(allStates[1])
    plot(times, [getfield(allStates[i], name) for i = 1:numStates], label=name)
end

legend()
