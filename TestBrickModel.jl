


custom_test_handler(r::Test.Success) = nothing;
custom_test_handler(r::Test.Failure) = error("test failed: $(r.expr)");
custom_test_handler(r::Test.Error)   = rethrow(r);

Test.with_handler(custom_test_handler) do



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
