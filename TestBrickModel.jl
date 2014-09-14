

##################


if (false)
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
end
