///|
fn main raise {
  // Input: a portable versioned model, not a built-in biological special case.
  let text = "{\"schema\":1,\"species\":[\"A\",\"B\"],\"volume\":20,\"reactions\":[{\"input\":[1,0],\"output\":[0,1],\"rate\":2},{\"input\":[0,1],\"output\":[1,0],\"rate\":1}]}"
  let net = @k.network_from_json(text)
  let times = @k.regular_grid(2.0, 20)
  let ode = net.integrate([1.0, 0.0], times, max_step=0.002)
  let states = Array::makei(21, a => [a, 20 - a])
  let generator = net.generator(states)
  let initial = Array::make(21, 0.0)
  initial[20] = 1.0
  let distribution = generator.evolve(initial, 2.0)
  let (mean, error) = generator.mean_counts(distribution)
  let p = 1.0 / 3.0 + 2.0 / 3.0 * @math.exp(-6.0)
  assert_true((mean[0] - 20.0 * p).abs() < error[0] + 1.0e-10)
  assert_true((ode.last()[0] - p).abs() < 1.0e-9)
  let path = net.simulate([20, 0], 2.0, seed=2026)
  for state in path.count_states() {
    assert_eq(state[0] + state[1], 20)
  }
  assert_eq(net.conservation_residual([1.0, 1.0]), [0.0, 0.0])
  let replay = @k.path_from_json(net, path.to_json().stringify())
  assert_eq(replay.count_states(), path.count_states())
  println(
    "biochemical: conserved=20 expected_A=\{mean[0]} ODE_A=\{ode.last()[0]} tail=\{distribution.tail_bound()}",
  )
  println(distribution.to_json().stringify())
}
