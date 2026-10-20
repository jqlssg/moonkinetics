///|
fn main raise {
  let net = @k.Network::new(
    ["A", "B"],
    [
      @k.Reaction::new([1, 0], [0, 1], 2.0),
      @k.Reaction::new([0, 1], [1, 0], 1.0),
    ],
    volume=20.0,
  )
  let path = net.simulate([20, 0], 2.0, seed=2026)
  for state in path.count_states() {
    assert_eq(state[0] + state[1], 20)
  }
  let replay = @k.path_from_json(net, path.to_json().stringify())
  assert_eq(replay.count_states(), path.count_states())
  println(path.to_json().stringify())
  println(net.log_likelihood(path))
}
