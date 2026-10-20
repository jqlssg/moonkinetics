///|
fn main raise {
  // S+I -> 2I, I -> R. Volume N gives propensity beta*S*I/N.
  // This is an educational well-mixed process, not a clinical predictor.
  let net = @k.Network::new(
    ["S", "I", "R"],
    [
      @k.Reaction::new([1, 1, 0], [0, 2, 0], 2.0),
      @k.Reaction::new([0, 1, 0], [0, 0, 1], 1.0),
    ],
    volume=100.0,
  )
  let initial = [95, 5, 0]
  let times = @k.regular_grid(8.0, 40)
  let path = net.simulate(initial, 8.0, seed=2026, max_events=200)
  let ode = net.integrate(net.concentrations(initial), times, max_step=0.002)
  let counts = path.count_states()
  for k in 0..<counts.length() {
    assert_eq(counts[k][0] + counts[k][1] + counts[k][2], 100)
    if k > 0 {
      assert_true(counts[k][0] <= counts[k - 1][0])
      assert_true(counts[k][2] >= counts[k - 1][2])
    }
  }
  for row in ode.samples() {
    assert_true((row[0] + row[1] + row[2] - 1.0).abs() < 1.0e-12)
  }
  let events = net.event_counts(path)
  let last = path.last_counts()
  assert_eq(events[0], 95 - last[0])
  assert_eq(events[1], last[2])
  let sampled = path.sample(times)
  assert_eq(sampled.last(), last.map(n => n.to_double()))
  println(
    "epidemic: conserved=100 infections=\{events[0]} recoveries=\{events[1]} final_I=\{last[1]} ODE_I=\{100.0*ode.last()[1]}",
  )
  println(path.to_json().stringify())
}
