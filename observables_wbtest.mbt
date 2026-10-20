///|
test "path exposure is exact including terminal interval" {
  let n = Network::new(["X"], [Reaction::new([0], [1], 2.0)])
  let p = Path::from_events(n, [1], [0.25, 0.75], [0, 0], 1.0)
  assert_eq(p.integrated_counts(), [2.0])
  assert_eq(n.integrated_propensities(p), [2.0])
  assert_eq(n.event_counts(p), [2])
}

///|
test "bounded CTMC observable returns retained expectation and tail error" {
  let n = Network::new(["A", "B"], [Reaction::new([1, 0], [0, 1], 1.0)])
  let g = n.generator([[1, 0], [0, 1]])
  let d = g.evolve([1.0, 0.0], 1.0)
  let (m, b) = g.mean_counts(d)
  assert_true((m[0] - @math.exp(-1.0)).abs() < b[0] + 1.0e-14)
  let (value, bound) = g.expectation(d, [-2.0, 2.0])
  assert_true((value - (2.0 - 4.0 * @math.exp(-1.0))).abs() < bound + 1.0e-13)
}
