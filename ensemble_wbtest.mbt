///|
test "seeded immigration ensemble matches Poisson mean and variance" {
  let n = Network::new(["X"], [Reaction::new([0], [1], 4.0)])
  let e = n.ensemble([0], [0.0, 1.0], 2000, max_events=100)
  let mean = e.means().last()[0]
  let variance = e.covariances()[1][0][0]
  assert_true((mean - 4.0).abs() < 0.25)
  assert_true((variance - 4.0).abs() < 0.5)
  assert_eq(e.runs(), 2000)
  assert_true(e.standard_errors()[1][0] > 0.0)
  let again = n.ensemble([0], [0.0, 1.0], 2000, max_events=100)
  assert_eq(again.means().samples(), e.means().samples())
}
