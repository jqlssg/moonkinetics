///|
test "hand-derived dimer propensity has falling factorial and volume" {
  let net = Network::new(
    ["A", "B"],
    [Reaction::new([2, 0], [0, 1], 3.0)],
    volume=2.0,
  )
  assert_eq(net.propensities([4, 0]), [18.0])
  assert_eq(net.fluxes([2.0, 0.0]), [12.0])
  assert_eq(net.drift([2.0, 0.0]), [-24.0, 12.0])
  assert_eq(net.count_drift([4, 0]), [-36.0, 18.0])
  assert_eq(net.propensities([1, 0]), [0.0])
}

///|
test "zero order scales by volume and null channel can be positive" {
  let net = Network::new(
    ["A"],
    [Reaction::new([0], [1], 2.0), Reaction::new([1], [1], 1.0)],
    volume=3.0,
  )
  assert_eq(net.propensities([4]), [6.0, 4.0])
  assert_eq(net.total_propensity([4]), 10.0)
  assert_eq(net.drift([2.0]), [2.0])
}
