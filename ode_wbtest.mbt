///|
test "RK4 reversible reaction agrees with analytic concentration and mass" {
  let net = Network::new(["A", "B"], [
    Reaction::new([1, 0], [0, 1], 2.0),
    Reaction::new([0, 1], [1, 0], 1.0),
  ])
  let t = net.integrate([1.0, 0.0], regular_grid(2.0, 20), max_step=0.002)
  let a = 1.0 / 3.0 + 2.0 / 3.0 * @math.exp(-6.0)
  assert_true((t.last()[0] - a).abs() < 1.0e-9)
  for row in t.samples() {
    assert_true((row[0] + row[1] - 1.0).abs() < 1.0e-12)
  }
}

///|
test "RK4 refuses negative stages and exhausted budgets" {
  let net = Network::new(["A"], [Reaction::new([1], [0], 100.0)])
  let bad = try {
    let _ = net.integrate([1.0], [0.0, 1.0], max_step=1.0)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
  let bad = try {
    let _ = net.integrate([1.0], [0.0, 1.0], max_step=0.001, max_steps=1)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
