///|
test "two-state continuous probabilities and segmented semigroup" {
  let n = Network::new(["A", "B"], [
    Reaction::new([1, 0], [0, 1], 2.0),
    Reaction::new([0, 1], [1, 0], 1.0),
  ])
  let g = n.generator([[1, 0], [0, 1]])
  for t in [0.0, 0.1, 1.0, 20.0] {
    let d = g.evolve([1.0, 0.0], t)
    let a = 1.0 / 3.0 + 2.0 / 3.0 * @math.exp(-3.0 * t)
    assert_true((d.masses()[0] - a).abs() < 2.0e-12)
    assert_true(d.tail_bound() <= 1.0e-12)
    assert_true((d.masses()[0] + d.masses()[1] - 1.0).abs() < 2.0e-12)
  }
  assert_true(g.evolve([1.0, 0.0], 20.0).segment_count() > 1)
  let bad = try {
    let _ = g.evolve([1.0, 0.0], 1.0, max_terms=1)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "absorbing and null channels have zero generator and exact identity" {
  let n = Network::new(["X"], [Reaction::new([1], [1], 100.0)])
  let g = n.generator([[1]])
  assert_eq(g.evolve([1.0], 1.0e20).masses(), [1.0])
  assert_eq(g.evolve([1.0], 1.0).tail_bound(), 0.0)
}
