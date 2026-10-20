///|
test "path likelihood includes unobserved final waiting and rate score" {
  let n = Network::new(["X"], [Reaction::new([0], [1], 2.0)])
  let p = Path::from_events(n, [0], [0.25, 0.75], [0, 0], 2.0)
  assert_true(
    (n.log_likelihood(p) - (2.0 * @math.ln(2.0) - 4.0)).abs() < 1.0e-12,
  )
  assert_eq(n.path_score(p), [-1.0])
  let empty = Path::from_events(n, [0], [], [], 2.0)
  assert_eq(n.log_likelihood(empty), -4.0)
  let h = 1.0e-5
  let score = (
      n.with_rate(0, 2.0 + h).log_likelihood(p) -
      n.with_rate(0, 2.0 - h).log_likelihood(p)
    ) /
    (2.0 * h)
  assert_true((score + 1.0).abs() < 1.0e-8)
}

///|
test "likelihood refuses dimension-compatible wrong stoichiometry" {
  let n = Network::new(["X"], [Reaction::new([0], [1], 1.0)])
  let p = Path::from_events(n, [0], [0.5], [0], 1.0)
  let other = Network::new(["X"], [Reaction::new([0], [2], 1.0)])
  let bad = try {
    let _ = other.log_likelihood(p)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
