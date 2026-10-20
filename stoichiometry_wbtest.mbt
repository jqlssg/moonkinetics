///|
test "checked stoichiometry, volume and atomic firing" {
  let n = Network::new(
    ["A", "B"],
    [Reaction::new([2, 0], [0, 1], 3.0)],
    volume=2.0,
  )
  assert_eq(n.stoichiometry(), [[-2], [1]])
  assert_eq(n.orders(), [2])
  let x = [3, 0]
  assert_eq(n.fire(x, 0), [1, 1])
  assert_eq(x, [3, 0])
  assert_true(!n.enabled([1, 0], 0))
  assert_eq(n.concentrations([4, 2]), [2.0, 1.0])
  assert_eq(n.rounded_counts([2.0, 1.0]), [4, 2])
  assert_eq(n.with_rate(0, 4.0).reaction(0).constant(), 4.0)
  let bad = try {
    let _ = n.fire([1, 0], 0)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
