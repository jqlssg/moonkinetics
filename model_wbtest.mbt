///|
test "model rejects invalid and owns its input arrays" {
  let a = [1, 0]
  let r = Reaction::new(a, [0, 1], 2.0)
  a[0] = 4
  assert_eq(r.reactants(), [1, 0])
  let net = Network::new(["A", "B"], [r])
  assert_eq(net.species_count(), 2)
  assert_eq(net.reaction(0).constant(), 2.0)
  let invalid = try {
    let _ = Network::new(["A", "A"], [])
    false
  } catch {
    _ => true
  }
  assert_true(invalid)
  let invalid = try {
    let _ = Reaction::new([5], [0], 1.0)
    false
  } catch {
    _ => true
  }
  assert_true(invalid)
}
