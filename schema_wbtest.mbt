///|
test "schema roundtrip, escaped names and exact integer rejection" {
  let net = Network::new(["A", "B,\""], [Reaction::new([1, 0], [0, 1], 2.0)])
  let parsed = network_from_json(net.to_json().stringify())
  assert_eq(parsed.species(), net.species())
  assert_eq(parsed.propensities([3, 0]), [6.0])
  let bad = try {
    let _ = network_from_json("{\"schema\":1,\"schema\":1}")
    false
  } catch {
    _ => true
  }
  assert_true(bad)
  let bad = try {
    let _ = guarded_json("{\"rate\":1,\"r\\u0061te\":2}")
    false
  } catch {
    _ => true
  }
  assert_true(bad)
  let bad = try {
    let _ = coefficients([1.5])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
