///|
test "minimal standard published seed-one vector and open uniform" {
  let r = Random::new(1)
  assert_eq(r.next_int(), 16807)
  assert_eq(r.next_int(), 282475249)
  assert_eq(r.next_int(), 1622650073)
  for _ in 0..<1000 {
    let u = r.uniform()
    assert_true(u > 0.0 && u < 1.0)
  }
  assert_eq(r.categorical([0.0, 5.0, 0.0]), 1)
  assert_true(r.exponential(2.0) > 0.0)
  let bad = try {
    let _ = Random::new(0)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
