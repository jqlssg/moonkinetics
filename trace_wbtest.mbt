///|
test "trace copies inputs, interpolates and integrates exact linear samples" {
  let values = [[1.0], [3.0], [5.0]]
  let t = Trace::new([0.0, 1.0, 2.0], values)
  values[0][0] = 999.0
  assert_eq(t.at(0.5), [2.0])
  assert_eq(t.last(), [5.0])
  assert_eq(t.integral(), [6.0])
  assert_eq(regular_grid(2.0, 2), [0.0, 1.0, 2.0])
  let bad = try {
    let _ = Trace::new([0.0, 0.0], [[1.0], [2.0]])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
