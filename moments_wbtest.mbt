///|
test "online and merged covariance match batch hand calculation" {
  let all = Moments::new(2)
  let left = Moments::new(2)
  let right = Moments::new(2)
  for row in [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]] {
    all.push(row)
  }
  left.push([1.0, 2.0])
  right.push([2.0, 4.0])
  right.push([3.0, 6.0])
  left.merge(right)
  assert_eq(all.average(), [2.0, 4.0])
  assert_eq(all.covariance(), [[1.0, 2.0], [2.0, 4.0]])
  assert_eq(left.covariance(), all.covariance())
  all.merge(all)
  assert_eq(all.sample_count(), 6)
  assert_eq(all.average(), [2.0, 4.0])
  let bad = try {
    all.push([0.0])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
  assert_eq(all.sample_count(), 6)
}
