///|
test "analytic Jacobian is valid at zero and agrees with finite difference" {
  let n = Network::new(["A", "B"], [
    Reaction::new([1, 0], [0, 1], 2.0),
    Reaction::new([2, 0], [0, 1], 3.0),
  ])
  assert_eq(n.jacobian([0.0, 0.0]), [[-2.0, 0.0], [2.0, 0.0]])
  let h = 1.0e-5
  let l = n.drift([2.0 - h, 1.0])
  let r = n.drift([2.0 + h, 1.0])
  let j = n.jacobian([2.0, 1.0])
  for i in 0..<2 {
    assert_true(((r[i] - l[i]) / (2.0 * h) - j[i][0]).abs() < 1.0e-7)
  }
  assert_eq(n.diffusion([1, 0]), [[2.0, -2.0], [-2.0, 2.0]])
}
