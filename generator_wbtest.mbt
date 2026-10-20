///|
test "CTMC sums parallel channels, ignores null events and rejects open sets" {
  let n = Network::new(["A", "B"], [
    Reaction::new([1, 0], [0, 1], 2.0),
    Reaction::new([1, 0], [0, 1], 3.0),
    Reaction::new([0, 1], [1, 0], 1.0),
    Reaction::new([1, 0], [1, 0], 7.0),
  ])
  let g = n.generator([[1, 0], [0, 1]])
  assert_eq(g.dense(), [[-5.0, 5.0], [1.0, -1.0]])
  assert_eq(g.derivative([1.0, 0.0]), [-5.0, 5.0])
  assert_eq(g.row_edges(0), [(1, 2.0), (1, 3.0)])
  let bad = try {
    let _ = n.generator([[1, 0]])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
  let bad = try {
    let _ = n.generator([[1, 0], [1, 0]])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
