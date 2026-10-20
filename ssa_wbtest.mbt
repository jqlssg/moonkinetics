///|
test "SSA exact seed, absorption and step rather than linear sampling" {
  let n = Network::new(["A"], [Reaction::new([1], [0], 1.0)])
  let path = n.simulate([1], 100.0)
  assert_eq(path.event_channels(), [0])
  assert_eq(path.last_counts(), [0])
  assert_true(
    (path.event_times()[0] + @math.ln(16807.0 / 2147483647.0)).abs() < 1.0e-12,
  )
  assert_eq(path.sample([0.0, path.event_times()[0], 100.0]).samples(), [
    [1.0],
    [0.0],
    [0.0],
  ])
  assert_eq(n.simulate([0], 1.0, max_events=0).last_counts(), [0])
}

///|
test "SSA budget failure is explicit, replay rejects impossible events" {
  let n = Network::new(["A"], [Reaction::new([0], [1], 100.0)])
  let bad = try {
    let _ = n.simulate([0], 100.0, max_events=1)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
  let bad = try {
    let _ = Path::from_events(n, [0], [0.5, 0.5], [0, 0], 1.0)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
  let p = Path::from_events(n, [0], [0.5], [0], 1.0)
  assert_eq(p.last_counts(), [1])
}
