///|
test "replayable JSON and safe numeric CSV exports" {
  let n = Network::new(["=unsafe name"], [Reaction::new([0], [1], 2.0)])
  let p = Path::from_events(n, [0], [0.25, 0.75], [0, 0], 1.0)
  let q = path_from_json(n, p.to_json().stringify())
  assert_eq(q.count_states(), p.count_states())
  assert_eq(q.event_times(), p.event_times())
  assert_true(
    p
    .to_csv()
    .has_prefix("event,time,channel,species_index,count\n0,0,-1,0,0\n"),
  )
  assert_true(
    p
    .sample([0.0, 1.0])
    .to_csv()
    .has_prefix("time,species_index,value\n0,0,0\n"),
  )
  match p.sample([0.0, 1.0]).to_json() {
    Object(_) => ()
    _ => fail("expected JSON object")
  }
}
