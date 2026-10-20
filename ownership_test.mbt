///|
test "public API defensive copies protect all state boundaries" {
  let inputs = [1, 0]
  let names = ["A", "B"]
  let r = @moonkinetics.Reaction::new(inputs, [0, 1], 1.0)
  let n = @moonkinetics.Network::new(names, [r])
  inputs[0] = 0
  names[0] = "changed"
  let species = n.species()
  species[0] = "changed again"
  assert_eq(n.species(), ["A", "B"])
  let reactants = n.reaction(0).reactants()
  reactants[0] = 0
  assert_eq(n.propensities([2, 0]), [2.0])
  let state = [2, 0]
  let p = @moonkinetics.Path::from_events(n, state, [0.5], [0], 1.0)
  state[0] = 0
  let states = p.count_states()
  states[0][0] = 99
  assert_eq(p.count_states(), [[2, 0], [1, 1]])
  let grid = [0.0, 1.0]
  let values = [[2.0, 0.0], [1.0, 1.0]]
  let trace = @moonkinetics.Trace::new(grid, values)
  grid[0] = 0.2
  values[0][0] = 99.0
  let samples = trace.samples()
  samples[0][0] = 200.0
  assert_eq(trace.samples()[0], [2.0, 0.0])
  let g = n.generator([[2, 0], [1, 1], [0, 2]])
  let counts = g.count_states()
  counts[0][0] = 99
  assert_eq(g.count_states()[0], [2, 0])
  let d = g.evolve([1.0, 0.0, 0.0], 1.0)
  let masses = d.masses()
  masses[0] = 100.0
  assert_true(d.masses()[0] < 1.0)
}
