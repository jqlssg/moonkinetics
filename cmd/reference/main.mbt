///|
fn main raise {
  for v in 0..<80 {
    let order = v % 5
    let input = [v % (order + 1), order - v % (order + 1)]
    let output = [v * 3 % 4, v * 7 % 5]
    let counts = [v % 9, v * 5 % 11]
    let rate = (v % 7 + 1).to_double() / 3.0
    let volume = (v % 4 + 1).to_double() / 2.0
    let net = @k.Network::new(
      ["A", "B"],
      [@k.Reaction::new(input, output, rate)],
      volume~,
    )
    let concentrations = net.concentrations(counts)
    let vector : Json = {
      "kind": "rate",
      "input": input,
      "output": output,
      "counts": counts,
      "rate": rate,
      "volume": volume,
      "concentrations": concentrations,
      "propensity": net.propensities(counts)[0],
      "flux": net.fluxes(concentrations)[0],
      "drift": net.drift(concentrations),
      "jacobian": net.jacobian(concentrations),
      "diffusion": net.diffusion(counts),
    }
    println(vector.stringify())
  }
  for n in [1, 2, 5, 10, 20] {
    let net = @k.Network::new(["A", "B"], [
      @k.Reaction::new([1, 0], [0, 1], 2.0),
      @k.Reaction::new([0, 1], [1, 0], 1.0),
    ])
    let g = net.generator(Array::makei(n + 1, a => [a, n - a]))
    let p = Array::make(n + 1, 0.0)
    p[n] = 1.0
    for t in [0.0, 0.1, 1.0, 10.0, 100.0] {
      let d = g.evolve(p, t)
      let vector : Json = {
        "kind": "binomial",
        "n": n,
        "time": t,
        "masses": d.masses(),
        "tail": d.tail_bound(),
      }
      println(vector.stringify())
    }
  }
  let net = @k.Network::new(
    ["A", "B"],
    [
      @k.Reaction::new([2, 0], [0, 1], 0.2),
      @k.Reaction::new([0, 1], [2, 0], 0.7),
    ],
    volume=3.0,
  )
  for seed in [1, 7, 42, 99991] {
    let path = net.simulate([12, 0], 4.0, seed~)
    let vector : Json = {
      "kind": "ssa",
      "seed": seed,
      "path": path.to_json(),
      "states": path.count_states(),
      "log_likelihood": net.log_likelihood(path),
      "score": net.path_score(path),
    }
    println(vector.stringify())
  }
}
