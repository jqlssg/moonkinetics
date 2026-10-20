///|
test "positive underflow is an error, not absorption" {
  let n = Network::new(
    ["A"],
    [Reaction::new([0], [1], 1.0e-320)],
    volume=1.0e-30,
  )
  let bad = try {
    let _ = n.propensities([0])
    false
  } catch {
    Numerical(_) => true
    _ => false
  }
  assert_true(bad)
  let n = Network::new(["A"], [Reaction::new([4], [0], 1.0)])
  let bad = try {
    let _ = n.fluxes([1.0e-100])
    false
  } catch {
    Numerical(_) => true
    _ => false
  }
  assert_true(bad)
}

///|
test "zero reactants short circuit otherwise overflowing products" {
  let n = Network::new(["A", "B"], [Reaction::new([3, 1], [0, 0], 1.0e100)])
  assert_eq(n.fluxes([1.0e100, 0.0]), [0.0])
  assert_eq(n.jacobian([0.0, 0.0]), [[0.0, 0.0], [0.0, 0.0]])
}

///|
test "failed moment update leaves accumulator unchanged" {
  let m = Moments::new(1)
  m.push([1.0])
  m.push([2.0])
  let bad = try {
    m.push([1.0e308])
    false
  } catch {
    Numerical(_) => true
    _ => false
  }
  assert_true(bad)
  assert_eq(m.sample_count(), 2)
  assert_eq(m.average(), [1.5])
  assert_eq(m.covariance(), [[0.5]])
}

///|
test "RNG long independent modular vector" {
  let r = Random::new(1)
  for _ in 0..<10000 {
    let _ = r.next_int()
  }
  assert_eq(r.seed_state(), 1043618065)
}

///|
test "event at horizon is right continuous and has no terminal exposure" {
  let n = Network::new(["A"], [Reaction::new([1], [0], 2.0)])
  let p = Path::from_events(n, [1], [1.0], [0], 1.0)
  assert_eq(p.sample([0.0, 1.0]).samples(), [[1.0], [0.0]])
  assert_true((n.log_likelihood(p) - (@math.ln(2.0) - 2.0)).abs() < 1.0e-14)
  assert_eq(p.integrated_counts(), [1.0])
}

///|
test "null channels remain observed but do not change CTMC" {
  let n = Network::new(["A"], [Reaction::new([0], [0], 2.0)])
  let p = Path::from_events(n, [3], [0.5], [0], 1.0)
  assert_eq(n.event_counts(p), [1])
  assert_eq(n.integrated_propensities(p), [2.0])
  assert_eq(n.generator([[3]]).dense(), [[-0.0]])
  assert_eq(n.generator([[3]]).evolve([1.0], 100.0).masses(), [1.0])
}

///|
test "CTMC semigroup with explicit truncation error and mass policy" {
  let n = Network::new(["A", "B"], [
    Reaction::new([1, 0], [0, 1], 2.0),
    Reaction::new([0, 1], [1, 0], 1.0),
  ])
  let g = n.generator([[1, 0], [0, 1]])
  let first = g.evolve([1.0, 0.0], 1.0)
  let twice = g.evolve(first.masses(), 2.0)
  let once = g.evolve([1.0, 0.0], 3.0)
  for i in 0..<2 {
    assert_true((twice.masses()[i] - once.masses()[i]).abs() < 3.0e-12)
  }
  let extra_mass = g.evolve([1.0 + 5.0e-11, 0.0], 3.0)
  let m = extra_mass.masses()
  assert_true(
    (1.0 + 5.0e-11 - m[0] - m[1]).abs() < extra_mass.tail_bound() + 1.0e-14,
  )
}

///|
test "reject negative count" {
  let bad = try {
    let n = Network::new(["A"], [])
    n.validate_counts([-1])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject count shape" {
  let bad = try {
    let n = Network::new(["A"], [])
    n.validate_counts([])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject nonfinite concentration" {
  let bad = try {
    let n = Network::new(["A"], [])
    n.validate_concentrations([0.0 / 0.0])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject infinite rate" {
  let bad = try {
    let _ = Reaction::new([1], [0], 1.0 / 0.0)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject negative volume" {
  let bad = try {
    let _ = Network::new(["A"], [], volume=-1.0)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject empty species" {
  let bad = try {
    let _ = Network::new([], [])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject state cap fire" {
  let bad = try {
    let n = Network::new(["A"], [Reaction::new([0], [1], 1.0)])
    let _ = n.fire([1000000000], 0)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject bad categorical mass" {
  let bad = try {
    let r = Random::new(1)
    let _ = r.categorical([0.0, 0.0])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject invalid random seed" {
  let bad = try {
    let _ = Random::new(0)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject ragged trace" {
  let bad = try {
    let _ = Trace::new([0.0, 1.0], [[1.0], []])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject duplicate grid" {
  let bad = try {
    let _ = Trace::new([0.0, 0.0], [[1.0], [1.0]])
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject impossible zero rate event" {
  let bad = try {
    let n = Network::new(["A"], [Reaction::new([1], [0], 0.0)])
    let _ = Path::from_events(n, [1], [0.5], [0], 1.0)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject invalid probability" {
  let bad = try {
    let n = Network::new(["A"], [])
    let _ = n.generator([[0]]).evolve([0.5], 1.0)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject uniformization budget" {
  let bad = try {
    let n = Network::new(["A"], [Reaction::new([1], [0], 1.0)])
    let _ = n.generator([[0], [1]]).evolve([0.0, 1.0], 100.0, max_terms=1)
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject nested duplicate JSON" {
  let bad = try {
    let _ = guarded_json("{\"a\":{\"x\":1,\"x\":2}}")
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject excessive nesting" {
  let bad = try {
    let _ = guarded_json("[[[[[[[]]]]]]]")
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "reject unknown model key" {
  let bad = try {
    let _ = network_from_json(
      "{\"schema\":1,\"species\":[\"A\"],\"volume\":1,\"reactions\":[],\"extra\":1}",
    )
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
