///|
/// Deterministic nonlinear seed mapping, not proof of independent/nonoverlap
/// streams. Run index is stable when scheduling is changed.
pub fn ensemble_seed(seed : Int32, run : Int32) -> Int32 raise KineticsError {
  let _ = Random::new(seed)
  if run < 0 || run >= 10000 {
    raise Invalid("ensemble run index 0..9999")
  }
  let x = (seed.to_int64() + run.to_int64() * 104729L) % 2147483646L + 1L
  let y = x * x % 2147483646L + 1L
  (y * y % 2147483646L + 1L).to_int()
}

///|
pub struct Ensemble {
  priv times : Array[Double]
  priv moments : Array[Moments]
}

///|
pub fn Network::ensemble(
  self : Network,
  initial : Array[Int32],
  times : Array[Double],
  runs : Int32,
  seed~ : Int32 = 1,
  max_events~ : Int32 = 10000,
) -> Ensemble raise KineticsError {
  self.validate_counts(initial)
  validate_grid(times)
  if times[0] != 0.0 {
    raise Invalid("ensemble grid starts at zero")
  }
  if runs < 2 || runs > 10000 {
    raise Invalid("ensemble needs 2..10000 runs")
  }
  if max_events < 0 || max_events > 100000 {
    raise Invalid("ensemble per-path budget")
  }
  let n = self.names.length()
  if times.length() > 2000000 / (n * n) {
    raise Budget("ensemble covariance cell budget")
  }
  if max_events > 10000000 / runs {
    raise Budget("ensemble aggregate event budget")
  }
  let moments : Array[Moments] = []
  for _ in times {
    moments.push(Moments::new(n))
  }
  let horizon = times[times.length() - 1]
  for r in 0..<runs {
    let path = self.simulate(
      initial,
      horizon,
      seed=ensemble_seed(seed, r),
      max_events~,
    )
    let samples = path.sample(times)
    for i in 0..<times.length() {
      moments[i].push(samples.values[i])
    }
  }
  { times: times.copy(), moments, }
}

///|
pub fn Ensemble::means(self : Ensemble) -> Trace raise KineticsError {
  let values : Array[Array[Double]] = []
  for m in self.moments {
    values.push(m.average())
  }
  Trace::new(self.times, values)
}

///|
pub fn Ensemble::covariances(
  self : Ensemble,
) -> Array[Array[Array[Double]]] raise KineticsError {
  let result : Array[Array[Array[Double]]] = []
  for m in self.moments {
    result.push(m.covariance())
  }
  result
}

///|
pub fn Ensemble::runs(self : Ensemble) -> Int32 {
  self.moments[0].count
}

///|
pub fn Ensemble::standard_errors(
  self : Ensemble,
) -> Array[Array[Double]] raise KineticsError {
  let result : Array[Array[Double]] = []
  for m in self.moments {
    let covariance = m.covariance()
    let row : Array[Double] = []
    for i in 0..<m.mean.length() {
      // Tiny negative variance would signal a numerical defect, not be hidden.
      nonnegative(covariance[i][i], "sample variance")
      row.push((covariance[i][i] / m.count.to_double()).sqrt())
    }
    result.push(row)
  }
  result
}
