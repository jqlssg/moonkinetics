///|
/// Exact integral of the right-continuous piecewise-constant count path.
pub fn Path::integrated_counts(
  self : Path,
) -> Array[Double] raise KineticsError {
  let result = Array::make(self.states[0].length(), 0.0)
  let mut previous = 0.0
  for k in 0..<(self.channels.length() + 1) {
    let end = if k == self.channels.length() {
      self.horizon
    } else {
      self.times[k]
    }
    for i in 0..<result.length() {
      result[i] += self.states[k][i].to_double() * (end - previous)
      finite(result[i], "integrated counts")
    }
    previous = end
  }
  result
}

///|
pub fn Network::event_counts(
  self : Network,
  path : Path,
) -> Array[Int32] raise KineticsError {
  verify_path(self, path)
  let counts = Array::make(self.reactions.length(), 0)
  for j in path.channels {
    counts[j] += 1
  }
  counts
}

///|
pub fn Network::integrated_propensities(
  self : Network,
  path : Path,
) -> Array[Double] raise KineticsError {
  verify_path(self, path)
  let result = Array::make(self.reactions.length(), 0.0)
  let mut previous = 0.0
  for k in 0..<(path.channels.length() + 1) {
    let end = if k == path.channels.length() {
      path.horizon
    } else {
      path.times[k]
    }
    let a = self.propensities(path.states[k])
    for j in 0..<a.length() {
      result[j] += a[j] * (end - previous)
      finite(result[j], "integrated propensity")
    }
    previous = end
  }
  result
}

///|
/// Observable expectation of the retained subdistribution, plus absolute
/// Poisson-tail error bound max|f|*tail. Rounding error is not included.
pub fn Generator::expectation(
  self : Generator,
  distribution : Distribution,
  observable : Array[Double],
) -> (Double, Double) raise KineticsError {
  if observable.length() != self.size() ||
    distribution.probabilities.length() != self.size() {
    raise Invalid("observable state shape")
  }
  let mut magnitude = 0.0
  for x in observable {
    finite(x, "observable")
    if x.abs() > magnitude {
      magnitude = x.abs()
    }
  }
  let value = vector_dot(distribution.probabilities, observable)
  let error = magnitude * distribution.tail
  finite(error, "observable tail bound")
  (value, error)
}

///|
pub fn Generator::mean_counts(
  self : Generator,
  distribution : Distribution,
) -> (Array[Double], Array[Double]) raise KineticsError {
  let means : Array[Double] = []
  let bounds : Array[Double] = []
  for i in 0..<self.states[0].length() {
    let observable = self.states.map(s => s[i].to_double())
    let (value, error) = self.expectation(distribution, observable)
    means.push(value)
    bounds.push(error)
  }
  (means, bounds)
}
