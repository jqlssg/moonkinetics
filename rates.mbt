///|
fn int_power(x : Double, exponent : Int32) -> Double {
  let mut y = 1.0
  for _ in 0..<exponent {
    y *= x
  }
  y
}

///|
pub fn Network::fluxes(
  self : Network,
  concentrations : Array[Double],
) -> Array[Double] raise KineticsError {
  self.validate_concentrations(concentrations)
  let values : Array[Double] = []
  for r in self.reactions {
    let mut value = r.rate
    for i in 0..<concentrations.length() {
      value *= int_power(concentrations[i], r.input[i])
      finite(value, "mass-action flux")
    }
    values.push(value)
  }
  values
}

///|
pub fn Network::drift(
  self : Network,
  concentrations : Array[Double],
) -> Array[Double] raise KineticsError {
  let fluxes = self.fluxes(concentrations)
  let result = Array::make(self.names.length(), 0.0)
  for j in 0..<fluxes.length() {
    let r = self.reactions[j]
    for i in 0..<result.length() {
      result[i] += (r.output[i] - r.input[i]).to_double() * fluxes[j]
      finite(result[i], "reaction drift")
    }
  }
  result
}

///|
/// Density-dependent convention, including zero-order immigration k*V.
pub fn Network::propensities(
  self : Network,
  state : Array[Int32],
) -> Array[Double] raise KineticsError {
  self.validate_counts(state)
  let values : Array[Double] = []
  for r in self.reactions {
    let mut order = 0
    let mut disabled = false
    for i in 0..<state.length() {
      order += r.input[i]
      if state[i] < r.input[i] {
        disabled = true
      }
    }
    if disabled || r.rate == 0.0 {
      values.push(0.0)
      continue
    }
    let mut value = r.rate * self.volume
    for i in 0..<state.length() {
      for k in 0..<r.input[i] {
        value *= (state[i] - k).to_double() / self.volume
        finite(value, "mass-action propensity")
      }
    }
    // The validated order bound ensures products remain small enough to check.
    if order > 4 {
      raise Invalid("unsupported reaction order")
    }
    finite(value, "propensity")
    values.push(value)
  }
  values
}

///|
pub fn Network::total_propensity(
  self : Network,
  state : Array[Int32],
) -> Double raise KineticsError {
  let mut total = 0.0
  for a in self.propensities(state) {
    total += a
    finite(total, "total propensity")
  }
  total
}

///|
/// Counts drift, not concentration drift; finite-copy corrections are retained.
pub fn Network::count_drift(
  self : Network,
  state : Array[Int32],
) -> Array[Double] raise KineticsError {
  let a = self.propensities(state)
  let result = Array::make(self.names.length(), 0.0)
  for j in 0..<a.length() {
    for i in 0..<result.length() {
      result[i] += (self.reactions[j].output[i] - self.reactions[j].input[i]).to_double() *
        a[j]
      finite(result[i], "count drift")
    }
  }
  result
}
