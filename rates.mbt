///|
// All factors are nonnegative; pre-scan zeros to avoid spurious overflow.
// Reject unrepresentable positive intermediates rather than invent absorption.
fn monomial(
  rate : Double,
  state : Array[Double],
  powers : Array[Int32],
) -> Double raise KineticsError {
  if rate == 0.0 {
    return 0.0
  }
  for i in 0..<state.length() {
    if powers[i] > 0 && state[i] == 0.0 {
      return 0.0
    }
  }
  let mut result = rate
  for i in 0..<state.length() {
    for _ in 0..<powers[i] {
      result *= state[i]
      finite(result, "mass-action monomial")
      if result == 0.0 {
        raise Numerical("mass-action monomial underflow")
      }
    }
  }
  result
}

///|
pub fn Network::fluxes(
  self : Network,
  concentrations : Array[Double],
) -> Array[Double] raise KineticsError {
  self.validate_concentrations(concentrations)
  let values : Array[Double] = []
  for r in self.reactions {
    let value = monomial(r.rate, concentrations, r.input)
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
    finite(value, "propensity scale")
    if value == 0.0 {
      raise Numerical("propensity scale underflow")
    }
    for i in 0..<state.length() {
      for k in 0..<r.input[i] {
        value *= (state[i] - k).to_double() / self.volume
        finite(value, "mass-action propensity")
        if value == 0.0 {
          raise Numerical("mass-action propensity underflow")
        }
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
