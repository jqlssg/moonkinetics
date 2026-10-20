///|
fn shifted(
  state : Array[Double],
  slope : Array[Double],
  h : Double,
) -> Array[Double] raise KineticsError {
  let result = Array::makei(state.length(), i => state[i] + h * slope[i])
  for x in result {
    nonnegative(x, "RK stage: reduce max_step for negative concentrations")
  }
  result
}

///|
fn rk4_step(
  net : Network,
  state : Array[Double],
  h : Double,
) -> Array[Double] raise KineticsError {
  let a = net.drift(state)
  let b = net.drift(shifted(state, a, h * 0.5))
  let c = net.drift(shifted(state, b, h * 0.5))
  let d = net.drift(shifted(state, c, h))
  let next = Array::makei(state.length(), i => {
    state[i] + h * (a[i] / 6.0 + b[i] / 3.0 + c[i] / 3.0 + d[i] / 6.0)
  })
  net.validate_concentrations(next)
  next
}

///|
/// Fixed maximum step RK4 baseline, exactly landing on requested output times.
/// No stiffness detection, clipping, error estimator or adaptive accuracy claim.
pub fn Network::integrate(
  self : Network,
  initial : Array[Double],
  times : Array[Double],
  max_step~ : Double = 0.01,
  max_steps~ : Int32 = 1000000,
) -> Trace raise KineticsError {
  self.validate_concentrations(initial)
  validate_grid(times)
  if times[0] != 0.0 {
    raise Invalid("ODE grid must start at zero")
  }
  finite(max_step, "max_step")
  if max_step <= 0.0 || max_steps < 1 || max_steps > 1000000 {
    raise Invalid("ODE step budget")
  }
  if times.length() > 2000000 / self.names.length() {
    raise Budget("ODE output budget")
  }
  let mut state = initial.copy()
  let values = [state.copy()]
  let mut time = 0.0
  let mut steps = 0
  for k in 1..<times.length() {
    let target = times[k]
    while time < target {
      if steps >= max_steps {
        raise Budget("ODE steps exhausted; no partial result")
      }
      let remaining = target - time
      let h = if remaining < max_step { remaining } else { max_step }
      let next_time = if h == remaining { target } else { time + h }
      if next_time <= time {
        raise Numerical("ODE clock cannot advance")
      }
      state = rk4_step(self, state, h)
      time = next_time
      steps += 1
    }
    values.push(state.copy())
  }
  Trace::new(times, values)
}
