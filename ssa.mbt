///|
/// Right-continuous count path: states[0] at zero, states[k+1] after event k.
pub struct Path {
  priv times : Array[Double]
  priv channels : Array[Int32]
  priv states : Array[Array[Int32]]
  priv horizon : Double
}

///|
pub fn Path::from_events(
  net : Network,
  initial : Array[Int32],
  times : Array[Double],
  channels : Array[Int32],
  horizon : Double,
) -> Path raise KineticsError {
  net.validate_counts(initial)
  nonnegative(horizon, "path horizon")
  if times.length() != channels.length() || times.length() > 100000 {
    raise Invalid("event arrays mismatch or exceed 100000")
  }
  if times.length() + 1 > 2000000 / net.names.length() {
    raise Budget("path cell budget")
  }
  let states = [initial.copy()]
  let mut state = initial.copy()
  let mut previous = 0.0
  for k in 0..<times.length() {
    finite(times[k], "event time")
    if times[k] <= previous || times[k] > horizon {
      raise Invalid("event times must increase within (0,horizon]")
    }
    let channel = channels[k]
    let _ = net.reaction(channel)
    if net.propensities(state)[channel] <= 0.0 {
      raise Invalid("event has zero propensity")
    }
    state = net.fire(state, channel)
    states.push(state.copy())
    previous = times[k]
  }
  { times: times.copy(), channels: channels.copy(), states, horizon, }
}

///|
pub fn Path::event_times(self : Path) -> Array[Double] {
  self.times.copy()
}

///|
pub fn Path::event_channels(self : Path) -> Array[Int32] {
  self.channels.copy()
}

///|
pub fn Path::count_states(self : Path) -> Array[Array[Int32]] {
  self.states.map(s => s.copy())
}

///|
pub fn Path::end_time(self : Path) -> Double {
  self.horizon
}

///|
pub fn Path::last_counts(self : Path) -> Array[Int32] {
  self.states[self.states.length() - 1].copy()
}

///|
/// Direct Gillespie, including null channels as observed reaction events.
/// Event/clock limits fail explicitly instead of returning a truncated path.
pub fn Network::simulate(
  self : Network,
  initial : Array[Int32],
  horizon : Double,
  seed~ : Int32 = 1,
  max_events~ : Int32 = 100000,
) -> Path raise KineticsError {
  self.validate_counts(initial)
  nonnegative(horizon, "SSA horizon")
  if max_events < 0 || max_events > 100000 {
    raise Invalid("SSA event budget")
  }
  let rng = Random::new(seed)
  let mut state = initial.copy()
  let mut time = 0.0
  let times : Array[Double] = []
  let channels : Array[Int32] = []
  let states = [state.copy()]
  while time < horizon {
    let a = self.propensities(state)
    let mut total = 0.0
    for x in a {
      total += x
      finite(total, "SSA total rate")
    }
    if total == 0.0 {
      break
    }
    let waiting = -@math.ln(rng.uniform()) / total
    if waiting > horizon - time {
      break
    }
    if waiting <= 0.0 || time + waiting <= time {
      raise Numerical("SSA clock cannot advance")
    }
    if times.length() >= max_events ||
      states.length() >= 2000000 / self.names.length() {
      raise Budget("SSA event/state budget exhausted; no partial result")
    }
    let channel = rng.categorical(a)
    state = self.fire(state, channel)
    time += waiting
    times.push(time)
    channels.push(channel)
    states.push(state.copy())
  }
  { times, channels, states, horizon, }
}

///|
/// Step sampling: an event exactly at t is already reflected at t.
pub fn Path::sample(
  self : Path,
  times : Array[Double],
) -> Trace raise KineticsError {
  validate_grid(times)
  if times[times.length() - 1] > self.horizon {
    raise Invalid("sample beyond path horizon")
  }
  if times.length() > 2000000 / self.states[0].length() {
    raise Budget("SSA sample cells")
  }
  let result : Array[Array[Double]] = []
  let mut event = 0
  for t in times {
    while event < self.times.length() && self.times[event] <= t {
      event += 1
    }
    result.push(self.states[event].map(n => n.to_double()))
  }
  Trace::new(times, result)
}
