///|
/// Rows are species, columns are reaction channels.
pub fn Network::stoichiometry(self : Network) -> Array[Array[Int32]] {
  Array::makei(self.names.length(), i => {
    Array::makei(self.reactions.length(), j => {
      self.reactions[j].output[i] - self.reactions[j].input[i]
    })
  })
}

///|
pub fn Network::enabled(
  self : Network,
  state : Array[Int32],
  channel : Int32,
) -> Bool raise KineticsError {
  self.validate_counts(state)
  let r = self.reaction(channel)
  for i in 0..<state.length() {
    if state[i] < r.input[i] {
      return false
    }
  }
  true
}

///|
/// Atomic count update: a disabled reaction never returns a partial state.
pub fn Network::fire(
  self : Network,
  state : Array[Int32],
  channel : Int32,
) -> Array[Int32] raise KineticsError {
  if !self.enabled(state, channel) {
    raise Invalid("reaction is disabled")
  }
  let r = self.reactions[channel]
  let next = state.copy()
  for i in 0..<state.length() {
    let n = state[i] - r.input[i]
    if r.output[i] > 1000000000 - n {
      raise Budget("count overflow policy")
    }
    next[i] = n + r.output[i]
  }
  next
}

///|
pub fn Network::concentrations(
  self : Network,
  counts : Array[Int32],
) -> Array[Double] raise KineticsError {
  self.validate_counts(counts)
  counts.map(n => n.to_double() / self.volume)
}

///|
/// Nearest integer conversion is explicit, never an implicit simulation step.
pub fn Network::rounded_counts(
  self : Network,
  concentrations : Array[Double],
) -> Array[Int32] raise KineticsError {
  self.validate_concentrations(concentrations)
  let result : Array[Int32] = []
  for c in concentrations {
    let x = c * self.volume
    finite(x, "scaled count")
    if x > 1000000000.0 {
      raise Budget("scaled count exceeds count policy")
    }
    result.push((x + 0.5).to_int())
  }
  result
}

///|
/// Canonical channel order. Species with identical net changes still remain
/// distinct channels because their event likelihoods can be different.
pub fn Network::orders(self : Network) -> Array[Int32] {
  self.reactions.map(r => {
    let mut n = 0
    for x in r.input {
      n += x
    }
    n
  })
}

///|
pub fn Network::with_rate(
  self : Network,
  channel : Int32,
  rate : Double,
) -> Network raise KineticsError {
  let r = self.reaction(channel)
  let reactions = self.reactions.copy()
  reactions[channel] = Reaction::new(r.input, r.output, rate)
  Network::new(self.names, reactions, volume=self.volume)
}
