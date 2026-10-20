///|
/// Checks compatibility by replay, not merely matching vector dimensions.
fn verify_path(net : Network, path : Path) -> Unit raise KineticsError {
  net.validate_counts(path.states[0])
  let mut state = path.states[0].copy()
  for k in 0..<path.channels.length() {
    let channel = path.channels[k]
    let _ = net.reaction(channel)
    if net.propensities(state)[channel] <= 0.0 {
      raise Invalid("path event impossible under this network")
    }
    state = net.fire(state, channel)
    if state != path.states[k + 1] {
      raise Invalid("path stoichiometry differs from network")
    }
  }
}

///|
/// Density of a fully observed channel-labelled path. Includes the final
/// no-event survival term; not the likelihood of sparsely sampled counts.
pub fn Network::log_likelihood(
  self : Network,
  path : Path,
) -> Double raise KineticsError {
  verify_path(self, path)
  let mut result = 0.0
  let mut previous = 0.0
  for k in 0..<path.channels.length() {
    let a = self.propensities(path.states[k])
    let mut total = 0.0
    for x in a {
      total += x
      finite(total, "likelihood total")
    }
    result += @math.ln(a[path.channels[k]]) - total * (path.times[k] - previous)
    finite(result, "path log likelihood")
    previous = path.times[k]
  }
  result -= self.total_propensity(path.last_counts()) *
    (path.horizon - previous)
  finite(result, "terminal survival likelihood")
  result
}

///|
/// Score d(log L)/dk_j for positive constants of linear rate channels.
pub fn Network::path_score(
  self : Network,
  path : Path,
) -> Array[Double] raise KineticsError {
  verify_path(self, path)
  for r in self.reactions {
    if r.rate <= 0.0 {
      raise Invalid("path score requires positive rate constants")
    }
  }
  let result = Array::make(self.reactions.length(), 0.0)
  let mut previous = 0.0
  for k in 0..<(path.channels.length() + 1) {
    let next = if k == path.channels.length() {
      path.horizon
    } else {
      path.times[k]
    }
    let a = self.propensities(path.states[k])
    for j in 0..<a.length() {
      result[j] -= a[j] / self.reactions[j].rate * (next - previous)
      finite(result[j], "path score")
    }
    if k < path.channels.length() {
      let j = path.channels[k]
      result[j] += 1.0 / self.reactions[j].rate
      finite(result[j], "event score")
    }
    previous = next
  }
  result
}
