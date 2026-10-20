///|
pub struct Distribution {
  priv probabilities : Array[Double]
  priv tail : Double
  priv terms : Int32
  priv segments : Int32
}

///|
pub fn Distribution::masses(self : Distribution) -> Array[Double] {
  self.probabilities.copy()
}

///|
/// Conservative Poisson truncation bound, not a floating-point rounding bound.
pub fn Distribution::tail_bound(self : Distribution) -> Double {
  self.tail
}

///|
pub fn Distribution::work_terms(self : Distribution) -> Int32 {
  self.terms
}

///|
pub fn Distribution::segment_count(self : Distribution) -> Int32 {
  self.segments
}

///|
fn uniform_step(
  g : Generator,
  p : Array[Double],
  lambda : Double,
) -> Array[Double] raise KineticsError {
  let result = Array::make(g.size(), 0.0)
  for i in 0..<g.size() {
    result[i] += p[i] * (1.0 - g.exits[i] / lambda)
    for edge in g.edges[i] {
      result[edge.0] += p[i] * (edge.1 / lambda)
    }
  }
  for x in result {
    nonnegative(x, "uniformized mass")
  }
  result
}

///|
/// exp(tQ) by Poisson uniformization. Time splitting keeps each mean <=16.
/// No renormalization hides the explicitly reported missing tail mass.
pub fn Generator::evolve(
  self : Generator,
  initial : Array[Double],
  time : Double,
  tolerance~ : Double = 1.0e-12,
  max_terms~ : Int32 = 100000,
  max_work~ : Int32 = 100000000,
) -> Distribution raise KineticsError {
  validate_probability(initial, self.size())
  nonnegative(time, "CTMC time")
  finite(tolerance, "Poisson tolerance")
  if tolerance < 1.0e-15 || tolerance > 1.0e-3 {
    raise Invalid("tolerance outside 1e-15..1e-3")
  }
  if max_terms < 1 ||
    max_terms > 1000000 ||
    max_work < 1 ||
    max_work > 100000000 {
    raise Invalid("uniformization work budget")
  }
  let mut lambda = 0.0
  let mut per_term = self.size()
  for i in 0..<self.size() {
    if self.exits[i] > lambda {
      lambda = self.exits[i]
    }
    per_term += self.edges[i].length()
  }
  if lambda == 0.0 || time == 0.0 {
    return { probabilities: initial.copy(), tail: 0.0, terms: 0, segments: 0, }
  }
  let mean = lambda * time
  finite(mean, "CTMC rate-time product")
  if mean > 160000.0 {
    raise Budget("uniformization at most 10000 time segments")
  }
  let mut segments = (mean / 16.0).to_int()
  if segments.to_double() * 16.0 < mean {
    segments += 1
  }
  if segments == 0 {
    segments = 1
  }
  let mu = mean / segments.to_double()
  let epsilon = tolerance / segments.to_double()
  let mut p = initial.copy()
  let mut tail = 0.0
  let mut terms = 0
  let mut work = 0
  for _ in 0..<segments {
    let mut term = p.copy()
    let mut weight = @math.exp(-mu)
    let out = term.map(x => x * weight)
    let mut k = 0
    while true {
      if terms >= max_terms || per_term > max_work - work {
        raise Budget("uniformization term/work budget exhausted")
      }
      terms += 1
      work += per_term
      let next_weight = weight * mu / (k + 1).to_double()
      let ratio = mu / (k + 2).to_double()
      if ratio < 1.0 {
        let bound = next_weight / (1.0 - ratio)
        if bound <= epsilon {
          tail += bound
          break
        }
      }
      term = uniform_step(self, term, lambda)
      weight = next_weight
      k += 1
      for i in 0..<out.length() {
        out[i] += weight * term[i]
      }
    }
    p = out
  }
  for x in p {
    nonnegative(x, "evolved probability")
  }
  { probabilities: p, tail, terms, segments, }
}
