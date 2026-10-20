///|
/// Park-Miller minimal standard using Schrage arithmetic, not cryptographic.
pub struct Random {
  priv mut state : Int32
}

///|
pub fn Random::new(seed : Int32) -> Random raise KineticsError {
  if seed <= 0 || seed >= 2147483647 {
    raise Invalid("seed must be 1..2147483646")
  }
  { state: seed, }
}

///|
pub fn Random::next_int(self : Random) -> Int32 {
  let hi = self.state / 127773
  let lo = self.state % 127773
  let next = 16807 * lo - 2836 * hi
  self.state = if next <= 0 { next + 2147483647 } else { next }
  self.state
}

///|
pub fn Random::uniform(self : Random) -> Double {
  self.next_int().to_double() / 2147483647.0
}

///|
pub fn Random::seed_state(self : Random) -> Int32 {
  self.state
}

///|
pub fn Random::exponential(
  self : Random,
  rate : Double,
) -> Double raise KineticsError {
  finite(rate, "exponential rate")
  if rate <= 0.0 {
    raise Invalid("exponential rate must be positive")
  }
  let value = -@math.ln(self.uniform()) / rate
  finite(value, "exponential waiting time")
  if value <= 0.0 {
    raise Numerical("waiting time underflow")
  }
  value
}

///|
pub fn Random::categorical(
  self : Random,
  weights : Array[Double],
) -> Int32 raise KineticsError {
  if weights.is_empty() || weights.length() > 100000 {
    raise Budget("categorical size")
  }
  let mut total = 0.0
  let mut last = -1
  for i in 0..<weights.length() {
    nonnegative(weights[i], "categorical weight")
    if weights[i] > 0.0 {
      last = i
    }
    total += weights[i]
    finite(total, "categorical sum")
  }
  if total <= 0.0 {
    raise Invalid("categorical needs positive mass")
  }
  let target = self.uniform() * total
  let mut cumulative = 0.0
  for i in 0..<weights.length() {
    cumulative += weights[i]
    if target < cumulative {
      return i
    }
  }
  // Floating sum rounding can only select the last positive weight.
  last
}
