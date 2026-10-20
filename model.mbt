///|
pub suberror KineticsError {
  Invalid(String)
  Numerical(String)
  Budget(String)
} derive(@debug.Debug)

///|
fn finite(x : Double, name : String) -> Unit raise KineticsError {
  if x.is_nan() || x.is_inf() {
    raise Numerical(name + " must be finite")
  }
}

///|
fn nonnegative(x : Double, name : String) -> Unit raise KineticsError {
  finite(x, name)
  if x < 0.0 {
    raise Invalid(name + " must be nonnegative")
  }
}

///|
pub struct Reaction {
  input : Array[Int32]
  output : Array[Int32]
  rate : Double
} derive(@debug.Debug)

///|
/// Rate convention: deterministic k*product(c^r), stochastic
/// k*V^(1-order)*product(falling(n,r)); no factorial divisor.
pub fn Reaction::new(
  input : Array[Int32],
  output : Array[Int32],
  rate : Double,
) -> Reaction raise KineticsError {
  if input.is_empty() ||
    input.length() != output.length() ||
    input.length() > 64 {
    raise Invalid("reaction shape: 1..64 matching species")
  }
  nonnegative(rate, "rate")
  if rate > 1.0e100 {
    raise Invalid("rate too large")
  }
  let mut order = 0
  for i in 0..<input.length() {
    if input[i] < 0 || input[i] > 4 || output[i] < 0 || output[i] > 1000000 {
      raise Invalid("stoichiometric coefficient outside policy")
    }
    order += input[i]
  }
  if order > 4 {
    raise Invalid("only total orders zero through four supported")
  }
  { input: input.copy(), output: output.copy(), rate, }
}

///|
pub fn Reaction::reactants(self : Reaction) -> Array[Int32] {
  self.input.copy()
}

///|
pub fn Reaction::products(self : Reaction) -> Array[Int32] {
  self.output.copy()
}

///|
pub fn Reaction::constant(self : Reaction) -> Double {
  self.rate
}

///|
pub struct Network {
  names : Array[String]
  reactions : Array[Reaction]
  volume : Double
} derive(@debug.Debug)

///|
pub fn Network::new(
  names : Array[String],
  reactions : Array[Reaction],
  volume~ : Double = 1.0,
) -> Network raise KineticsError {
  if names.is_empty() || names.length() > 64 || reactions.length() > 1024 {
    raise Budget("network policy: 1..64 species, at most 1024 reactions")
  }
  finite(volume, "volume")
  if volume < 1.0e-30 || volume > 1.0e30 {
    raise Invalid("volume outside 1e-30..1e30")
  }
  for i in 0..<names.length() {
    if names[i].is_empty() || names[i].length() > 128 {
      raise Invalid("invalid species name")
    }
    for j in 0..<i {
      if names[i] == names[j] {
        raise Invalid("duplicate species name")
      }
    }
  }
  for r in reactions {
    if r.input.length() != names.length() {
      raise Invalid("reaction/network dimension mismatch")
    }
  }
  { names: names.copy(), reactions: reactions.copy(), volume, }
}

///|
pub fn Network::species(self : Network) -> Array[String] {
  self.names.copy()
}

///|
pub fn Network::species_count(self : Network) -> Int32 {
  self.names.length()
}

///|
pub fn Network::reaction_count(self : Network) -> Int32 {
  self.reactions.length()
}

///|
pub fn Network::system_volume(self : Network) -> Double {
  self.volume
}

///|
pub fn Network::reaction(
  self : Network,
  index : Int32,
) -> Reaction raise KineticsError {
  if index < 0 || index >= self.reactions.length() {
    raise Invalid("reaction index")
  }
  self.reactions[index]
}

///|
pub fn Network::validate_counts(
  self : Network,
  state : Array[Int32],
) -> Unit raise KineticsError {
  if state.length() != self.names.length() {
    raise Invalid("count state shape")
  }
  for n in state {
    if n < 0 || n > 1000000000 {
      raise Invalid("count outside 0..1000000000")
    }
  }
}

///|
pub fn Network::validate_concentrations(
  self : Network,
  state : Array[Double],
) -> Unit raise KineticsError {
  if state.length() != self.names.length() {
    raise Invalid("concentration state shape")
  }
  for x in state {
    nonnegative(x, "concentration")
    if x > 1.0e100 {
      raise Invalid("concentration too large")
    }
  }
}

///|
pub extend KineticsError with @debug.Debug::{to_repr}

///|
pub extend Reaction with @debug.Debug::{to_repr}

///|
pub extend Network with @debug.Debug::{to_repr}
