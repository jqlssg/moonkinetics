///|
pub struct Generator {
  priv states : Array[Array[Int32]]
  priv edges : Array[Array[(Int32, Double)]]
  priv exits : Array[Double]
}

///|
fn count_key(state : Array[Int32]) -> String {
  state.map(n => n.to_string()).join(",")
}

///|
/// Constructs a closed CTMC from an explicit user-supplied state set. A missing
/// positive-rate successor is an error, never an undocumented reflecting wall.
pub fn Network::generator(
  self : Network,
  states : Array[Array[Int32]],
) -> Generator raise KineticsError {
  if states.is_empty() || states.length() > 4096 {
    raise Budget("CTMC needs 1..4096 explicit states")
  }
  if states.length() > 2000000 / self.names.length() {
    raise Budget("CTMC state cells")
  }
  let indices : Map[String, Int32] = Map([])
  for i in 0..<states.length() {
    self.validate_counts(states[i])
    let key = count_key(states[i])
    if indices.contains(key) {
      raise Invalid("duplicate CTMC state")
    }
    indices[key] = i
  }
  let rows : Array[Array[(Int32, Double)]] = []
  let exits : Array[Double] = []
  let mut edge_count = 0
  for i in 0..<states.length() {
    let a = self.propensities(states[i])
    let edges : Array[(Int32, Double)] = []
    let mut total = 0.0
    for j in 0..<a.length() {
      if a[j] == 0.0 {
        continue
      }
      let next = self.fire(states[i], j)
      let target = match indices.get(count_key(next)) {
        Some(index) => index
        None =>
          raise Invalid(
            "CTMC state set is not closed under positive-rate reactions",
          )
      }
      if target == i {
        continue
      }
      edges.push((target, a[j]))
      total += a[j]
      finite(total, "CTMC exit rate")
      edge_count += 1
      if edge_count > 1000000 {
        raise Budget("CTMC edge budget")
      }
    }
    rows.push(edges)
    exits.push(total)
  }
  { states: states.map(s => s.copy()), edges: rows, exits, }
}

///|
pub fn Generator::size(self : Generator) -> Int32 {
  self.states.length()
}

///|
pub fn Generator::count_states(self : Generator) -> Array[Array[Int32]] {
  self.states.map(s => s.copy())
}

///|
pub fn Generator::exit_rates(self : Generator) -> Array[Double] {
  self.exits.copy()
}

///|
/// Parallel channels are retained as separate edges; dense() sums them.
pub fn Generator::row_edges(
  self : Generator,
  row : Int32,
) -> Array[(Int32, Double)] raise KineticsError {
  if row < 0 || row >= self.size() {
    raise Invalid("generator row index")
  }
  self.edges[row].copy()
}

///|
pub fn Generator::dense(
  self : Generator,
) -> Array[Array[Double]] raise KineticsError {
  if self.size() > 512 {
    raise Budget("dense generator limited to 512 states")
  }
  let matrix = Array::makei(self.size(), _ => Array::make(self.size(), 0.0))
  for i in 0..<self.size() {
    matrix[i][i] = -self.exits[i]
    for edge in self.edges[i] {
      matrix[i][edge.0] += edge.1
    }
  }
  matrix
}

///|
fn validate_probability(
  p : Array[Double],
  size : Int32,
) -> Unit raise KineticsError {
  if p.length() != size {
    raise Invalid("probability shape")
  }
  let mut sum = 0.0
  for x in p {
    nonnegative(x, "probability")
    sum += x
  }
  if (sum - 1.0).abs() > 1.0e-10 {
    raise Invalid("probability mass must sum to one")
  }
}

///|
/// Forward master-equation derivative for a row probability vector p.
pub fn Generator::derivative(
  self : Generator,
  probabilities : Array[Double],
) -> Array[Double] raise KineticsError {
  validate_probability(probabilities, self.size())
  let result = Array::make(self.size(), 0.0)
  for i in 0..<self.size() {
    result[i] -= probabilities[i] * self.exits[i]
    for edge in self.edges[i] {
      result[edge.0] += probabilities[i] * edge.1
    }
  }
  for x in result {
    finite(x, "probability derivative")
  }
  result
}
