///|
fn vector_dot(
  a : Array[Double],
  b : Array[Double],
) -> Double raise KineticsError {
  if a.length() != b.length() {
    raise Invalid("dot product shape")
  }
  let mut sum = 0.0
  for i in 0..<a.length() {
    finite(a[i], "dot weight")
    finite(b[i], "dot value")
    sum += a[i] * b[i]
    finite(sum, "dot product")
  }
  sum
}

///|
/// Residual for a proposed linear conserved quantity w^T S.
pub fn Network::conservation_residual(
  self : Network,
  weights : Array[Double],
) -> Array[Double] raise KineticsError {
  if weights.length() != self.names.length() {
    raise Invalid("conservation shape")
  }
  let s = self.stoichiometry()
  let result = Array::make(self.reactions.length(), 0.0)
  for j in 0..<result.length() {
    for i in 0..<weights.length() {
      finite(weights[i], "conservation weight")
      result[j] += weights[i] * s[i][j].to_double()
      finite(result[j], "conservation residual")
    }
  }
  result
}

///|
/// Floating-point basis of ker(S^T), not an exact integer certificate.
pub fn Network::conservation_basis(
  self : Network,
  tolerance~ : Double = 1.0e-10,
) -> Array[Array[Double]] raise KineticsError {
  finite(tolerance, "tolerance")
  if tolerance <= 0.0 || tolerance > 1.0e-3 {
    raise Invalid("conservation tolerance")
  }
  let n = self.names.length()
  let a = self.reactions.map(r => {
    Array::makei(n, i => (r.output[i] - r.input[i]).to_double())
  })
  let pivots : Array[Int32] = []
  let mut row = 0
  for col in 0..<n {
    let mut best = row
    let mut magnitude = tolerance
    for k in row..<a.length() {
      if a[k][col].abs() > magnitude {
        best = k
        magnitude = a[k][col].abs()
      }
    }
    if magnitude == tolerance {
      continue
    }
    let temp = a[row]
    a[row] = a[best]
    a[best] = temp
    let pivot = a[row][col]
    for j in col..<n {
      a[row][j] /= pivot
    }
    for k in 0..<a.length() {
      if k == row {
        continue
      }
      let factor = a[k][col]
      for j in col..<n {
        a[k][j] -= factor * a[row][j]
      }
    }
    pivots.push(col)
    row += 1
    if row == a.length() {
      break
    }
  }
  let basis : Array[Array[Double]] = []
  for col in 0..<n {
    if pivots.contains(col) {
      continue
    }
    let vector = Array::make(n, 0.0)
    vector[col] = 1.0
    for k in 0..<pivots.length() {
      vector[pivots[k]] = -a[k][col]
    }
    for x in self.conservation_residual(vector) {
      if x.abs() > tolerance * 100.0 {
        raise Numerical("unstable conservation basis")
      }
    }
    basis.push(vector)
  }
  basis
}

///|
pub fn conserved_value(
  weights : Array[Double],
  values : Array[Double],
) -> Double raise KineticsError {
  vector_dot(weights, values)
}
