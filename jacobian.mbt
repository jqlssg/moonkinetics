///|
/// Analytic derivative without dividing by concentration: valid at zero.
pub fn Network::jacobian(
  self : Network,
  state : Array[Double],
) -> Array[Array[Double]] raise KineticsError {
  self.validate_concentrations(state)
  let n = self.names.length()
  let matrix = Array::makei(n, _ => Array::make(n, 0.0))
  for r in self.reactions {
    for col in 0..<n {
      let power = r.input[col]
      if power == 0 {
        continue
      }
      let mut d = r.rate * power.to_double()
      for i in 0..<n {
        d *= int_power(state[i], r.input[i] - (if i == col { 1 } else { 0 }))
        finite(d, "flux derivative")
      }
      for row in 0..<n {
        matrix[row][col] += (r.output[row] - r.input[row]).to_double() * d
        finite(matrix[row][col], "Jacobian")
      }
    }
  }
  matrix
}

///|
/// Noise covariance S*diag(a)*S^T in count units per time.
pub fn Network::diffusion(
  self : Network,
  counts : Array[Int32],
) -> Array[Array[Double]] raise KineticsError {
  let a = self.propensities(counts)
  let s = self.stoichiometry()
  let n = self.names.length()
  let result = Array::makei(n, _ => Array::make(n, 0.0))
  for i in 0..<n {
    for k in 0..<n {
      let mut sum = 0.0
      for j in 0..<a.length() {
        sum += s[i][j].to_double() * s[k][j].to_double() * a[j]
        finite(sum, "diffusion")
      }
      result[i][k] = sum
    }
  }
  result
}
