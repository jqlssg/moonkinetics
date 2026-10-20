///|
/// Online vector moments; pushes and merges commit only after all checks.
pub struct Moments {
  priv mut count : Int32
  priv mut mean : Array[Double]
  priv mut m2 : Array[Array[Double]]
}

///|
pub fn Moments::new(dimension : Int32) -> Moments raise KineticsError {
  if dimension < 1 || dimension > 64 {
    raise Invalid("moment dimension 1..64")
  }
  {
    count: 0,
    mean: Array::make(dimension, 0.0),
    m2: Array::makei(dimension, _ => Array::make(dimension, 0.0)),
  }
}

///|
pub fn Moments::push(
  self : Moments,
  sample : Array[Double],
) -> Unit raise KineticsError {
  if sample.length() != self.mean.length() {
    raise Invalid("moment sample shape")
  }
  if self.count >= 100000000 {
    raise Budget("moment count budget")
  }
  for x in sample {
    finite(x, "moment sample")
  }
  let count = self.count + 1
  let delta = Array::makei(sample.length(), i => sample[i] - self.mean[i])
  let mean = Array::makei(sample.length(), i => {
    self.mean[i] + delta[i] / count.to_double()
  })
  let m2 = self.m2.map(row => row.copy())
  for i in 0..<sample.length() {
    finite(mean[i], "updated mean")
    for j in 0..<sample.length() {
      m2[i][j] += delta[i] * (sample[j] - mean[j])
      finite(m2[i][j], "updated moment")
    }
  }
  self.count = count
  self.mean = mean
  self.m2 = m2
}

///|
pub fn Moments::sample_count(self : Moments) -> Int32 {
  self.count
}

///|
pub fn Moments::average(self : Moments) -> Array[Double] raise KineticsError {
  if self.count == 0 {
    raise Invalid("mean requires samples")
  }
  self.mean.copy()
}

///|
pub fn Moments::covariance(
  self : Moments,
) -> Array[Array[Double]] raise KineticsError {
  if self.count < 2 {
    raise Invalid("sample covariance requires two samples")
  }
  self.m2.map(row => row.map(x => x / (self.count - 1).to_double()))
}

///|
/// Parallel Welford merge, including merging an accumulator with itself.
pub fn Moments::merge(
  self : Moments,
  other : Moments,
) -> Unit raise KineticsError {
  if self.mean.length() != other.mean.length() {
    raise Invalid("moment merge shape")
  }
  if other.count == 0 {
    return
  }
  if other.count > 100000000 - self.count {
    raise Budget("merged moment count budget")
  }
  if self.count == 0 {
    self.count = other.count
    self.mean = other.mean.copy()
    self.m2 = other.m2.map(row => row.copy())
    return
  }
  let count = self.count + other.count
  let n1 = self.count.to_double()
  let n2 = other.count.to_double()
  let delta = Array::makei(self.mean.length(), i => other.mean[i] - self.mean[i])
  let mean = Array::makei(self.mean.length(), i => {
    self.mean[i] + delta[i] * (n2 / count.to_double())
  })
  let m2 = self.m2.map(row => row.copy())
  for i in 0..<mean.length() {
    finite(mean[i], "merged mean")
    for j in 0..<mean.length() {
      m2[i][j] += other.m2[i][j] +
        delta[i] * delta[j] * (n1 * (n2 / count.to_double()))
      finite(m2[i][j], "merged covariance")
    }
  }
  self.count = count
  self.mean = mean
  self.m2 = m2
}
