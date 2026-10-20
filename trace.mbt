///|
fn validate_grid(times : Array[Double]) -> Unit raise KineticsError {
  if times.is_empty() || times.length() > 100001 {
    raise Budget("sample grid needs 1..100001 points")
  }
  for i in 0..<times.length() {
    nonnegative(times[i], "sample time")
    if i > 0 && times[i] <= times[i - 1] {
      raise Invalid("sample times must strictly increase")
    }
  }
}

///|
/// Continuous sampled trajectory; linear interpolation is not SSA interpolation.
pub struct Trace {
  priv times : Array[Double]
  priv values : Array[Array[Double]]
}

///|
pub fn Trace::new(
  times : Array[Double],
  values : Array[Array[Double]],
) -> Trace raise KineticsError {
  validate_grid(times)
  if values.length() != times.length() ||
    values[0].is_empty() ||
    values[0].length() > 64 {
    raise Invalid("trace shape")
  }
  if values.length() > 2000000 / values[0].length() {
    raise Budget("trace cell budget")
  }
  let n = values[0].length()
  for row in values {
    if row.length() != n {
      raise Invalid("ragged trace")
    }
    for x in row {
      nonnegative(x, "trace value")
    }
  }
  { times: times.copy(), values: values.map(row => row.copy()), }
}

///|
pub fn Trace::sample_times(self : Trace) -> Array[Double] {
  self.times.copy()
}

///|
pub fn Trace::samples(self : Trace) -> Array[Array[Double]] {
  self.values.map(row => row.copy())
}

///|
pub fn Trace::last(self : Trace) -> Array[Double] {
  self.values[self.values.length() - 1].copy()
}

///|
pub fn Trace::at(
  self : Trace,
  time : Double,
) -> Array[Double] raise KineticsError {
  finite(time, "interpolation time")
  if time < self.times[0] || time > self.times[self.times.length() - 1] {
    raise Invalid("time outside trace")
  }
  let mut lo = 0
  let mut hi = self.times.length() - 1
  while lo < hi {
    let mid = lo + (hi - lo + 1) / 2
    if self.times[mid] <= time {
      lo = mid
    } else {
      hi = mid - 1
    }
  }
  if lo == self.times.length() - 1 || time == self.times[lo] {
    return self.values[lo].copy()
  }
  let weight = (time - self.times[lo]) / (self.times[lo + 1] - self.times[lo])
  Array::makei(self.values[lo].length(), i => {
    (1.0 - weight) * self.values[lo][i] + weight * self.values[lo + 1][i]
  })
}

///|
/// Trapezoid integral of the piecewise-linear sampled trajectory.
pub fn Trace::integral(self : Trace) -> Array[Double] raise KineticsError {
  let result = Array::make(self.values[0].length(), 0.0)
  for k in 1..<self.times.length() {
    let dt = self.times[k] - self.times[k - 1]
    for i in 0..<result.length() {
      result[i] += dt * (0.5 * self.values[k][i] + 0.5 * self.values[k - 1][i])
      finite(result[i], "trace integral")
    }
  }
  result
}

///|
pub fn regular_grid(
  end : Double,
  intervals : Int32,
) -> Array[Double] raise KineticsError {
  finite(end, "end time")
  if end <= 0.0 || intervals < 1 || intervals > 100000 {
    raise Invalid("regular grid bounds")
  }
  let times = Array::makei(intervals + 1, i => {
    end * (i.to_double() / intervals.to_double())
  })
  validate_grid(times)
  times
}
