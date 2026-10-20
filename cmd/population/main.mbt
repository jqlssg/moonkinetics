///|
fn main raise {
  // Immigration-death population, concentration birth b and per-capita death d.
  let b = 4.0
  let d = 0.5
  let net = @k.Network::new(["population"], [
    @k.Reaction::new([0], [1], b),
    @k.Reaction::new([1], [0], d),
  ])
  let grid = @k.regular_grid(2.0, 10)
  let ensemble = net.ensemble([10], grid, 2000, seed=731, max_events=200)
  let means = ensemble.means()
  let final_mean = means.last()[0]
  let variance = ensemble.covariances()[10][0][0]
  let q = @math.exp(-d * 2.0)
  let expected = 10.0 * q + b / d * (1.0 - q)
  let expected_variance = 10.0 * q * (1.0 - q) + b / d * (1.0 - q)
  let se = ensemble.standard_errors()[10][0]
  assert_true((final_mean - expected).abs() < 6.0 * se)
  assert_true((variance - expected_variance).abs() < 1.0)
  let path = net.simulate([10], 2.0, seed=42)
  let likelihood = net.log_likelihood(path)
  let score = net.path_score(path)
  assert_true(!likelihood.is_nan() && !likelihood.is_inf())
  assert_eq(score.length(), 2)
  let exposure = net.integrated_propensities(path)
  assert_eq(exposure[0], 8.0)
  println(
    "population: runs=\{ensemble.runs()} mean=\{final_mean} analytic=\{expected} variance=\{variance} analytic_variance=\{expected_variance} standard_error=\{se}",
  )
  println(
    "path_log_likelihood=\{likelihood} birth_score=\{score[0]} death_score=\{score[1]}",
  )
  println(means.to_csv())
}
