///|
test "weighted mass invariant for dimer and empty network" {
  let n = Network::new(["A", "B"], [Reaction::new([2, 0], [0, 1], 3.0)])
  assert_eq(n.conservation_residual([1.0, 2.0]), [0.0])
  assert_eq(n.conservation_basis(), [[0.5, 1.0]])
  assert_eq(conserved_value([1.0, 2.0], [4.0, 1.0]), 6.0)
  let e = Network::new(["X", "Y"], [])
  assert_eq(e.conservation_basis(), [[1.0, 0.0], [0.0, 1.0]])
}
