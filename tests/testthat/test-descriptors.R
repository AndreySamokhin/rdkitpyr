test_that("CalculateExactMass()", {
  expect_equal(
    CalculateExactMass(test_compounds$smiles),
    test_compounds$exact_mass
  )
  expect_identical(CalculateExactMass("C???"), NA_real_)
})

