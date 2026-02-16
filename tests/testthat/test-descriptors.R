test_that("CalculateExactMass()", {
  expect_equal(
    CalculateExactMass(test_compounds$smiles),
    test_compounds$exact_mass
  )
  expect_identical(CalculateExactMass("C???"), NA_real_)
})


test_that("CalculateMolecularWeight()", {
  expect_equal(
    CalculateMolecularWeight(test_compounds$smiles),
    test_compounds$molecular_weight,
    tolerance = 0.001
  )
  expect_identical(CalculateMolecularWeight("C???"), NA_real_)
})

