test_that("CalculateAllDescriptors()", {
  desc <- CalculateAllDescriptors(test_compounds$smiles)
  expect_true(is.data.frame(desc))
  expect_identical(nrow(desc), 10L)
  expect_true(ncol(desc) > 100L)
  expect_equal(desc$MolWt, CalculateMolecularWeight(test_compounds$smiles))
  expect_equal(desc$ExactMolWt, CalculateExactMass(test_compounds$smiles))

  expect_identical(
    as.integer(CalculateAllDescriptors(c("C???", "C"))$MolWt),
    c(NA_integer_, 16L)
  )
  expect_error(
    CalculateAllDescriptors(c("C???", "CC???", "CCC???")),
    "All input molecules are invalid."
  )
})


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

