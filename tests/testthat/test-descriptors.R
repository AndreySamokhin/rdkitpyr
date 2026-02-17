test_that("CalculateAllDescriptors()", {
  desc <- CalculateAllDescriptors(test_compounds$smiles)
  expect_true(is.data.frame(desc))
  expect_identical(nrow(desc), 10L)
  expect_true(ncol(desc) > 100L)
  expect_equal(desc$MolWt, CalculateMolecularWeight(test_compounds$smiles))
  expect_equal(desc$ExactMolWt, CalculateExactMass(test_compounds$smiles))

  desc <- CalculateAllDescriptors(c("C???", "C"))
  expect_identical(desc$MolWt[[1L]], NA_real_)
  expect_identical(as.integer(desc$MolWt[[2L]]), 16L)

  desc <- CalculateAllDescriptors(c("C???", "CC???", "CCC???"))
  expect_true(is.data.frame(desc))
  expect_true(all(is.na(desc)))
  expect_true(all(c("MolWt", "MolLogP", "TPSA", "MolMR") %in% names(desc)))
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

