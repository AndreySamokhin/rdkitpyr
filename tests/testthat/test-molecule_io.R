test_that("ConvertToSmiles(), convert InChI to SMILES", {
  expect_identical(
    ConvertToSmiles(test_compounds$inchi),
    test_compounds$smiles
  )
  expect_identical(ConvertToSmiles("InChI=1S/CH4/h1H4"), "C")
  expect_identical(ConvertToSmiles("InChI=1S/CH4/h1H4_???"), NA_character_)
  expect_identical(ConvertToSmiles(NA_character_), NA_character_)
})


test_that("ConvertToSmiles(), convert SMILES to SMILES", {
  expect_identical(
    ConvertToSmiles(c("OCCC", "CC(=O)C", "c1cc(O)ccc1", "c1c(C)c(O)ccc1")),
    c("CCCO", "CC(C)=O", "Oc1ccccc1", "Cc1ccccc1O")
  )
  expect_identical(
    ConvertToSmiles(c("F\\C=C\\F", "O[C@H](F)C"), isomeric = FALSE),
    c("FC=CF", "CC(O)F")
  )
  expect_identical(
    ConvertToSmiles(c("CCCCC", "c1ccoc1", "c1ccccc1"), kekule = TRUE),
    c("CCCCC", "C1=COC=C1", "C1=CC=CC=C1")
  )
  expect_identical(
    ConvertToSmiles(c("CCCO", "CC(C)=O", "c1ccccc1"), explicit_bonds = TRUE),
    c("C-C-C-O", "C-C(-C)=O", "c1:c:c:c:c:c:1")
  )
  expect_identical(
    ConvertToSmiles(c("CCCO", "c1ccccc1"), explicit_hydrogens = TRUE),
    c("[CH3][CH2][CH2][OH]", "[cH]1[cH][cH][cH][cH][cH]1")
  )
  expect_identical(ConvertToSmiles("C"), "C")
  expect_identical(ConvertToSmiles("C_???"), NA_character_)
})


test_that("ConvertToInchi(), convert InChI to InChI", {
  expect_identical(ConvertToInchi("InChI=1S/CH4/h1H4"), "InChI=1S/CH4/h1H4")
  expect_identical(ConvertToInchi("InChI=1S/CH4/h1H4_???"), NA_character_)
})


test_that("ConvertToInchi(), convert SMILES to InChI", {
  expect_identical(
    ConvertToInchi(test_compounds$smiles),
    test_compounds$inchi
  )
  expect_identical(ConvertToInchi("C"), "InChI=1S/CH4/h1H4")
  expect_identical(ConvertToInchi("C_???"), NA_character_)
  expect_identical(ConvertToInchi(NA_character_), NA_character_)
})


test_that("ConvertToInchikey(), convert InChI to InChIKey", {
  expect_identical(
    ConvertToInchikey(test_compounds$inchi),
    test_compounds$inchikey
  )
  expect_identical(
    ConvertToInchikey("InChI=1S/CH4/h1H4"),
    "VNWKTOKETHGBQD-UHFFFAOYSA-N"
  )
  expect_identical(ConvertToInchikey("InChI=1S/CH4/h1H4_???"), NA_character_)
  expect_identical(ConvertToInchikey(NA_character_), NA_character_)
})


test_that("ConvertToInchikey(), convert SMILES to InChIKey", {
  expect_identical(
    ConvertToInchikey(test_compounds$smiles),
    test_compounds$inchikey
  )
  expect_identical(ConvertToInchikey("C"), "VNWKTOKETHGBQD-UHFFFAOYSA-N")
  expect_identical(ConvertToInchikey("C_???"), NA_character_)
})


