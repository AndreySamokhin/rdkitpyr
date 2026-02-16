test_that("CalculateMaccsFingerprints()", {
  skip_on_cran() # RDKit (Python module) is unlikely available on CRAN machines
  SkipIfNoRdkit()
  fps <- CalculateMaccsFingerprints(c("C", "C???"))
  expect_true(is.integer(fps))
  expect_true(is.matrix(fps))
  expect_identical(ncol(fps), 167L)
  expect_identical(nrow(fps), 2L)
  expect_true(all(fps[1L, ] %in% c(0L, 1L)))
  expect_true(all(is.na(fps[2L, ])))
  expect_identical(attr(fps, "valid"), !is.na(fps[, 1L]))
})

