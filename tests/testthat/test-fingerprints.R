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


test_that("CalculateRdkitFingerprints()", {
  skip_on_cran() # RDKit (Python module) is unlikely available on CRAN machines
  SkipIfNoRdkit()
  fps <- CalculateRdkitFingerprints(c("C", "C???"), fp_size = 128L)
  expect_true(is.integer(fps))
  expect_true(is.matrix(fps))
  expect_identical(ncol(fps), 128L)
  expect_identical(nrow(fps), 2L)
  expect_true(all(fps[1L, ] %in% c(0L, 1L)))
  expect_true(all(is.na(fps[2L, ])))
  expect_identical(attr(fps, "valid"), !is.na(fps[, 1L]))
  expect_no_error(CalculateRdkitFingerprints(
    c("C", "CC", "CCC"),
    min_path = 1.0,
    max_path = 7.0,
    fp_size = 128.0,
    n_bits_per_hash = 2.0,
    min_size = 128.0
  ))
})


test_that("CalculateMorganFingerprints()", {
  skip_on_cran() # RDKit (Python module) is unlikely available on CRAN machines
  SkipIfNoRdkit()
  fps <- CalculateMorganFingerprints(c("C", "C???"), fp_size = 128L)
  expect_true(is.integer(fps))
  expect_true(is.matrix(fps))
  expect_identical(ncol(fps), 128L)
  expect_identical(nrow(fps), 2L)
  expect_true(all(fps[1L, ] %in% c(0L, 1L)))
  expect_true(all(is.na(fps[2L, ])))
  expect_identical(attr(fps, "valid"), !is.na(fps[, 1L]))
  expect_no_error(CalculateMorganFingerprints(
    c("C", "CC", "CCC"),
    radius = 3.0,
    fp_size = 128.0
  ))
})

