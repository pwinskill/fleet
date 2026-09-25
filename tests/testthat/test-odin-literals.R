test_that("every literal in the odin source prints the same on every platform", {
  # odin2 prints a fractional literal with deparse(control = "digits17"), whose
  # output differs between platforms unless the number is a short binary
  # fraction: 1e-12 prints as 1e-12 on Windows arm64 and 9.9999999999999998e-13
  # on Linux, so the regenerated C++ differs and CI's generated-code check fails
  # although nothing changed. Whole numbers print with format(), the same on all
  # of them, so such a constant is written as a quotient, 1 / 1e12.
  src <- system.file("odin", "malaria_daily.R", package = "fleet")
  pd <- utils::getParseData(parse(src, keep.source = TRUE))
  x <- suppressWarnings(as.numeric(pd$text[pd$token == "NUM_CONST"]))
  frac <- unique(x[is.finite(x) & x != round(x)])
  unstable <- frac[frac * 2^20 != round(frac * 2^20)]
  expect_equal(unstable, numeric(0),
               info = "write these as quotients of whole numbers, e.g. 1 / 1e12")
})
