test_that("symmetry wrt favor_positive", {
  mpm <- money_priming_sub()

  mpm <- mpm |> mutate(yi_flipped = -yi)

  # using only 5 studies to speed up tests; too-long tests fail CRAN checks
  res1 <- phacking_meta(mpm$yi[1:5], mpm$vi[1:5], parallelize = FALSE)
  res2 <- phacking_meta(mpm$yi_flipped[1:5], mpm$vi[1:5], parallelize = FALSE,
                        favor_positive = FALSE)

  expect_equal(res1$stats$mode, res2$stats$mode, tolerance = 0.01)
})
