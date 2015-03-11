context("trivial examples")

pseg3 <- function(...){
  y <- c(...)
  w <- rep(1, length(y))
  ## TODO: when we write a real solver it should replace cDPA, which
  ## does not pass all these tests.
  fit <- cDPA(y, rep(1, length(y)), 3)
  ends <- getPath(fit)
  e <- ends[3,]
  if(any(is.na(e))){
    stop("no feasible model with 3 segments")
  }
  e
}
  
test_that("infeasible models", {
  expect_error({
    pseg3(3, 2, 1)
  }, "no feasible model with 3 segments")
  expect_error({
    pseg3(1, 2, 3)
  }, "no feasible model with 3 segments")
})

test_that("feasible models", {
  ends <- pseg3(1, 3, 2)
  expect_equal(ends, 1:3)
  ends <- pseg3(2, 3, 1)
  expect_equal(ends, 1:3)
})

## A real solver for the Peaks problem would be able to pass this
## test. The DP solver in this package constraints too much, so it
## does not pass.

## test_that("feasible model 4 data points", {
##   ends <- pseg3(97, 101, 104, 103)
##   expect_equal(ends, c(2, 3, 4))
## })
