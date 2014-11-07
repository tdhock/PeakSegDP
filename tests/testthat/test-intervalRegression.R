context("Interval regression")

test_that("target interval in the middle, uniform cost", {
  cost <- c(2,2,1,1,1,
            0,0,1,1,1,
            0,0,0,1)
  size <- rep(1, length(cost))
  largest <- largestContinuousMinimum(cost, size)
  expect_equal(largest$start, 11)
  expect_equal(largest$end, 13)

  cost <- c(2,2,1,1,1,
            0,0,0,1,1,
            1,0,0,1)
  size <- rep(1, length(cost))
  largest <- largestContinuousMinimum(cost, size)
  expect_equal(largest$start, 6)
  expect_equal(largest$end, 8)
})

test_that("target interval in the middle, non-uniform cost", {
  cost <- c(2,2,1,1,1,
            0,1,1,1,1,
            0,0,0,1)
  size <- rep(1, length(cost))
  largest <- largestContinuousMinimum(cost, size)
  expect_equal(largest$start, 11)
  expect_equal(largest$end, 13)

  size[6] <- 4
  largest <- largestContinuousMinimum(cost, size)
  expect_equal(largest$start, 6)
  expect_equal(largest$end, 6)
})

test_that("target interval at the start", {
  cost <- c(0,0,0,0,0,0,0,0,1,1,1,1,0,1,1)
  size <- rep(1, length(cost))
  largest <- largestContinuousMinimum(cost, size)
  expect_equal(largest$start, 1)
  expect_equal(largest$end, 8)
})

test_that("target interval at the end", {
  cost <- c(1,1,0,0,1,1,1,0,0,0,0,0,0,0,0)
  size <- rep(1, length(cost))
  largest <- largestContinuousMinimum(cost, size)
  expect_equal(largest$start, 8)
  expect_equal(largest$end, 15)
})
