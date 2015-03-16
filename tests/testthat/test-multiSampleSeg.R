context("multiSampleSeg")

count <-
  as.integer(c(rep(c(0, 1), each=5),
               rep(c(10, 11), each=10),
               rep(c(0, 1), each=5)))

test_that("with 1 sample we refine peaks" , {
  chromEnd <- seq_along(count)
  profiles <-
    data.frame(chromStart=chromEnd-1L,
               chromEnd,
               sample.id="sample1",
               count)
  peak <- multiSampleSegHeuristic(profiles)
  expect_equal(peak$chromStart, 10)
  expect_equal(peak$chromEnd, 30)
})
