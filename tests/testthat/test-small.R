context("small data set")

library(PeakSegDP)
data(H3K36me3.AM.immune.19)

## profiles <- do.call(rbind, H3K36me3.AM.immune.19)
## ggplot()+
##   geom_step(aes(chromStart/1e3, count), data=profiles)+
##   facet_grid(cell.type + sample.id ~ .)+
##   theme_bw()+
##   theme(panel.margin=grid::unit(0, "cm"))

test_that("infeasible model generates warning", {
  one <- H3K36me3.AM.immune.19[["McGill0101"]]
  expect_warning({
    fit <- PeakSegDP(one, maxPeaks=1L)
  }, "infeasible model")
})


