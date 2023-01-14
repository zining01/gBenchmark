library(testthat)

test_that(desc = "check that coverage can be simulated",
          code = {
              ## no bias is used
              gr = gBenchmark::simulate_coverage(gr = system.file("extdata", "ovcar.gr.rds", package = "gBenchmark"))
              expect_true(any(!is.na(gr$cov)))
              expect_true(all(gr$background == 1))

              ## check that a nonzero bias can be used
              gr = gBenchmark::simulate_coverage(gr = system.file("extdata", "ovcar.gr.rds", package = "gBenchmark"), bias = system.file("extdata", "bias.gr.rds", package = "gBenchmark"))
              expect_true(any(gr$background != 1, na.rm = TRUE))
          })
