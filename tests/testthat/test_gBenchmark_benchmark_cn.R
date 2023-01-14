library(testthat)

test_that(desc = "check that coverage can be simulated",
          code = {
              res = gBenchmark::benchmark_cn(x = system.file("extdata", "ovcar.gr.rds", package = "gBenchmark"), y = system.file("extdata", "ovcar.random.gr.rds", package = "gBenchmark"))
              expect_true(res$pearson.cn > 0)
              expect_true(res$spearman.cn > 0)
              expect_true(res$rmse > 0)

              res = gBenchmark::benchmark_cn(x = system.file("extdata", "ovcar.gr.rds", package = "gBenchmark"), y = system.file("extdata", "bias.gr.rds", package = "gBenchmark"), y.field = "background")
              expect_true(res$pearson.cn < 0.1)
              expect_true(res$spearman.cn < 0.1)
              expect_true(res$rmse > 1)
          })
