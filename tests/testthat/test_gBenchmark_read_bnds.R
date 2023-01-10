library(testthat)

data.table::setDTthreads(1)

test_that(desc = "read breakends",
          code = {
              ## check that breakends can be read from gGraph
              res = gBenchmark:::read_bnds(x = system.file("extdata", "gg.rds", package = "gBenchmark"))
              expect_true(length(res) == 10)

              ## check that strand information is removed if desired
              res = gBenchmark:::read_bnds(x = system.file("extdata", "gg.rds", package = "gBenchmark"),
                                           ignore.strand = TRUE)
              expect_true(length(res) == 10)
              expect_true(all(as.character(strand(res)) == "*"))

              ## check that breakends can be read from segments
              res = gBenchmark:::read_bnds(x = system.file("extdata", "gr.rds", package = "gBenchmark"), field = "cn")
              expect_true(length(res) == 2)
              
              ## check that breakends can be read from segments provided as data.table
              res = gBenchmark:::read_bnds(x = system.file("extdata", "dt.rds", package = "gBenchmark"), field = "cn")
              expect_true(length(res) == 2)

              ## check that breakends can be read from Junctions
              res = gBenchmark:::read_bnds(x = system.file("extdata", "jj.rds", package = "gBenchmark"))
              expect_true(length(res) == 10)

              ## check that breakends can be read from GRangesList
              res = gBenchmark:::read_bnds(x = system.file("extdata", "grl.rds", package = "gBenchmark"))
              expect_true(length(res) == 10)

              ## check that breakends can be read from
              res = gBenchmark:::read_bnds(system.file("extdata", "svaba.subset.vcf", package = "gBenchmark"))
              expect_true(length(res) == 40)
          })
