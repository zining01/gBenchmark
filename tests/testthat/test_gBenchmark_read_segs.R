library(testthat)

data.table::setDTthreads(1)

test_that(desc = "read_segs total CN",
          code = {
              ## segments can be read from .rds file containing GRanges
              gr = gBenchmark:::read_segs(x = system.file("extdata", "gr.rds", package = "gBenchmark"),
                                          field = "cn")
              expect_true(length(gr) == 3)

              ## segments can be read from .rds file containing data.table
              gr = gBenchmark:::read_segs(x = system.file("extdata", "dt.rds", package = "gBenchmark"),
                                          field = "cn")
              expect_true(length(gr) == 3)

              ## segments can be read from .txt file containing data.table
              gr = gBenchmark:::read_segs(x = system.file("extdata", "dt.txt", package = "gBenchmark"),
                                          field = "cn")
              expect_true(length(gr) == 3)

              ## segments can be read from .txt file containing data.table
              ## with nonstandard column names
              gr = gBenchmark:::read_segs(x = system.file("extdata", "dt2.txt", package = "gBenchmark"),
                                          field = "cn")
              expect_true(length(gr) == 3)

              ## segments can be created from GRanges
              og.gr = readRDS(system.file("extdata", "gr.rds", package = "gBenchmark"))
              gr = gBenchmark:::read_segs(og.gr, field = "cn")
              expect_true(length(gr) == 3)

              og.dt = data.table::as.data.table(og.gr)
              gr = gBenchmark:::read_segs(og.dt, field = "cn")
              expect_true(length(gr) == 3)
          })
              
test_that(desc = "read_segs allelic CN (without cnloh)",
          code = {
              gr = gBenchmark:::read_segs(x = system.file("extdata", "agr.rds", package = "gBenchmark"),
                                          field = c("majorCN", "minorCN"),
                                          allelic = TRUE, cnloh = FALSE)
              expect_true(length(gr) == 6)
              gr = gBenchmark:::read_segs(x = system.file("extdata", "adt.rds", package = "gBenchmark"),
                                          field = c("majorCN", "minorCN"),
                                          allelic = TRUE, cnloh = FALSE)
              expect_true(length(gr) == 6)
              gr = gBenchmark:::read_segs(x = system.file("extdata", "adt.txt", package = "gBenchmark"),
                                          field = c("majorCN", "minorCN"),
                                          allelic = TRUE, cnloh = FALSE)
              expect_true(length(gr) == 6)
              gBenchmark:::read_segs(x = system.file("extdata", "adt2.txt", package = "gBenchmark"),
                                     field = c("majorCN", "minorCN"),
                                     allelic = TRUE, cnloh = FALSE)
              expect_true(length(gr) == 6)

              og.agr = readRDS(system.file("extdata", "agr.rds", package = "gBenchmark"))
              gr = gBenchmark:::read_segs(x = og.agr,
                                          field = c("majorCN", "minorCN"),
                                          allelic = TRUE,
                                          cnloh = FALSE)
              expect_true(length(gr) == 6)

              og.adt = as.data.table(og.agr)
              gr = gBenchmark:::read_segs(x = og.adt,
                                          field = c("majorCN", "minorCN"),
                                          allelic = TRUE,
                                          cnloh = FALSE)
              expect_true(length(gr) == 6)
          })

test_that(desc = "read_segs cnloh",
          code = {
              gr = gBenchmark:::read_segs(x = system.file("extdata", "agr.rds", package = "gBenchmark"),
                                          field = c("majorCN", "minorCN"),
                                          allelic = TRUE, cnloh = TRUE)
              expect_true(length(gr) == 3)
              gr = gBenchmark:::read_segs(x = system.file("extdata", "adt.rds", package = "gBenchmark"),
                                          field = c("majorCN", "minorCN"),
                                          allelic = TRUE, cnloh = TRUE)
              expect_true(length(gr) == 3)
              gr = gBenchmark:::read_segs(x = system.file("extdata", "adt.txt", package = "gBenchmark"),
                                          field = c("majorCN", "minorCN"),
                                          allelic = TRUE, cnloh = TRUE)
              expect_true(length(gr) == 3)
              gBenchmark:::read_segs(x = system.file("extdata", "adt2.txt", package = "gBenchmark"),
                                     field = c("majorCN", "minorCN"),
                                     allelic = TRUE, cnloh = TRUE)
              expect_true(length(gr) == 3)

              og.agr = readRDS(system.file("extdata", "agr.rds", package = "gBenchmark"))
              gr = gBenchmark:::read_segs(x = og.agr,
                                          field = c("majorCN", "minorCN"),
                                          allelic = TRUE,
                                          cnloh = TRUE)
              expect_true(length(gr) == 3)

              og.adt = as.data.table(og.agr)
              gr = gBenchmark:::read_segs(x = og.adt,
                                          field = c("majorCN", "minorCN"),
                                          allelic = TRUE,
                                          cnloh = TRUE)
              expect_true(length(gr) == 3)
          })

