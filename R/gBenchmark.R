#' @import data.table
#' @import gUtils
#' @import rtracklayer
#'
#' @importFrom gGnome gG
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom IRanges IRanges

#' @name benchmark_cn
#' @title benchmark_cn
#'
#' @description
#' compare copy number (CN) predictions
#'
#' @details
#' checks consistency between total CN, allelic CN, or CNLOH predictions
#'
#' inputs can be one of the follow objects, or a file containing one of the following objects:
#' - bigwig
#' - data.table coercible to GRanges
#' - GRanges
#' - gGnome gGraph
#'
#' note that the comparison for cnloh is non-symmetric and y is assumed to be the ground truth
#' 
#' @param x character vector containing a file path, data.table, GRanges, or gGraph
#' @param y character vector containing a file path, data.table, GRanges, or gGraph
#' @param x.field (character) field(s) in x containing CN data (default "cn")
#' @param y.field (character) field(s) in y containing CN data (default "cn")
#' @param allelic (logical) allele analysis? default FALSE
#' @param cnloh (logical) cnloh analysis? default FALSE
#' @param tile.width (numeric) tile width for averaging copy number estimates (default 10 Kbp)
#' @param min.width (numeric) minimum width segment to retain (default 0, retains all segments)
#' @param genome (character) one of hg19 or h38, default hg19
#' @param std.only (character) only consider standard chromosomes? e.g. 1-22, X, Y. default TRUE
#' @param verbose (logical) default FALSE
#'
#' @return data.table with benchmarking results
#'
#' @export
benchmark_cn = function(x, y,
                        x.field = "cn",
                        y.field = "cn",
                        allelic = FALSE,
                        cnloh = FALSE,
                        tile.width = 1e4,
                        min.width = 0,
                        genome = "hg19",
                        std.only = TRUE,
                        verbose = FALSE)
{
    ## get chromosome sizes
    if (genome != "hg19" & genome != "hg38")
    {
        stop("genome must be one of hg19 or hg38")
    }
    chrom.sizes.dt = data.table::fread(file = system.file("extdata", paste0(genome, ".chrom.sizes"),
                                                          package = "gBenchmark"),
                                       header = FALSE,
                                       col.names = c("seqnames", "size"))
    if (std.only)
    {
        chrom.sizes.dt = chrom.sizes.dt[grepl("^(chr)*[0-9XY]+$", seqnames),]
    }
    chrom.sizes = chrom.sizes.dt[, size]
    names(chrom.sizes) = chrom.sizes.dt[, seqnames]

    ## tile eligible regions
    tiles.gr = gUtils::gr.tile(gUtils::si2gr(si = chrom.sizes), width = tile.width)

    ## grab input data
    gr1 = read_segs(x, field = x.field, allelic = allelic, cnloh = cnloh)
    gr2 = read_segs(y, field = y.field, allelic = allelic, cnloh = cnloh)

    ## overlap with tiles
    if (cnloh | (!allelic))
    {
        GenomicRanges::mcols(tiles.gr)[, "score1"] = gUtils::gr.val(query = tiles.gr,
                                                                    target = gr1,
                                                                    val = "score",
                                                                    na.rm = TRUE)$score
        GenomicRanges::mcols(tiles.gr)[, "score2"] = gUtils::gr.val(query = tiles.gr,
                                                                    target = gr2,
                                                                    val = "score",
                                                                    na.rm = TRUE)$score
    }
    else
    {
        tiles2.gr = tiles.gr[, c()]
        GenomicRanges::mcols(tiles.gr)[, "score1"] = gUtils::gr.val(query = tiles.gr,
                                                                    target = gr1 %Q% (allele == "major"),
                                                                    val = "score",
                                                                    na.rm = TRUE)$score
        GenomicRanges::mcols(tiles.gr)[, "score2"] = gUtils::gr.val(query = tiles.gr,
                                                                    target = gr2 %Q% (allele == "major"),
                                                                    val = "score",
                                                                    na.rm = TRUE)$score
        GenomicRanges::mcols(tiles2.gr)[, "score1"] = gUtils::gr.val(query = tiles2.gr,
                                                                     target = gr1 %Q% (allele == "minor"),
                                                                     val = "score",
                                                                     na.rm = TRUE)$score
        GenomicRanges::mcols(tiles2.gr)[, "score2"] = gUtils::gr.val(query = tiles2.gr,
                                                                     target = gr2 %Q% (allele == "minor"),
                                                                     val = "score",
                                                                     na.rm = TRUE)$score
        GenomicRanges::mcols(tiles.gr)[, "allele"] = "major"
        GenomicRanges::mcols(tiles2.gr)[, "allele"] = "minor"
        tiles.gr = c(tiles.gr, tiles2.gr)
    }

    ## rmse/pearson/spearman correlation if integer CN
    ## precision/recall/F1 if cnloh
    if (!cnloh)
    {
        res = data.table(pearson.cn = stats::cor(x = GenomicRanges::mcols(tiles.gr)[, "score1"],
                                                 y = GenomicRanges::mcols(tiles.gr)[, "score2"],
                                                 use = "na.or.complete",
                                                 method = "pearson"),
                         spearman.cn = stats::cor(x = GenomicRanges::mcols(tiles.gr)[, "score1"],
                                                  y = GenomicRanges::mcols(tiles.gr)[, "score2"],
                                                  use = "na.or.complete",
                                                  method = "spearman"),
                         rmse = sqrt(mean((GenomicRanges::mcols(tiles.gr)[, "score1"] -
                                           GenomicRanges::mcols(tiles.gr)[, "score2"])^2,
                                          na.rm = TRUE)))
    }
    else
    {
        tp = sum(GenomicRanges::mcols(tiles.gr)[, "score1"] > 0 &
                 GenomicRanges::mcols(tiles.gr)[, "score2"] > 0,
                 na.rm = TRUE)
        fn = sum(GenomicRanges::mcols(tiles.gr)[, "score1"] == 0 &
                 GenomicRanges::mcols(tiles.gr)[, "score2"] > 0,
                 na.rm = TRUE)
        fp = sum(GenomicRanges::mcols(tiles.gr)[, "score1"] > 0 &
                 GenomicRanges::mcols(tiles.gr)[, "score2"] == 0,
                 na.rm = TRUE)
        tn = sum(GenomicRanges::mcols(tiles.gr)[, "score1"] == 0 &
                 GenomicRanges::mcols(tiles.gr)[, "score2"] == 0,
                 na.rm = TRUE)
        res = data.table(precision.cnloh = ifelse(tp > 0, tp / (tp + fp), 0),
                         recall.cnloh = ifelse(tp > 0, tp / (tp + fn), 0))
        res[, f1.cnloh := ifelse(precision.cnloh > 0 & recall.cnloh > 0,
                                 2 * precision.cnloh * recall.cnloh / (precision.cnloh + recall.cnloh),
                                 0)]
    }
    
    return(res)
}

#' @name read_segs
#' @title read_segs
#'
#' @description
#' read segments and copy number from input file
#'
#' @param x character vector containing a file path, data.table, GRanges, bigWig, or gGraph
#' @param field field(s) in x containing copy number data (default 'cn')
#' @param allelic (logical) allelic analysis? default FALSE
#' @param cnloh (logical) CNLOH analysis? default FALSE
#' @param xy.sub (logical) if '23' and '24' are provided as seqnames, replace them with X and Y? default TRUE
#'
#' @return GRanges with field "score"
#' score represents copy number if "cnloh" is TRUE
#' additional field "allele" with values "major" and "minor" if "allelic" is TRUE
read_segs = function(x,
                     field = "cn",
                     allelic = FALSE,
                     cnloh = FALSE,
                     xy.sub = TRUE)
{
    if (inherits(x, "character"))
    {
        if (!check_file(x)) { stop("invalid file supplied: ", x) }
        if (grepl(".rds$", x))
        {
            x = readRDS(x)
        }
        else if (grepl("(.csv$)|(.txt$)|(.tsv$)", x))
        {
            x = data.table::fread(x)
        }
        else if (grepl(".bw$", x))
        {
            x = rtracklayer::import.bw(con = x)
        }
        else
        {
            stop("invalid file type supplied: ", x)
        }
    }

    ## check if x is a gGraph... if so, use just the segments
    if (inherits(x, "list"))
    {
        x = gGnome::gG(jabba = x)
    }

    if (inherits(x, "gGraph"))
    {
        x = gUtils::gr.stripstrand(x$nodes$gr)
    }

    if (inherits(x, "data.table"))
    {
        x = copy(x)
        ## make sure field is in the column names
        if (!all(field %in% names(x))) { stop(field, " is not in column names!") }
        ## make sure column names conform to something that can be coerced to GRanges
        sn = grep("^([Cc]hr*)|(seqnames)", colnames(x), value = T)
        st = st = grep("^[Ss]tart", colnames(x), value = T)
        ed = grep("^[Ee]nd", colnames(x), value = T)
        if (length(sn) != 1 | length(st) != 1 | length(ed) != 1)
        {
            stop("supplied data has invalid column names!")
        }
        ## coerce to expected column names
        setnames(x, old = c(sn, st, ed), new = c("seqnames", "start", "end"))
        if (xy.sub)
        {
            dt[seqnames == "23", seqnames := "X"]
            dt[seqnames == "24", seqnames := "Y"]
        }

        if (allelic | cnloh)
        {
            setnames(x, field, c("major_cn", "minor_cn"))
            if (!cnloh)
            {
                x = rbind(x[, .(seqnames, start, end, strand = "*", score = major_cn, allele = "major")],
                          x[, .(seqnames, start, end, strand = "*", score = minor_cn, allele = "minor")])
                out = GRanges(seqnames = x[, seqnames],
                              ranges = IRanges(start = x[, start],
                                               end = x[, end]),
                              strand = "*",
                              score = x[, score],
                              allele = x[, allele])
            }
            else
            {
                out = GRanges(seqnames = x[, seqnames],
                              ranges = IRanges(start = x[, start],
                                               end = x[, end]),
                              strand = "*",
                              score = x[, pmax(major_cn, minor_cn) == 2 & pmin(major_cn, minor_cn) == 0])
            }
        }
        else
        {
            out = GRanges(seqnames = x[, seqnames],
                        ranges = IRanges(start = x[, start],
                                         end = x[, end]),
                        strand = "*",
                        score = x[, base::get(field)])
        }
    } else if (inherits(x, "GRanges"))
    {
        if (!all(field %in% names(values(x)))) { stop(field, " contains invalid values!") }
        if (allelic | cnloh)
        {
            if (!cnloh)
            {
                major.gr = x[, c()]
                GenomicRanges::mcols(major.gr)[, "score"] = GenomicRanges::mcols(x)[, field[1]]
                GenomicRanges::mcols(major.gr)[, "allele"] = "major"
                minor.gr = x[, c()]
                GenomicRanges::mcols(minor.gr)[, "score"] = GenomicRanges::mcols(x)[, field[2]]
                GenomicRanges::mcols(minor.gr)[, "allele"] = "minor"
                out = c(major.gr, minor.gr)
            } else
            {
                out = x[, c()]
                GenomicRanges::mcols(out)[, "score"] = pmax(GenomicRanges::mcols(x)[, field[1]],
                                                            GenomicRanges::mcols(x)[, field[2]]) == 2 &
                    pmin(GenomicRanges::mcols(x)[, field[1]],
                         GenomicRanges::mcols(x)[, field[2]]) == 0
                GenomicRanges::mcols(out)[, "score"] = as.numeric(GenomicRanges::mcols(out)[, "score"])
            }
        }
        else
        {
            out = x[, c()]
            GenomicRanges::mcols(out)[, "score"] = GenomicRanges::mcols(x)[, field]
        }
    } else
    {
        stop("invalid data type!")
    }

    ## at this point x should be a GRanges
    return(out)
}

#' @name check_file
#' @title check_file
#'
#' @description
#' Check whether supplied file exists and is nonempty
#'
#' @param fn
#'
#' @return logical, TRUE if file exists and is nonempty
check_file = function(fn = NULL) {
    if (is.null(fn)) {
        return(FALSE)
    }

    if (!file.exists(fn)) {
        return(FALSE)
    }

    if (!file.info(fn)$size) {
        return(FALSE)
    }

    return(TRUE)
}
