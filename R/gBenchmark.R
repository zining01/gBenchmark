#' @import data.table
#' @import gUtils
#' @import rtracklayer
#'
#' @importFrom gGnome gG
#' @importFrom GenomicRanges GRanges, mcols
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
#' @param x character vector containing a file path, data.table, GRanges, or gGraph
#' @param y character vector containing a file path, data.table, GRanges, or gGraph
#' @param x.field (character) field(s) in x containing CN data (default "cn")
#' @param y.field (character) field(s) in y containing CN data (default "cn")
#' @param allelic (logical) allele analysis? default FALSE
#' @param cnloh (logical) cnloh analysis? default FALSE
#' @param tile.width (numeric) tile width for averaging copy number estimates (default 10 Kbp)
#' @param min.width (numeric) minimum width segment to retain (default 0, retains all segments)
#' @param verbose (logical) default FALSE
benchmark_cn = function(x, y,
                        x.field = "cn",
                        y.field = "cn",
                        allelic = FALSE,
                        cnloh = FALSE,
                        tile.width = 1e4,
                        min.width = 0,
                        verbose)
{
    
    return()
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
        x = x$nodes$gr
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
            }
        }
        else
        {
            out = x[, field]
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
