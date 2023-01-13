#' @import data.table
#' @import gUtils
#' @import rtracklayer
#'
#' @importFrom gGnome gG jJ
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom IRanges IRanges
#' @importFrom stats rpois rbinom

#' @name simulate_coverage
#' @title simulate_coverage
#'
#' @description
#' simulate binned read counts given expected depth
#'
#' @details
#' Given a genomic segments with estimated tumor copy number (CN) and genomic bin size, simulate read counts per genomic bin.
#' This can be done for either a diploid genome (diploid = TRUE) or haploid genome (haploid = TRUE).
#'
#' Optionally, a multiplicative bias (representing replication timing, batch effect, etc.) can be supplied to make the coverage appear more realistic.
#' This should be supplied as a GRanges as a numeric vector in field "background" (or some other field specified in bias.field)
#' 
#' The expected coverage and read size can be tuned by setting basecov and readsize, respectively.
#'
#' By default, read counts are simulated from a Poisson distribution, but if overdispersion (numeric) is supplied will simulate from a Negative Binomial distribution.
#'
#' @param gr (GRanges, gGraph, or file containing one of these things) segment copy number, with numeric cn supplied in field "cn"
#' @param field (character) field containing copy number in gr, default cn
#' @param bins (GRanges) genomic bins
#' @param binsize (numeric) width of genomic bins (bp) for simulating read counts. ignored if bins if supplied. default 1e3.
#' @param diploid (logical) simulate diploid genome? default TRUE
#' @param bias (GRanges) mulitplicative bias for binned means
#' @param bias.field (character) column in bias GRanges metadata containing multiplicative bias. default "background"
#' @param basecov (numeric) expected coverage depth. default 60.
#' @param readsize (numeric) expected read size (bp). default 150.
#' @param purity (numeric) sample purity, between zero and one. default 0.95.
#' @param poisson (numeric) simulate actual integer reads? if NA will just output relative copy number.
#' @param overdispersion (numeric) if supplied, should be > 1 and will be applied to simulate from a negative binomial distribution. default NA.
#' @param normalize (logical) unit normalize resulting reads by dividing by mean? default TRUE
#'
#' @return GRanges with coverage per genomic bin in field "cov"
#'
#' @export
simulate_coverage = function(gr,
                             field = "cn",
                             bins = NULL,
                             binsize = 1000,
                             diploid = TRUE,
                             bias = NULL,
                             bias.field = "background",
                             basecov = 60,
                             readsize = 150,
                             purity = 0.95,
                             poisson = TRUE,
                             overdispersion = NA,
                             normalize = TRUE,
                             genome = "hg19")
{
    gr = read_segs(gr, field = field, allelic = FALSE, cnloh = FALSE)
    chrom.sizes = grab_chrom_sizes(genome)

    ## generate tiled genome
    if (is.null(bins) || !inherits(bins, "GRanges"))
    {
        bins = gUtils::gr.tile(gUtils::si2gr(chrom.sizes), width = binsize)
    }

    ## compute expected bin coverage
    bincov = (basecov * (readsize + binsize - 1)) / readsize

    ## get expected normal cn
    ncn = 1
    if (diploid) { ncn = 2 }

    ## compute relative CN
    GenomicRanges::mcols(bins)[, "cn"] = gUtils::gr.val(query = bins,
                                                        target = gr,
                                                        val = "score",
                                                        mean = TRUE,
                                                        na.rm = TRUE)$score

    GenomicRanges::mcols(bins)[, "cn.rel"] = (purity * GenomicRanges::mcols(bins)[, "cn"] +
                                              ncn * (1 - purity)) /
        (purity * mean(GenomicRanges::mcols(bins)[, "cn"], na.rm = TRUE) +
         ncn * (1 - purity))

    ## add multiplcative bias, if supplied
    GenomicRanges::mcols(bins)[, "background"] = 1
    if (!is.null(bias) && inherits(bias, "GRanges"))
    {
        GenomicRanges::mcols(bias)[, "background"] = GenomicRanges::mcols(bias)[, bias.field]
        GenomicRanges::mcols(bins)[, "background"] = gUtils::gr.val(query = bins,
                                                                    target = bias,
                                                                    val = "background",
                                                                    mean = TRUE,
                                                                    na.rm = TRUE)$background
    }

    ## simulate from Poisson distribution
    GenomicRanges::mcols(bins)[, "cov"] = GenomicRanges::mcols(bins)[, "cn.rel"] * GenomicRanges::mcols(bins)[, "background"] * bincov
    if (poisson)
    {
        if (is.na(overdispersion))
        {
            GenomicRanges::mcols(bins)[, "cov"] = stats::rpois(n = length(bins), lambda = GenomicRanges::mcols(bins)[, "cov"])
        }
        else
        {
            GenomicRanges::mcols(bins)[, "cov"] = stats::rbinom(n = length(bins), mu =  GenomicRanges::mcols(bins)[, "cov"], size = 1 / overdispersion)
        }
    }

    ## normalize if desired
    if (normalize)
    {
        GenomicRanges::mcols(bins)[, "cov"] = GenomicRanges::mcols(bins)[, "cov"]  / mean(GenomicRanges::mcols(bins)[, "cov"], na.rm = TRUE)
    }
    return(bins)
}
                             
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
#' the "vanilla" analysis is comparison to total copy number and will return a data.table with Pearson correlation, Spearman correlation, and root mean square error (RMSE) computed over all genomic tiles
#'
#' one variation of this analysis is comparison of major and minor allelic copy number. In this case, the input for x.field and y.field is expected to be a length two character vector, where the first element is the major allele CN field and the second element is the minor allele CN field. this will be run if allelic = TRUE.
#'
#' another variation is comparison of copy-neutral LOH loci. In this case, the input for x.field and y.field is also a length two character vector (same as described above). this will be run if cnloh = TRUE. However, the output will be precision, recall, and F1 score of genomic tiles with cnLOH. Note that this comparison is non-symmetric and y is assumed to be the ground truth.
#' 
#' "unmappable" genomic tiles may also be masked by supplying a GRanges file path for the argument mask. coverage masks for hg19 and hg38 are included in extdata.
#'
#' @examples
#' ## compare the copy number information stored in a GRanges (x) with a data.table (y)
#' res = gBenchmark::benchmark_cn(x = system.file("extdata", "gr.rds", package = "gBenchmark"),
#'                                y = system.file("extdata", "dt.txt", package = "gBenchmark"),
#'                                x.field = "cn",
#'                                y.field = "cn",
#'                                tile.width = 1e3)
#' 
#' ## compare the copy number information stored in a gGraph (x) to a GRanges (y)
#' res = gBenchmark::benchmark_cn(x = system.file("extdata", "gg.rds", package = "gBenchmark"),
#'                                y = system.file("extdata", "gr.txt", package = "gBenchmark"),
#'                                x.field = "cn",
#'                                y.field = "cn",
#'                                tile.width = 1e3)
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
#' @param mask (character) path to coverage mask
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
                        mask = system.file("extdata", "hg19.mask.rds", package = "gBenchmark"),
                        verbose = FALSE)
{
    chrom.sizes = grab_chrom_sizes(genome = genome, std.only = std.only)
    

    ## tile eligible regions
    tiles.gr = gUtils::gr.tile(gUtils::si2gr(si = chrom.sizes), width = tile.width)

    ## grab input data
    gr1 = read_segs(x, field = x.field, allelic = allelic, cnloh = cnloh)
    gr2 = read_segs(y, field = y.field, allelic = allelic, cnloh = cnloh)

    ## remove tiles with > 90% overlap with the mask
    if (check_file(mask))
    {
        mask.gr = readRDS(mask)
        tiles.gr = tiles.gr %Q% (tiles.gr %O% mask.gr)
    }

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

#' @name benchmark_bnd
#' @title benchmark_bnd
#'
#' @description
#' compare two sets of breakend calls
#'
#' @details
#' given two sets of breakend calls x and y, with y being the 'ground truth' or 'gold standard'
#' calculate precision, recall, and f1 score of x
#'
#' x and y can be (signed) GRangesList, GRanges, gGraph (from which junctions/loose ends are extracted), or Junctions (or a file path containing any of these objects). in each of these cases, x.field and y.field do not have to be supplied.
#'
#' in addition, x and y can be a segmentation output (such as GRanges) with some field containing copy number (supplied as x.field or y.field). using the copy number change points, a set of signed breakends are inferred.
#'
#' if ignore.strand == TRUE, the orientation of breakends is ignored and two breakends in potentially opposite directions may be considered overlapping.
#'
#' @examples
#' ## compare a gGraph (x) with a set of junction calls provided as vcf (y)
#' res = gBenchmark:::benchmark_bnd(x = system.file("extdata", "gg.rds", package = "gBenchmark"),
#'                                 y = system.file("extdata", "svaba.subset.vcf", package = "gBenchmark"))
#' 
#' ## compare a gGraph (x) with a set of junction calls (provided as grl) (y)
#' res = gBenchmark:::benchmark_bnd(x = system.file("extdata", "gg.rds", package = "gBenchmark"),
#'                                 y = system.file("extdata", "grl.rds", package = "gBenchmark"))
#' 
#' ## compare a GRanges (x) with a set of junction calls (provided as grl) (y) and ignore the strand
#' res = gBenchmark:::benchmark_bnd(x = system.file("extdata", "gr.rds", package = "gBenchmark"),
#'                                 y = system.file("extdata", "grl.rds", package = "gBenchmark"),
#'                                 ignore.strand = TRUE)
#' 
#'
#' @param x character vector containing a file path, GRangesList, gGraph, GRanges, or Junction object (sv calls)
#' @param y character vector containing a file path, GRangesList, gGraph, GRanges, or Junction object (ground truth)
#' @param x.field (character) default NA
#' @param y.field (character) default NA
#' @param ignore.strand (logical) perform unsigned comparison of breakends? default FALSE
#' @param max.dist (numeric) maximum distance between breakends to be considered a match, default 1 Kbp
#' @param genome (character) one of hg19, hg3. default hg19
#' @param std.only (logical) standard chromosomes only? default TRUE
#' @param mask (character) path to coverage mask
#' @param verbose (logical) message stuff? default FALSE
benchmark_bnd = function(x, y,
                         x.field = NA_character_,
                         y.field = NA_character_,
                         ignore.strand = FALSE,
                         max.dist = 1e3,
                         genome = "hg19",
                         std.only = TRUE,
                         mask = system.file("extdata", "hg19.mask.rds", package = "gBenchmark"),
                         verbose = FALSE)
{
    chrom.sizes = grab_chrom_sizes(genome = genome, std.only = std.only)
    bnd1 = read_bnds(x, field = x.field, ignore.strand = ignore.strand)
    bnd2 = read_bnds(y, field = y.field, ignore.strand = ignore.strand)

    ## keep only things overlapping with desired chromosomes
    bnd1 = bnd1 %Q% (bnd1 %^% si2gr(chrom.sizes))
    bnd2 = bnd2 %Q% (bnd2 %^% si2gr(chrom.sizes))

    ## exclude coverage mask
    if (check_file(mask))
    {
        mask.gr = readRDS(mask)
        bnd1 = bnd1 %Q% (!bnd1 %^% mask.gr)
        bnd2 = bnd2 %Q% (!bnd2 %^% mask.gr)
    }

    ov.dt = data.table(row = numeric(), col = numeric(), dist = numeric())
    if (length(bnd1) && length(bnd2))
    {
        ov = gUtils::gr.dist(gr1 = bnd1, gr2 = bnd2, ignore.strand = ignore.strand)
        ov.dt = data.table(which(ov < max.dist, arr.ind = TRUE))
        ov.dt[, dist := ov[cbind(row, col)]]
    }

    res = data.table(n.bnd.x = length(bnd1),
                     n.bnd.y = length(bnd2),
                     tp = length(unique(ov.dt[, row])))

    res[, fp := n.bnd.x - tp]
    res[, fn := n.bnd.y - length(unique(ov.dt[, col]))]
    res[, precision := ifelse(tp > 0, tp / (tp + fp), 0)]
    res[, recall := ifelse(tp > 0, tp / (tp + fn), 0)]
    res[, f1 := ifelse(precision > 0 & recall > 0, 2 * precision * recall / (precision + recall), 0)]

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
#' @param nofield (logical) ignore field? (e.g. we just want segmentation) default FALSE
#' @param xy.sub (logical) if '23' and '24' are provided as seqnames, replace them with X and Y? default TRUE
#'
#' @return GRanges with field "score"
#' score represents copy number if "cnloh" is TRUE
#' additional field "allele" with values "major" and "minor" if "allelic" is TRUE
read_segs = function(x,
                     field = "cn",
                     allelic = FALSE,
                     cnloh = FALSE,
                     nofield = FALSE,
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
        if (!x[, .N]) { stop("empty table given as input!") }
        if (nofield)
        {
            x[, score := NA]
            field = "score"
        }
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
            x[seqnames == "23", seqnames := "X"]
            x[seqnames == "24", seqnames := "Y"]
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
        if (!length(x)) { stop("empty GRanges given as input!") }
        if (nofield)
        {
            GenomicRanges::mcols(x)[, "score"] = NA
            field = "score"
        }
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

#' @name seg2bnd
#' @title seg2bnd
#'
#' @description
#' convert segments with 'score' field representing CN to breakends
#'
#' @param gr
#'
#' @return signed GRanges representing breakends
seg2bnd = function(gr)
{
    dt = as.data.table(gr)[order(start),][order(seqnames),]
    dt[, ":="(prev.seqnames = data.table::shift(seqnames, type = "lag"),
              next.seqnames = data.table::shift(seqnames, type = "lead"),
              prev.score = data.table::shift(score, type = "lag"),
              next.score = data.table::shift(score, type = "lead"))]
    cncp.dt = rbind(dt[(seqnames == prev.seqnames) & (score > prev.score),
                       .(seqnames, start, end = start, strand = "+")],
                    dt[(seqnames == next.seqnames) & (score > next.score),
                       .(seqnames, start = end, end, strand = "-")])

    cncp.gr = GRanges()
    if (cncp.dt[, .N])
    {
        cncp.gr = GRanges(seqnames = cncp.dt[, seqnames],
                          ranges = IRanges::IRanges(start = cncp.dt[, start],
                                                    width = 1),
                          strand = cncp.dt[, strand])
    }
    return(cncp.gr)
}

#' @name read_bnds
#' @title read_bnds
#'
#' @description
#' grab breakends from input file or object
#'
#' @param x character vector containing file path, or GRanges, gGraph, Junction
#' @param field (character) copy number field, default NA
#' @param ignore.strand (logical) return output with no strand information (default FALSE)
#'
#' @return GRanges (possibly unsigned) with breakend locations
read_bnds = function(x, field = NA_character_, ignore.strand = FALSE)
{
    ## if field is supplied, assume that we have to read copy number
    if (!is.na(field))
    {
        x = read_segs(x, field = field)
        x = seg2bnd(x)
    }
    else if (inherits(x, "character"))
    {
        if (!check_file(x)) { stop("invalid file supplied: ", x) }
        if (grepl(".rds$", x))
        {
            x = readRDS(x)
        }
        else if (grepl("(.vcf$)|(.vcf.gz$)|(.bedpe$)|(.bedpe.gz$)", x))
        {
            x = gGnome::jJ(rafile = x)
        }
        else if (grepl("(.txt$)|(.tsv$)|(.csv$)", x))
        {
            x = data.table::fread(file = x)
        }
        else
        {
            stop("invalid file type supplied: ", x)
        }
    }

    ## check if x is the correct R object type...
    if (inherits(x, "gGraph"))
    {
        out = GRanges(seqlengths = seqlengths(x$gr))
        if (length(x$junctions) && length(x$junctions[type == "ALT"]))
        {
            out = stack(x$junctions[type == "ALT"]$grl)[, c()]
        }
        if (length(x$loose) && length(x$loose %Q% (!terminal)))
        {
            loose.gr = x$loose %Q% (!terminal)
            out = c(out, loose.gr[, c()])
        }
    }
    else if (inherits(x, "GRangesList"))
    {
        out = stack(x)[, c()]
    }
    else if (inherits(x, "Junction"))
    {
        out = stack(x$grl)[, c()]
    }
    else if (inherits(x, "data.table"))
    {
        if ("seqnames" %in% names(x) &&
            "start" %in% names(x) &&
            "end" %in% names(x) &&
            "strand" %in% names(x))
        {
            x = copy(x)[start > 1,]
            out = GRanges(seqnames = x[, seqnames],
                          ranges = IRanges::IRanges(start = x[, start],
                                                    end = x[, end]),
                          strand = x[, strand])
        }
        else
        {
            stop("supplied data table contains invalid names")
        }
    }
    else if (inherits(x, "GRanges"))
    {
        out = gUtils::gr.start(x[, c()] %Q% (start(x) > 1))
    }
    else
    {
        stop("file contains invalid object!")
    }
    
    if (ignore.strand) { out = gUtils::gr.stripstrand(out) }
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

#' @name grab_chrom_sizes
#' @title grab_chrom_sizes
#'
#' @description
#' helper function to get seqlengths
#'
#' @param genome (character) hg19 or hg38?
#' @param std.only (logical) only "standard" contigs
#'
#' @return named numeric vector
grab_chrom_sizes = function(genome, std.only = TRUE)
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

    return(chrom.sizes)
}
