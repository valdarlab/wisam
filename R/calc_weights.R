#' Calculate wISAM Strain Weights
#'
#' Computes the per-strain regression weights that \code{\link{wisam}} would use
#' for a given weighting scheme, without running a genome scan. Useful when you
#' only need the weights themselves (e.g. to feed into another scan function such
#' as \code{miqtl::scan.h2lmm}) rather than wISAM's own heteroscedastic scan.
#'
#' @param y A N length vector of phenotypes for each individual organism, where N is the total number of individuals. Not required if \code{sample_vars} and \code{counts} are both supplied directly.
#' @param strains Either a factor/vector of length N giving each individual's strain, or a s by N incidence matrix that maps every individual to a strain (same format accepted by \code{\link{wisam}}). If \code{sample_vars}/\code{counts} are supplied directly, this may instead just be a length-s vector of strain names/IDs (or omitted).
#' @param weights A string specifying the weighting scheme to use. One of \code{"none"}, \code{"samplevars"}, \code{"zshrink"}, \code{"limma"}, \code{"limma_nr"}, \code{"vashr"}, \code{"vashr_single"}, \code{"counts"}, or \code{"user"}.
#' @param user_weights A s length vector of weights for each strain, used if \code{weights = "user"}.
#' @param use_individuals If \code{TRUE}, return one weight per individual (length N) instead of one per strain (length s). Ignored if \code{sample_vars}/\code{counts} are supplied directly.
#' @param sample_vars An optional s length vector of already-computed per-strain sample variances, e.g. when replicates have already been summarized upstream (as in this project's \code{protdata}, which carries \code{NOISE}/\code{NUM.OBS} columns rather than raw per-replicate values). If supplied together with \code{counts}, the per-strain grouping/aggregation of \code{y}/\code{strains} is skipped entirely and \code{y}/\code{strains} are not required (\code{weights = "zshrink"} is not supported in this mode, since it needs the raw per-individual data).
#' @param counts An optional s length vector of already-computed per-strain replicate counts, paired with \code{sample_vars} (see above).
#' @param verbose If \code{TRUE}, print which weighting scheme is being used.
#'
#' @return A list containing:
#' \itemize{
#' \item{\code{weights}: the multiplicative strain weights (what \code{wisam(..., export_weights = TRUE)} returns as \code{weights}).}
#' \item{\code{ds}: the regression weights actually passed to the scan (\code{1 / (weights * counts)} for shrinkage-based schemes; see \code{\link{wisam}} for details of each scheme).}
#' \item{\code{strains}: the strain identifiers, in the same order as \code{weights}/\code{ds}.}
#' \item{\code{sample_vars}: the observed per-strain sample variance of \code{y}.}
#' \item{\code{counts}: the number of individuals contributing to each strain.}
#' \item{\code{removed_strains}: strains dropped because they had zero variance or a single observation (only applies to \code{weights = "samplevars"}; \code{NULL} otherwise). Drop the corresponding rows from \code{G}/\code{K}/\code{X} yourself before calling \code{\link{wisam}} with \code{weights = "samplevars"} on the same data.}
#' }
#'
#' @seealso \code{\link{wisam}}, which calls this same weighting logic internally before running its scan.
#'
#' @import tidyverse
#' @import stringr
#' @importFrom limma squeezeVar
#'
#' @export
calc_weights <- function(y, strains, weights = "none", user_weights = NULL,
                          use_individuals = FALSE, sample_vars = NULL, counts = NULL,
                          verbose = FALSE){

  if (is.null(user_weights) & weights == "user") { stop('Must provide user_weights (vector of weights) if weights = "user".') }
  if (!(weights %in% c("none", "samplevars", "zshrink", "limma", "limma_nr", "vashr", "vashr_single", "counts", "user"))) {
    stop('Weights must be one of "none", "samplevars", "zshrink", "limma_nr", "limma", "vashr", "vashr_single", "counts", and "user"')
  }

  precomputed <- !is.null(sample_vars) && !is.null(counts)

  if (precomputed) {
    if (weights == "zshrink") { stop('weights = "zshrink" needs raw per-individual y/strains and cannot be used with precomputed sample_vars/counts.') }
    if (length(sample_vars) != length(counts)) { stop("sample_vars and counts must be the same length.") }
    noise <- sample_vars
    names <- if (!missing(strains) && !is.null(strains) && is.null(dim(strains))) strains else seq_along(noise)
  } else {
    if (missing(y)) { stop('Must provide y (vector of phenotypes), or sample_vars and counts directly.') }
    if (missing(strains)) { stop('Must provide strains, or sample_vars and counts directly.') }

    if (is.null(dim(strains))) {
      if (length(y) != length(strains)) { stop("Input dimensions don't match.") }
    } else {
      if (length(y) != ncol(strains)) { stop("Input dimensions don't match.") }
    }

    if (is.null(names(y))) {
      warning("y vector does not have names - please be sure that the y vector is in the same order as strains.")
    }
    if (!is.null(dim(strains))) {
      strains <- apply(strains, 2, function(x) x * c(1:nrow(strains))) %>% colSums()
    }

    pheno_long <- data.frame(y = y, strains = strains)
    pheno_means <- pheno_long %>% dplyr::group_by(strains) %>% dplyr::summarise(mean = mean(y),
                                                                                 noise = var(y),
                                                                                 counts = dplyr::n())
    pheno_means <- as.data.frame(pheno_means)
    pheno_long$noise <- pheno_means[match(pheno_long$strains, pheno_means$strains), "noise"]
    pheno_long$counts <- pheno_means[match(pheno_long$strains, pheno_means$strains), "counts"]

    if (use_individuals) {
      y <- pheno_long$y
      noise <- pheno_long$noise
      counts <- pheno_long$counts
      names <- pheno_long$strains
    } else {
      y <- pheno_means$mean
      noise <- pheno_means$noise
      counts <- pheno_means$counts
      names <- pheno_means$strains
    }
  }

  removed_strains <- NULL
  ## check for strains with 0 variance and take them out for sample variance estimates
  if (weights %in% c("samplevars")) {
    ind <- which(noise == 0 | counts == 1)
    if (length(ind) > 0) {
      removed_strains <- names[ind]
      noise <- noise[-ind]
      counts <- counts[-ind]
      names <- names[-ind]
      print(ind)
    }
  }

  ######### WEIGHTS
  sample_vars <- noise
  switch(EXPR = weights,
         samplevars = {
           if (verbose) print("caluclating samplevars")
           weights <- mean(sample_vars) / sample_vars
           ds <- 1 / (weights * counts)
         },
         zshrink = {
           if (verbose) print("calculating shrinkage estimates with integrated conditional likelihood ebayes")
           shrink_estimates <- zshrink(pheno_long, strains)$VarEstimated
           if (use_individuals) shrink_estimates <- rep(shrink_estimates, pheno_means$counts)
           weights <- (mean(shrink_estimates) / shrink_estimates)
           ds <- 1 / (weights * counts)
         }, limma_nr = {
           if (verbose) print("calculating shrinkage estimates with non-robust limma")
           shrink_estimates <- squeezeVar(sample_vars, counts - 1, robust = FALSE)$var.post
           weights <- mean(shrink_estimates) / shrink_estimates
           ds <- 1 / (weights * counts)
         },
         limma = {
           if (verbose) print("calculating shrinkage estimates with robust limma")
           shrink_estimates <- squeezeVar(sample_vars, counts - 1, robust = TRUE)$var.post
           weights <- (mean(shrink_estimates) / shrink_estimates)
           ds <- 1 / (weights * counts)
         },
         vashr = {
           if (verbose) print("calculating shrinkage estimates with vashr")
           shrink_estimates <- vashr::vash(sample_vars, df = counts[1] - 1)$sd.post
           weights <- mean(shrink_estimates) / shrink_estimates
           ds <- 1 / (weights * counts)
         },
         vashr_single = {
           #vashr cannot handle NA
           sample_vars[is.na(sample_vars)] <- 0
           if (verbose) print("calculating shrinkage estimates with vashr and one prior")
           shrink_estimates <- vashr::vash(sample_vars, singlecomp = TRUE, df = counts[1] - 1)$sd.post
           weights <- mean(shrink_estimates) / shrink_estimates
           ds <- 1 / (weights * counts)
         },
         none = {
           if (verbose) print("no weights used")
           weights <- rep(1, length(y))
           ds <- weights
         },
         counts = {
           if (verbose) print("using counts as weights")
           weights <- counts
           ds <- 1 / (weights)
         },
         user = {
           weights <- user_weights
           ds <- 1 / (weights * counts)
         },
         stop("unkown weights input")
  )

  list(weights = weights,
       ds = ds,
       strains = names,
       sample_vars = sample_vars,
       counts = counts,
       removed_strains = removed_strains)
}
