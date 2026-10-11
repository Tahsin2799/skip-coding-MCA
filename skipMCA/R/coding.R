# Eligibility, item levels and the indicator matrix under each skip coding.

CODINGS <- c("fuzzy_average", "subset", "category", "no", "drop")

# Items that share one filter. NA marks a skip: every eligible respondent
# answered every item and every skipped respondent answered none.
filter_status <- function(data, items, eligible = NULL) {
  if (!length(items)) stop("`items` must name at least one filtered item.", call. = FALSE)
  unknown <- setdiff(items, names(data))
  if (length(unknown)) stop("Unknown items: ", paste(unknown, collapse = ", "), call. = FALSE)
  na <- vapply(items, function(v) is.na(data[[v]]), logical(nrow(data)))
  na <- matrix(na, nrow = nrow(data))
  if (is.null(eligible)) {
    eligible <- !na[, 1]
  } else if (is.character(eligible) && length(eligible) == 1L) {
    eligible <- data[[eligible]]
  }
  eligible <- as.logical(eligible)
  if (length(eligible) != nrow(data) || anyNA(eligible)) {
    stop("`eligible` must give one TRUE or FALSE per row.", call. = FALSE)
  }
  if (any(na[eligible, , drop = FALSE])) {
    stop("Some eligible respondents have missing answers. skipMCA reads NA as a skip, ",
         "so item nonresponse among the eligible must be handled first.", call. = FALSE)
  }
  if (!all(na[!eligible, , drop = FALSE])) {
    stop("Some skipped respondents have answers; the items must share one filter.", call. = FALSE)
  }
  eligible
}

check_substantive <- function(data, substantive) {
  unknown <- setdiff(substantive, names(data))
  if (length(unknown)) stop("Unknown substantive items: ", paste(unknown, collapse = ", "), call. = FALSE)
  if (any(vapply(substantive, function(v) anyNA(data[[v]]), logical(1)))) {
    stop("Substantive items must have no missing values.", call. = FALSE)
  }
}

observed_levels <- function(x) {
  x <- x[!is.na(x)]
  if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x)))
}

# Answer distribution of the eligible on each filtered item, over the
# development answer levels.
asked_profile <- function(data, items, levels, eligible) {
  stats::setNames(lapply(items, function(v) {
    x <- as.character(data[[v]][eligible])
    tab <- table(factor(x, levels = levels[[v]]))
    as.numeric(tab) / sum(tab)
  }), items)
}

# Everything fixed in development that a coding needs.
make_spec <- function(data, items, substantive, eligible, no_level) {
  lv <- stats::setNames(lapply(c(substantive, items), function(v) observed_levels(data[[v]])),
                        c(substantive, items))
  if (is.null(no_level)) no_level <- vapply(items, function(v) lv[[v]][1], character(1))
  if (length(no_level) == 1L) no_level <- rep(no_level, length(items))
  no_level <- stats::setNames(as.character(no_level), items)
  bad <- !mapply(function(v, l) l %in% lv[[v]], items, no_level)
  if (any(bad)) stop("`no_level` is not an observed answer of: ", paste(items[bad], collapse = ", "),
                     call. = FALSE)
  list(items = items, substantive = substantive, levels = lv, no_level = no_level,
       profile = asked_profile(data, items, lv, eligible), p = mean(eligible))
}

crisp <- function(x, levels, prefix) {
  z <- matrix(0, length(x), length(levels), dimnames = list(NULL, paste(prefix, levels, sep = "=")))
  j <- match(as.character(x), levels)
  ok <- !is.na(j)
  z[cbind(which(ok), j[ok])] <- 1
  z
}

# Indicator columns of one filtered item under a coding, and which columns
# the analysis keeps. Rows with an answer development never saw are flagged.
item_block <- function(x, eligible, coding, levels, profile, no_level, prefix) {
  z <- crisp(ifelse(eligible, as.character(x), NA), levels, prefix)
  unseen <- eligible & !(as.character(x) %in% levels)
  skipped <- !eligible
  keep <- rep(TRUE, ncol(z))
  if (coding == "fuzzy_average") {
    z[skipped, ] <- matrix(profile, sum(skipped), length(levels), byrow = TRUE)
  } else if (coding == "no") {
    z[skipped, match(no_level, levels)] <- 1
  } else {
    z <- cbind(z, as.numeric(skipped))
    colnames(z)[ncol(z)] <- paste0(prefix, "=<skipped>")
    keep <- c(keep, coding == "category")
  }
  list(z = z, keep = keep, unseen = unseen)
}

# The full indicator matrix of an index under a coding: substantive items
# crisp, filtered items as the coding says, nothing for dropped items.
index_block <- function(data, spec, eligible, coding, profile = spec$profile,
                        include_substantive = TRUE, include_items = TRUE) {
  parts <- list()
  unseen <- rep(FALSE, nrow(data))
  if (include_substantive) {
    for (v in spec$substantive) {
      z <- crisp(data[[v]], spec$levels[[v]], v)
      unseen <- unseen | rowSums(z) == 0
      parts[[length(parts) + 1]] <- list(z = z, keep = rep(TRUE, ncol(z)))
    }
  }
  if (include_items && coding != "drop") {
    for (v in spec$items) {
      b <- item_block(data[[v]], eligible, coding, spec$levels[[v]], profile[[v]],
                      spec$no_level[[v]], v)
      unseen <- unseen | b$unseen
      parts[[length(parts) + 1]] <- b
    }
  }
  list(z = do.call(cbind, lapply(parts, `[[`, "z")),
       keep = unlist(lapply(parts, `[[`, "keep")),
       q = length(parts), unseen = unseen)
}
