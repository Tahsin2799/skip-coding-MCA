# Post hoc extensions of the production simulation, asked for by the 7 October
# editor review. Nothing in methodology/v2 or outputs/methodology/production_v03
# is written; the production generator, arms and metric code are sourced
# unchanged.
#
#   sweep       advice loading on V of 0, .25, .5, .75 and 1 under both filters
#               (review item B4: where imputation starts to beat dropping)
#   devsize     development samples of 2,000 and 4,000 under the diagnosis
#               filter (A6: are the imputers' losses small-sample failures?)
#   devcontact  development contact of 20 and 80 percent under both filters
#               (A1: separate a development-to-target shift from the target's
#               own contact level)
#
# Every configuration also scores imputation with Firth-penalized item models
# and FLIC intercept correction (A6), missing fuzzy average, regularized
# iterative MCA (in a separate process, see prepare_rimca), and four principal
# component indices on indicator variables (the extension of the theory to
# asset-type indices). Configurations with a development sample of 1,000 reuse
# the production seed streams of the Q3 reference cell, so the reference
# configuration reproduces production exactly and every other loading or
# development contact level sees the same people and the same uniforms.
#
#   source("JCE/extensions/extensions.R")   # after the production sources

# ---- Configurations -----------------------------------------------------

REFERENCE_STREAM <- "Q3_nd1000_nt20000_r60_q00_base"

extension_grid <- function(studies = c("sweep", "devsize", "devcontact")) {
  row <- function(study, gate, loading = 0.25, n_dev = 1000L, dev_setting = "c50") {
    tibble(study = study, gate = gate, loading = loading, n_dev = n_dev,
           dev_setting = dev_setting)
  }
  parts <- list()
  if ("sweep" %in% studies) {
    g <- tidyr::expand_grid(gate = GATES, loading = c(0, 0.25, 0.5, 0.75, 1))
    parts$sweep <- row("sweep", g$gate, g$loading)
  }
  if ("devsize" %in% studies) parts$devsize <- row("devsize", "diagnosis", n_dev = c(2000L, 4000L))
  if ("devcontact" %in% studies) {
    g <- tidyr::expand_grid(gate = GATES, dev_setting = c("c20", "c80"))
    parts$devcontact <- row("devcontact", g$gate, dev_setting = g$dev_setting)
  }
  g <- dplyr::bind_rows(parts)
  g$population <- "Q3"
  g$n_target <- 20000L
  g$config_id <- sprintf("ext_%s_l%03.0f_nd%d_dev%s_g%s", g$study, g$loading * 100, g$n_dev,
                         g$dev_setting, substr(g$gate, 1, 3))
  # Development size changes the draw, so those configurations get their own
  # stream; the rest share the production stream of the reference cell.
  g$stream <- ifelse(g$n_dev == 1000L, "production", "extension")
  g$stream_id <- ifelse(g$n_dev == 1000L, REFERENCE_STREAM,
                        sprintf("Q3_nd%d_nt20000_r60_q00_ext", g$n_dev))
  g$spec <- lapply(seq_len(nrow(g)), function(i) {
    modifyList(pilot_spec(population = "Q3", gate = g$gate[i], r2 = 0.6),
               list(advice_v_loading = g$loading[i]))
  })
  stopifnot(!anyDuplicated(g$config_id))
  g
}

# generate_replicate() with the development contact level as an argument. With
# dev_setting = "c50" it is the production generator, draw for draw.
generate_extension <- function(pool, spec, n_dev, n_target, dev_setting, scl) {
  dev <- draw_population(pool, n_dev, "dev")
  tgt <- draw_population(pool, n_target, "tgt")
  for (nm in c("dev", "tgt")) {
    d <- get(nm)
    d$V <- true_vulnerability(d, spec, scl)
    d$A <- cluster_access(d$cluster, d$V, spec$access_rho)
    d <- impose_outcome(d, spec)
    assign(nm, draw_uniforms(d, length(spec$advice_probabilities)))
  }
  list(development = impose_access(dev, dev_setting, spec),
       targets = lapply(spec$settings, \(s) impose_access(tgt, s, spec)) |>
         setNames(spec$settings))
}

# ---- Firth-penalized item models with intercept correction ---------------

# FLIC (Puhr et al. 2017) through logistf::flic(): Firth's penalized likelihood
# (Firth 1993; Heinze and Schemper 2002) for the slopes, then the intercept
# re-estimated by maximum likelihood so that the mean predicted probability
# equals the observed event rate. The estimates replace the maximum likelihood
# coefficients of the production glm object, which complete_items() then
# predicts from unchanged. Aliased columns stay out, as in the production fit.
flic_refit <- function(fit) {
  ok <- !is.na(stats::coef(fit))
  x <- stats::model.matrix(fit)[, ok, drop = FALSE]
  stopifnot(colnames(x)[1] == "(Intercept)")
  # An intercept-only model has no slopes to penalize; FLIC then equals the
  # maximum likelihood fit.
  if (ncol(x) == 1L) return(fit)
  df <- data.frame(y = fit$y, x[, -1, drop = FALSE])
  names(df) <- c("y", paste0("x", seq_len(ncol(x) - 1L)))
  # flic() rebuilds its logistf() call with substitute() and evaluates it in its
  # own frame, so the data are passed by value through do.call().
  fit$coefficients[ok] <- stats::coef(do.call(logistf::flic, list(formula = y ~ ., data = df)))
  fit
}

# The production item models with FLIC estimates in place of maximum
# likelihood. The glm object is kept as a skeleton so that complete_items()
# predicts from it unchanged.
item_models_flic <- function(d) {
  mod <- suppressWarnings(item_models(d))
  mod$fits <- lapply(mod$fits, flic_refit)
  mod
}

# The production imputation family with FLIC item models.
prepare_flic <- function(dev_raw, spec) {
  mod <- item_models_flic(dev_raw)
  comp <- complete_items(dev_raw, mod, spec$m_impute)
  mps <- lapply(comp$frames, mca_dev, vars = c(INDEX_INPUTS, advice_names(dev_raw)))
  std <- lapply(mps, function(mp) {
    sgn <- if (dev_tiers(mp$coord_dev, dev_raw)$flip) -1 else 1
    list(sgn = sgn, mu = mean(mp$coord_dev), sd = stats::sd(mp$coord_dev))
  })
  z <- function(coords) Reduce(`+`, Map(\(x, st) st$sgn * (x - st$mu) / st$sd, coords, std)) / length(coords)
  out <- list(dev = dev_raw, s_dev = z(lapply(mps, `[[`, "coord_dev")), n_fix = comp$n_fix,
              family = "impute_flic", project = function(tgt, arm) {
                tm <- if (arm == "impute_flic_adaptive") item_models_flic(tgt) else mod
                tc <- complete_items(tgt, tm, spec$m_impute)
                list(score = z(Map(mca_coords, mps, tc$frames)), n_fix = tc$n_fix)
              })
  out$dt <- dev_tiers(out$s_dev, out$dev)
  out$dm <- dev_metrics(out$dev, out$dt)
  out
}

# ---- Principal component indices on indicator variables ------------------

# Asset-index practice: every category of each input as a 0/1 indicator, the
# first principal component of their correlation matrix as the score (Filmer
# and Pritchett 2001). A filtered binary item enters as its "yes" indicator,
# with the skip coded 0 ("no"), replaced by the share of those asked answering
# yes (mean imputation, the PCA counterpart of missing fuzzy average), given
# indicators of its own ("category": yes, no and skipped), or left out ("drop").
PCA_CODINGS <- c(pca_no = "no", pca_fuzzy = "fuzzy", pca_category = "category", pca_drop = "drop")

pca_matrix <- function(d, lv, items, coding, prof) {
  zx <- do.call(cbind, lapply(INDEX_INPUTS, function(v) {
    m <- outer(as.character(d[[v]]), lv[[v]], `==`) * 1
    colnames(m) <- paste(v, lv[[v]], sep = "=")
    m
  }))
  if (coding == "drop") return(zx)
  za <- do.call(cbind, lapply(items, function(v) {
    x <- d[[v]]
    sk <- is.na(x)
    yes <- as.numeric(x %in% 1L)
    m <- switch(coding,
      no = cbind(yes),
      fuzzy = {
        yes[sk] <- prof[[v]]
        cbind(yes)
      },
      category = cbind(yes, no = as.numeric(x %in% 0L), skipped = as.numeric(sk)))
    colnames(m) <- paste(v, colnames(m), sep = "=")
    m
  }))
  cbind(zx, za)
}

pca_first <- function(z) {
  mu <- colMeans(z)
  sdv <- apply(z, 2, stats::sd)
  keep <- sdv > 0
  zs <- sweep(sweep(z[, keep, drop = FALSE], 2, mu[keep]), 2, sdv[keep], "/")
  e <- eigen(stats::cor(z[, keep, drop = FALSE]), symmetric = TRUE)
  list(keep = keep, mu = mu, sd = sdv, v = e$vectors[, 1], eig = e$values,
       score = drop(zs %*% e$vectors[, 1]))
}

pca_score <- function(pf, z) {
  zs <- sweep(sweep(z[, pf$keep, drop = FALSE], 2, pf$mu[pf$keep]), 2, pf$sd[pf$keep], "/")
  drop(zs %*% pf$v)
}

asked_mean <- function(d, items) {
  setNames(lapply(items, function(v) {
    x <- d[[v]][!is.na(d[[v]])]
    if (length(x)) mean(x == 1L) else 0.5
  }), items)
}

prepare_pca <- function(dev, arm) {
  coding <- PCA_CODINGS[[arm]]
  items <- advice_names(dev)
  lv <- lapply(setNames(INDEX_INPUTS, INDEX_INPUTS), \(v) levels(droplevels(dev[[v]])))
  pf <- pca_first(pca_matrix(dev, lv, items, coding, asked_mean(dev, items)))
  out <- list(dev = dev, s_dev = pf$score, n_fix = 0L, family = arm, fit = pf,
              project = function(tgt, arm) {
                ok <- Reduce(`&`, lapply(INDEX_INPUTS, \(v) as.character(tgt[[v]]) %in% lv[[v]]))
                s <- pca_score(pf, pca_matrix(tgt, lv, items, coding, asked_mean(tgt, items)))
                s[!ok] <- NA_real_
                list(score = s, n_fix = 0L)
              })
  out$dt <- dev_tiers(out$s_dev, out$dev)
  out$dm <- dev_metrics(out$dev, out$dt)
  out
}

# Theory for principal components, from development data: the squared
# correlation of eligibility with each filtered column, empirical and by the
# closed forms of supplement A6.8, the inertia ratio lambda_A / lambda_X of the
# two blocks' correlation matrices (for principal components an eigenvalue is
# already a sum over columns), and the R-squared of the first component with the
# filtered block's own component and with eligibility.
pca_diagnostics <- function(dev) {
  items <- advice_names(dev)
  g <- as.numeric(observed_rows(dev))
  p <- mean(g)
  lv <- lapply(setNames(INDEX_INPUTS, INDEX_INPUTS), \(v) levels(droplevels(dev[[v]])))
  prof <- asked_mean(dev, items)
  zx <- pca_matrix(dev, lv, items, "drop", prof)
  px <- pca_first(zx)
  bind_rows(lapply(c("no", "fuzzy", "category"), function(cd) {
    z <- pca_matrix(dev, lv, items, cd, prof)
    za <- z[, setdiff(colnames(z), colnames(zx)), drop = FALSE]
    r2 <- vapply(seq_len(ncol(za)), \(j) if (stats::sd(za[, j]) > 0) stats::cor(za[, j], g)^2 else 0, numeric(1))
    pi_k <- unlist(prof)[sub("=.*", "", colnames(za))]
    lvl <- sub(".*=", "", colnames(za))
    formula <- if (cd == "fuzzy") {
      ifelse(lvl == "skipped", 1, 0)
    } else {
      ifelse(lvl == "skipped", 1,
             ifelse(lvl == "yes", pi_k * (1 - p) / (1 - p * pi_k),
                    (1 - pi_k) * (1 - p) / (1 - p * (1 - pi_k))))
    }
    pa <- pca_first(za)
    pf <- pca_first(z)
    tibble(coding = cd, p_dev = p, eta2_mean = mean(r2),
           eta2_max_diff = max(abs(r2 - formula)),
           lambda_x = px$eig[1], lambda_a = pa$eig[1], inertia_ratio = pa$eig[1] / px$eig[1],
           r2_block = stats::cor(pf$score, pa$score)^2,
           r2_eligibility = stats::cor(pf$score, g)^2)
  }))
}

# ---- Regularized iterative MCA ---------------------------------------------

# Regularized iterative MCA (Josse et al. 2012) through missMDA::imputeMCA()
# with its default two dimensions: the skipped cells of the indicator matrix are
# imputed by a low-rank reconstruction, and the index is the first axis of the
# MCA of the completed fuzzy indicator matrix (FactoMineR::MCA with tab.disj,
# as missMDA documents). The target is imputed the same way on its own data and
# placed by the transition formula on the development column standard
# coordinates. missMDA requires FactoMineR 2.17 or later while production ran
# 2.13, so this arm runs in its own process with that library first on the path
# (run_extensions.R --arms rimca --lib).
rimca_frame <- function(d, lv) {
  items <- advice_names(d)
  out <- data.frame(row.names = seq_len(nrow(d)))
  for (v in INDEX_INPUTS) out[[v]] <- factor(as.character(d[[v]]), levels = lv[[v]])
  for (v in items) {
    out[[v]] <- factor(ifelse(is.na(d[[v]]), NA, paste0(v, "_", d[[v]])), levels = paste0(v, "_", 0:1))
  }
  out
}

prepare_rimca <- function(dev) {
  lv <- lapply(setNames(INDEX_INPUTS, INDEX_INPUTS), \(v) levels(droplevels(dev[[v]])))
  don <- rimca_frame(dev, lv)
  imp <- missMDA::imputeMCA(don, ncp = 2)
  res <- FactoMineR::MCA(don, tab.disj = imp$tab.disj, ncp = 1, graph = FALSE)
  gamma <- as.matrix(res$var$coord)[, 1] / sqrt(res$eig[1, 1])
  q <- ncol(don)
  cols <- colnames(imp$tab.disj)
  s_dev <- as.numeric((imp$tab.disj[, cols] / q) %*% gamma[cols])
  out <- list(dev = dev, s_dev = s_dev, n_fix = 0L, family = "rimca",
              check = max(abs(abs(s_dev) - abs(as.matrix(res$ind$coord)[, 1]))),
              project = function(tgt, arm) {
                ok <- Reduce(`&`, lapply(INDEX_INPUTS, \(v) as.character(tgt[[v]]) %in% lv[[v]]))
                dt <- rimca_frame(tgt[ok, , drop = FALSE], lv)
                td <- missMDA::imputeMCA(dt, ncp = 2)$tab.disj
                z <- matrix(0, nrow(td), length(cols), dimnames = list(NULL, cols))
                common <- intersect(cols, colnames(td))
                z[, common] <- td[, common]
                s <- rep(NA_real_, nrow(tgt))
                s[ok] <- as.numeric((z / q) %*% gamma[cols])
                list(score = s, n_fix = 0L)
              })
  out$dt <- dev_tiers(out$s_dev, out$dev)
  out$dm <- dev_metrics(out$dev, out$dt)
  out
}

# ---- One replicate ------------------------------------------------------

EXT_ARMS <- c(ALL_ARMS, "fuzzy_average", "impute_flic_frozen", "impute_flic_adaptive",
              names(PCA_CODINGS))

ext_family <- function(arm) {
  if (arm %in% c("impute_flic_frozen", "impute_flic_adaptive")) "impute_flic"
  else if (arm == "rimca") "rimca"
  else if (arm %in% c(ALL_ARMS)) prep_family(arm)
  else arm
}

ext_prepare <- function(dev, fam, spec) {
  switch(fam,
    fuzzy_average = fuzzy_average_prep(dev),
    impute_flic = prepare_flic(dev, spec),
    rimca = prepare_rimca(dev),
    pca_no = , pca_fuzzy = , pca_category = , pca_drop = prepare_pca(dev, fam),
    prepare_family(dev, fam, spec))
}

# Events per parameter of the development item models, for A6.
item_model_size <- function(dev) {
  mod <- suppressWarnings(item_models(dev))
  ev <- vapply(mod$fits, \(f) min(sum(f$y), sum(1 - f$y)), numeric(1))
  tibble(dev_eligible = nrow(mod$obs), item_parameters = mean(vapply(mod$fits, \(f) length(stats::coef(f)), numeric(1))),
         smaller_class = mean(ev),
         events_per_parameter = mean(ev / vapply(mod$fits, \(f) length(stats::coef(f)), numeric(1))))
}

# FLIC item models are scored only under the diagnosis filter, where the
# development models rest on about 84 respondents (review item A6). Under the
# contact filter the item models have about 500 development and up to 16,000
# target respondents, where penalization changes little and refitting the
# target models with FLIC costs most of a replicate's time.
ext_arms <- function(cfg) {
  if (cfg$gate == "contact") setdiff(EXT_ARMS, c("impute_flic_frozen", "impute_flic_adaptive")) else EXT_ARMS
}

run_ext_replicate <- function(cfg, rep, arms = ext_arms(cfg)) {
  spec <- cfg$spec[[1]]
  key <- paste(cfg$stream, cfg$stream_id, rep, sep = "|")
  set.seed(digest2int(paste0(key, "|gen")))
  sim <- generate_extension(pools[[cfg$population]], spec, cfg$n_dev, cfg$n_target,
                            cfg$dev_setting, scales[[cfg$population]])
  preps <- list()
  for (fam in unique(vapply(arms, ext_family, character(1)))) {
    set.seed(digest2int(paste0(key, "|prep|", fam)))
    preps[[fam]] <- tryCatch(ext_prepare(sim$development, fam, spec), error = \(e) e)
  }
  rows <- list()
  for (arm in arms) {
    prep <- preps[[ext_family(arm)]]
    for (s in spec$settings) {
      set.seed(digest2int(paste0(key, "|score|", arm, "|", s)))
      rows[[length(rows) + 1]] <- if (inherits(prep, "error")) {
        tibble(setting = s, arm = arm, failed = TRUE, stage = "development",
               message = conditionMessage(prep))
      } else tryCatch(
        score_prepared(prep, arm, sim$targets[[s]], spec) |>
          mutate(setting = s, arm = arm, failed = FALSE, stage = NA_character_, message = NA_character_),
        error = \(e) tibble(setting = s, arm = arm, failed = TRUE, stage = "target",
                            message = conditionMessage(e)))
    }
  }
  rates <- tibble(
    setting = spec$settings,
    contact_rate = vapply(spec$settings, \(s) mean(sim$targets[[s]]$contacted), numeric(1)),
    recognized_rate = vapply(spec$settings, \(s) mean(sim$targets[[s]]$recognized), numeric(1)),
    eligible_rate = vapply(spec$settings, \(s) mean(sim$targets[[s]]$eligible), numeric(1)))
  dev <- sim$development
  # Diagnostics belong to the main pass only; the regularized iterative MCA pass
  # runs under a different FactoMineR.
  main <- "drop_items" %in% arms
  diag_mca <- if (main) tryCatch(coding_diagnostics(dev), error = \(e) tibble()) else tibble()
  diag_pca <- if (main) tryCatch(pca_diagnostics(dev), error = \(e) tibble()) else tibble()
  size <- if (main) tryCatch(item_model_size(dev), error = \(e) tibble(dev_eligible = sum(observed_rows(dev))))
          else tibble()
  list(metrics = bind_rows(rows) |> left_join(rates, by = "setting") |>
         mutate(config_id = cfg$config_id, rep = rep),
       diag_mca = diag_mca |> mutate(config_id = cfg$config_id, rep = rep),
       diag_pca = diag_pca |> mutate(config_id = cfg$config_id, rep = rep),
       size = size |> mutate(config_id = cfg$config_id, rep = rep))
}
