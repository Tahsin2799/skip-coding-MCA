# Generator contracts. Run: Rscript methodology/v2/tests/test_generator.R
#
# Every assertion here corresponds to a defect that reached a completed run, or
# to a property the design depends on. They are cheap. Run them before any job.

suppressPackageStartupMessages({library(dplyr); library(digest)})
root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
for (f in c("harmonize_donors", "plasmode_generate", "index_arms",
            "quadrant_populations", "pilot_v2_config")) {
  source(file.path(root, sprintf("methodology/v2/%s.R", f)))
}

harmonized <- lapply(readRDS(file.path(root, "outputs/methodology/donors_harmonized_geo.rds")),
                     prefix_levels)
keep_cols <- function(d) droplevels(as.data.frame(
  d[, c(INDEX_INPUTS, "age", "sex", "cluster")]))
pools <- lapply(build_quadrants(harmonized,
  quadrant_cells(file.path(root, "methodology/v2/quadrant_cells.csv"))), keep_cols)
pools[["bd_allage"]] <- keep_cols(build_age_pool(harmonized, "bangladesh", c(18L, Inf)))
for (nm in names(DONOR_POPULATIONS)) {
  pools[[nm]] <- keep_cols(build_age_pool(harmonized, DONOR_POPULATIONS[[nm]], c(18L, 49L)))
}
stopifnot("every declared population can be built" =
            setequal(names(pools), ALL_POPULATIONS))

POP <- "Q1"
spec <- pilot_spec(population = POP)
scl <- signal_scale(pools[[POP]], spec)
set.seed(digest2int("test_generator_1"))
sim <- generate_replicate(pools[[POP]], spec, 1000L, 20000L, scl)
dev <- sim$development

# ---- Rung scores come from the canonical ladder --------------------------
# Until v03 the spacing came from nlevels() after droplevels(), so a rung absent
# from a pool rescaled every other rung: charcoal is empty in all four donors,
# which put kero_coal at .50 instead of 1/3. A construct shared across
# populations needs spacing that does not depend on which categories appeared.
one <- pools[[POP]][1, ]
fuel_rung <- vapply(c("fuel_clean", "fuel_kero_coal", "fuel_biomass"), function(lv) {
  one$fuel <- factor(lv, levels = levels(pools[[POP]]$fuel))
  unname(rung_matrix(one)[, "fuel"])
}, numeric(1))
stopifnot("the fuel ladder keeps its canonical spacing" =
            isTRUE(all.equal(unname(fuel_rung), c(0, 1 / 3, 1))))
stopifnot("charcoal is absent from the pool but still holds its rung" =
            !"fuel_charcoal" %in% levels(pools[[POP]]$fuel) &&
            length(RUNG_LADDERS$fuel) == 4L)
one$education <- factor("education_none", levels = levels(pools[[POP]]$education))
stopifnot("education runs the other way, so none is the worst rung" =
            rung_matrix(one)[, "education"] == 1)

# ---- V is standardized on the pool, not on the draw ----------------------
# The development and target halves must mean the same thing by V, or the tier
# cuts carried from one to the other are cuts on two different scales. v2
# standardized the signal inside each draw, so a 1,000-row development half and
# a 20,000-row target were on two scales that differed by the sampling error of
# about fifty clusters.
raw_signal <- as.numeric(rung_matrix(dev, spec$bmi_risk) %*% spec$loadings)
set.seed(4)
v_pool <- true_vulnerability(dev, spec, scl)
set.seed(4)
v_raw <- true_vulnerability(dev, spec, c(center = 0, scale = 1))
stopifnot("V uses the constants it is handed, not the draw's own" =
            isTRUE(all.equal(v_pool - v_raw,
                             sqrt(spec$r2) * ((raw_signal - scl[["center"]]) / scl[["scale"]] -
                                                raw_signal))))
stopifnot("the pool constants are a property of the pool and the loadings" =
            isTRUE(all.equal(scl, signal_scale(pools[[POP]], spec))) &&
            !isTRUE(all.equal(unname(scl),
                              unname(signal_scale(pools[["Q2"]], spec)))))
stopifnot("both halves land near unit scale" =
            abs(stats::sd(dev$V) - 1) < 0.2 && abs(stats::sd(sim$targets$c50$V) - 1) < 0.1)

# ---- No column may shadow a model term ----------------------------------
# `impose_outcome` once wrote the contact noise into a column called `z`, and
# `glm(y ~ z + zbar, data = dev)` silently resolved it from the data frame, so
# the Mundlak within slope regressed the outcome on N(0, 1) for a whole pilot.
reserved <- c("z", "w", "zbar", "score", "tier", "lp")
for (f in c(list(dev), sim$targets)) {
  stopifnot("no generated column shadows a model term" = !any(reserved %in% names(f)))
}
stopifnot("the contact noise is named eps_contact" = "eps_contact" %in% names(dev))

# ---- The outcome is imposed before access -------------------------------
ys <- lapply(sim$targets, \(t) t$y)
stopifnot("y is identical across access settings" =
            all(vapply(ys[-1], identical, logical(1), ys[[1]])))
stopifnot("V is identical across access settings" =
            all(vapply(sim$targets[-1], \(t) identical(t$V, sim$targets[[1]]$V), logical(1))))

# ---- Access settings are common random numbers ---------------------------
# Contact, recognition and advice are thresholds on uniforms drawn once per
# person, so raising access can only add contacts. v2 redrew each setting with
# rbinom, which made individual movement across the grid partly noise and made
# any migration result a mixture of the access change and the redraw.
c20 <- sim$targets$c20
c50 <- sim$targets$c50
c80 <- sim$targets$c80
stopifnot("contact is nested across the access grid" =
            all(c20$contacted <= c50$contacted) && all(c50$contacted <= c80$contacted))
stopifnot("recognition is nested too" =
            all(c20$recognized <= c50$recognized) && all(c50$recognized <= c80$recognized))
stopifnot("the tie-break is a property of the person, not of the setting" =
            identical(c20$tiebreak, c80$tiebreak))
stopifnot("the tie-break is independent of the truth" =
            abs(stats::cor(c50$tiebreak, c50$V)) < 0.05 &&
            abs(stats::cor(c50$tiebreak, c50$age)) < 0.05)

# ---- The two gates ------------------------------------------------------
# Under the contact gate the advice items are asked of everyone who reached
# care; under the diagnosis gate only of people told they have the condition,
# which is what BDHS 2017-18 Q317B-Q318 does. The population must not change
# with the gate: the contrast between them is paired.
dx <- pilot_spec(population = POP, gate = "diagnosis")
set.seed(digest2int("test_generator_1"))
sim_dx <- generate_replicate(pools[[POP]], dx, 1000L, 20000L, scl)
stopifnot("the gate does not change who is in the population" =
            identical(sim$targets$c50$V, sim_dx$targets$c50$V) &&
            identical(sim$targets$c50$y, sim_dx$targets$c50$y) &&
            identical(sim$targets$c50$contacted, sim_dx$targets$c50$contacted) &&
            identical(sim$targets$c50$recognized, sim_dx$targets$c50$recognized))
items <- advice_names(dev)
stopifnot("four advice items" = length(items) == 4L)
for (v in items) {
  stopifnot("under the contact gate, advice is NA exactly when not contacted" =
              identical(is.na(sim$targets$c50[[v]]), sim$targets$c50$contacted == 0L))
  stopifnot("under the diagnosis gate, advice is NA exactly when not recognized" =
              identical(is.na(sim_dx$targets$c50[[v]]), sim_dx$targets$c50$recognized == 0L))
}
stopifnot("eligibility is the gate" =
            identical(sim$targets$c50$eligible, sim$targets$c50$contacted) &&
            identical(sim_dx$targets$c50$eligible, sim_dx$targets$c50$recognized))
stopifnot("every observed respondent under the diagnosis gate is recognized" =
            all(sim_dx$targets$c50$recognized[!is.na(sim_dx$targets$c50$advice1)] == 1L))
stopifnot("factual receipt still needs contact under both gates" =
            all(sim$targets$c50$advice1[sim$targets$c50$contacted == 0L] %in% NA_integer_))
stopifnot("nobody is recognized without contact" =
            all(dev$recognized[dev$contacted == 0L] == 0L))

# ---- Realized rates, for every population --------------------------------
# The intercepts are solved offline and pasted in, so the only thing standing
# between a stale constant and a wasted run is this assertion. Checked on a seed
# the solver never used.
nominal <- c(c20 = 0.20, c50 = 0.50, c80 = 0.80)
for (nm in names(pools)) {
  sp <- pilot_spec(population = nm)
  stopifnot("every population has a solved outcome intercept" =
              nm %in% names(sp$outcome_intercept) && !is.na(sp$outcome_intercept[[nm]]))
  set.seed(digest2int(paste0("test_rates_", nm)))
  s2 <- generate_replicate(pools[[nm]], sp, 1000L, 20000L)
  realized <- vapply(sp$settings, \(k) mean(s2$targets[[k]]$contacted), numeric(1))
  stopifnot("realized contact is within .02 of nominal at every setting" =
              all(abs(realized[names(nominal)] - nominal) < 0.02))
  stopifnot("realized outcome prevalence is within .02 of the target" =
              abs(mean(s2$targets$c50$y) - OUTCOME_PREVALENCE) < 0.02)
}

# ---- Every variant still realizes the grid it is read against ------------
# A variant that quietly lands at .18 is a confounded comparison rather than a
# stress test. Five of them change a term inside p_contact or the construction
# of V and carry their own solved intercepts; the rest must hit the grid on the
# baseline set, and this is what says so.
for (nm in names(c(STRESS, CONTROLS))) {
  sp <- modifyList(pilot_spec(population = "Q3"), variant_overrides(nm))
  n_dev <- if (nm %in% names(VARIANT_NDEV)) VARIANT_NDEV[[nm]] else 1000L
  set.seed(digest2int(paste0("test_variant_", nm)))
  s3 <- generate_replicate(pools[["Q3"]], sp, n_dev, 20000L)
  realized <- vapply(sp$settings, \(k) mean(s3$targets[[k]]$contacted), numeric(1))
  stopifnot("a variant realizes the nominal contact grid" =
              all(abs(realized[names(nominal)] - nominal) < 0.02))
  stopifnot("and the target outcome prevalence" =
              abs(mean(s3$targets$c50$y) - OUTCOME_PREVALENCE) < 0.02)
  stopifnot("no variant leaves an unsolved intercept" =
              !anyNA(unlist(sp[c("outcome_intercept", "contact_intercepts")])))
}
stopifnot("the item-count variants change the item count" =
            length(modifyList(spec, STRESS$items_two)$advice_probabilities) == 2L &&
            length(modifyList(spec, STRESS$items_eight)$advice_probabilities) == 8L)

# ---- The mechanism controls do what they say ------------------------------
flat <- modifyList(pilot_spec(population = "Q3"), CONTROLS$recog_flat)
set.seed(digest2int("test_control_recog"))
s_flat <- generate_replicate(pools[["Q3"]], flat, 1000L, 20000L)
t <- s_flat$targets$c50
rec <- c(case = mean(t$recognized[t$y == 1L & t$contacted == 1L]),
         non = mean(t$recognized[t$y == 0L & t$contacted == 1L]))
stopifnot("recog_flat removes the disease-to-recognition link" =
            abs(rec[["case"]] - rec[["non"]]) < 0.03)

proxy <- modifyList(pilot_spec(population = "Q3"), CONTROLS$proxy_block)
set.seed(digest2int("test_control_proxy"))
s_proxy <- generate_replicate(pools[["Q3"]], proxy, 1000L, 20000L)
stopifnot("the proxy block is asked of everyone" =
            !anyNA(s_proxy$targets$c20$advice1) && all(s_proxy$targets$c20$eligible == 1L))
stopifnot("and does not move with access" =
            identical(s_proxy$targets$c20$advice1, s_proxy$targets$c80$advice1))
stopifnot("while still carrying information about need" =
            stats::cor(s_proxy$targets$c50$advice1, s_proxy$targets$c50$V) > 0.2)

strong <- modifyList(pilot_spec(population = "Q3"), CONTROLS$advice_strong)
set.seed(digest2int("test_control_strong"))
s_strong <- generate_replicate(pools[["Q3"]], strong, 1000L, 20000L)
obs <- \(d) d$advice1[d$eligible == 1L]
stopifnot("advice_strong makes the items more informative about V" =
            stats::cor(obs(s_strong$targets$c50), s_strong$targets$c50$V[s_strong$targets$c50$eligible == 1L]) >
            stats::cor(obs(sim$targets$c50), sim$targets$c50$V[sim$targets$c50$eligible == 1L]))

mono <- modifyList(pilot_spec(population = "Q3"), CONTROLS$bmi_monotone)
u <- pools[["Q3"]][1, ]
u$bmi <- factor("bmi_underweight", levels = levels(pools[["Q3"]]$bmi))
stopifnot("bmi_monotone puts underweight at the bottom of the risk ladder" =
            rung_matrix(u, mono$bmi_risk)[, "bmi"] == 0 &&
            rung_matrix(u, spec$bmi_risk)[, "bmi"] == 0.70)

# ---- Recognition varies with age in the sourced direction ----------------
# Foti and Hardy 2025: adults 18-39 with care access are diagnosed more slowly
# than those 40 and over. That is recognition given contact, so a positive
# q1_age_slope must raise recognition with age among contacted cases.
tilt <- pilot_spec(q1_age_slope = 1.5, population = POP)
a <- impose_access(dev, "c50", pilot_spec(q1_age_slope = 0, population = POP))
b <- impose_access(dev, "c50", tilt)
rec_by_age <- function(d) {
  ok <- d$contacted == 1L & d$y == 1L
  old <- d$age[ok] > stats::median(d$age[ok])
  c(young = mean(d$recognized[ok][!old]), old = mean(d$recognized[ok][old]))
}
rb <- rec_by_age(b)
ra <- rec_by_age(a)
stopifnot("with a positive slope, older cases are recognized more often" =
            rb[["old"]] > rb[["young"]])
stopifnot("and the flat spec has no gradient of that size" =
            (rb[["old"]] - rb[["young"]]) > abs(ra[["old"]] - ra[["young"]]))
stopifnot("a zero slope is deterministic given the drawn uniforms" =
            identical(a$recognized, impose_access(dev, "c50",
              pilot_spec(q1_age_slope = 0, population = POP))$recognized))

# ---- Seeds reproduce ----------------------------------------------------
set.seed(digest2int("test_generator_1"))
again <- generate_replicate(pools[[POP]], spec, 1000L, 20000L, scl)
stopifnot("the same seed reproduces the replicate" =
            identical(sim$development$y, again$development$y),
          "and the target too" =
            identical(sim$targets$c50$contacted, again$targets$c50$contacted))
set.seed(digest2int("test_generator_2"))
other <- generate_replicate(pools[[POP]], spec, 1000L, 20000L, scl)
stopifnot("a different seed gives a different replicate" =
            !identical(sim$development$y, other$development$y))

# ---- Cluster draw and row ids -------------------------------------------
stopifnot("development and target clusters never collide" =
            !any(unique(dev$cluster) %in% unique(sim$targets$c50$cluster)))
stopifnot("cluster labels are the dev/tgt tagged form" =
            all(grepl("^dev_", unique(dev$cluster))))
stopifnot("target rows carry stable ids" =
            identical(sim$targets$c20$row_id, sim$targets$c80$row_id) &&
            identical(sim$targets$c50$row_id, seq_len(nrow(sim$targets$c50))))

# ---- The age term is a property of age, not of the population -------------
stopifnot("age enters per decade, from the published odds ratio" =
            isTRUE(all.equal(spec$age_log_or_per_decade, log(1.74))),
          "and there is no population-scaled age slope left" =
            is.null(spec$age_slope))

# ---- Rung matrix completeness is enforced, not silently dropped ----------
bad <- dev
bad$water[1] <- NA
stopifnot("an incomplete rung matrix aborts" =
            inherits(try(true_vulnerability(bad, spec, scl), silent = TRUE), "try-error"))

cat("test_generator.R: all contracts hold\n")
