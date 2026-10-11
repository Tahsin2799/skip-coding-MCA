# Supplementary tables and technical appendix for the JCE manuscript.
#
#   Rscript JCE/supplement/supplement_tables.R
#
# S1 describes the design grid. S2 lists the contrasts the decision rule left
# undecided, with where each falls short; its strata come from the cellwise
# check recorded in methodology/v2/PRODUCTION_RESULTS_2026-09-25.md section 3,
# and its counts are asserted against the claims file here. S3, S4 and S9 are
# computed from the displays file. S10 and the missing fuzzy average note to S9
# come from the post hoc production-linked theory summary.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

ds <- readRDS("outputs/methodology/production_v03_displays.rds")
review <- readRDS("JCE/analysis/review_analysis.rds")
claims <- read.csv("outputs/methodology/production_v03_claims.csv")
theory <- readRDS("JCE/theory/theory_link_summary_r1000.rds")

METHODS <- c(
  oracle = "Oracle (reference)",
  simple_score = "Deprivation count",
  impute_adaptive = "Imputation, target models",
  impute_frozen = "Imputation, development models",
  drop_items = "Items dropped",
  raw_items = "Supervised model (reference)",
  zero = "Skip coded as \"no\"",
  passive_mca = "Subset MCA",
  complete_case = "Complete-case analysis",
  category = "Missing single (skip as its own category)"
)
SHORT <- c(oracle = "Oracle", simple_score = "Count", impute_adaptive = "Imp-T",
           impute_frozen = "Imp-D", drop_items = "Dropped", raw_items = "Superv",
           zero = "No", passive_mca = "Subset", complete_case = "CC",
           category = "Single")
METRICS <- c(cluster_rho_all = "Cluster rank correlation",
             cl10_misalloc = "Inclusion error",
             b10_hn_capture = "High-need share selected")
VARIANTS <- c(
  q1_low = "Diagnosis given disease .20",
  q1_high = "Diagnosis given disease .80",
  q0_high = "False-positive diagnosis .15",
  rho_zero = "Access unrelated to deprivation",
  rho_strong = "Access correlation -.5",
  kappa_zero = "No between-cluster spread in access",
  items_two = "Two advice items",
  items_eight = "Eight advice items",
  loadings_flat = "Equal loadings",
  ndev_small = "Development sample of 500",
  impute_m20 = "Twenty imputations",
  recog_flat = "Diagnosis unrelated to disease",
  proxy_block = "Always-observed proxies instead of advice",
  advice_strong = "Advice strongly informative about need",
  bmi_monotone = "Monotone BMI truth"
)

md_table <- function(d, digits = 3) {
  fmt <- function(x) {
    if (!is.numeric(x)) return(x)
    out <- formatC(round(x, digits) + 0, format = "f", digits = digits)
    # Proportions drop the leading zero; percentages (one decimal) keep it.
    if (digits >= 2) out <- sub("^(-?)0\\.", "\\1.", out)
    sub("^-((0?)\\.0+)$", "\\1", out)
  }
  d <- as.data.frame(lapply(d, fmt), check.names = FALSE)
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, \(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

# S2 counts are checked against the claims file, so the typed strata cannot
# drift from the verdicts they describe.
S2 <- tribble(
  ~metric, ~lead, ~other, ~clear, ~where, ~size,
  "cluster_rho_all", "zero", "passive_mca", 70, "Bangladesh 18 to 95, contact filter, 80% contact, age slope .5", ".017 to .190",
  "cl10_misalloc", "oracle", "drop_items", 60, "Bangladesh 18 to 95, r² .8, both filters", ".001 to .026",
  "cl10_misalloc", "simple_score", "drop_items", 60, "Bangladesh 18 to 95, r² .8, both filters", ".001 to .017",
  "cl10_misalloc", "drop_items", "raw_items", 65, "Pooled stratum Q3, r² .8, contact filter and one diagnosis-filter scenario", ".005 to .079",
  "cl10_misalloc", "impute_adaptive", "impute_frozen", 63, "Bangladesh 18 to 95, contact filter, 50% and 80% contact", ".0002 to .025",
  "cl10_misalloc", "raw_items", "passive_mca", 65, "Bangladesh 18 to 95, contact filter", ".008 to .267",
  "cl10_misalloc", "zero", "passive_mca", 68, "Bangladesh 18 to 95, contact filter, 80% contact", ".007 to .104",
  "b10_hn_capture", "impute_adaptive", "category", 70, "Pooled stratum Q3, diagnosis filter, 20% contact", ".006 to .173",
  "b10_hn_capture", "drop_items", "impute_adaptive", 70, "Bangladesh 18 to 95, contact filter, r² .4", ".002 to .074",
  "b10_hn_capture", "impute_adaptive", "passive_mca", 69, "Pooled stratum Q3, diagnosis filter, 20% contact", ".005 to .168",
  "cl10_misalloc", "oracle", "simple_score", 60, "Bangladesh 18 to 95, r² .8, both filters; identical in every replicate", ".0002 to .008"
)
undecided <- claims |> filter(claim %in% c("unresolved", "equivalent"))
stopifnot(nrow(undecided) == nrow(S2))
chk <- S2 |>
  mutate(key = paste(metric, pmin(lead, other), pmax(lead, other))) |>
  left_join(undecided |> mutate(key = paste(metric, pmin(better, worse), pmax(better, worse)),
                                cleared = pmax(pos, neg)) |> select(key, cleared, claim),
            by = "key")
stopifnot(!anyNA(chk$cleared), all(chk$cleared == chk$clear))

s2 <- chk |>
  transmute(Metric = METRICS[metric],
            `Leading method` = METHODS[lead],
            `Other method` = METHODS[other],
            Verdict = if_else(claim == "equivalent", "equivalent", "undecided"),
            `Scenarios clearing 3 MCSE (of 72)` = as.character(clear),
            `Where the others sit` = where,
            `Difference where it clears` = size)

s3_base <- ds$variants |>
  filter(block %in% c("stress", "controls")) |>
  mutate(variant = factor(VARIANTS[variant], levels = VARIANTS),
         variant = if_else(gate == "diagnosis",
                           factor(paste0(variant, ", diagnosis filter")), factor(variant)))
s3a <- s3_base |>
  filter(setting == "c20") |>
  mutate(method = factor(SHORT[as.character(arm)], levels = SHORT)) |>
  select(Scenario = variant, method, cluster_rho_all) |>
  pivot_wider(names_from = method, values_from = cluster_rho_all) |>
  select(Scenario, all_of(unname(SHORT))) |>
  arrange(match(sub(", diagnosis filter", "", Scenario), VARIANTS), Scenario)
s3b <- s3_base |>
  filter(setting %in% c("c20", "c80")) |>
  mutate(method = factor(SHORT[as.character(arm)], levels = SHORT)) |>
  select(Scenario = variant, method, setting, cluster_rho_all) |>
  pivot_wider(names_from = setting, values_from = cluster_rho_all) |>
  mutate(change = c80 - c20) |>
  select(Scenario, method, change) |>
  pivot_wider(names_from = method, values_from = change) |>
  select(Scenario, all_of(unname(SHORT))) |>
  arrange(match(sub(", diagnosis filter", "", Scenario), VARIANTS), Scenario)

s4 <- ds$coverage_gap |>
  filter(setting %in% c("c20", "c80")) |>
  mutate(col = paste(if_else(gate == "contact", "Contact", "Diagnosis"),
                     sub("c", "", setting))) |>
  select(arm, col, capture_all, flag_rate) |>
  pivot_wider(names_from = col, values_from = c(flag_rate, capture_all),
              names_glue = "{.value}, {col}%") |>
  mutate(arm = as.character(arm), Method = METHODS[arm]) |>
  arrange(match(arm, names(METHODS))) |>
  select(Method, starts_with("flag_rate"), starts_with("capture_all")) |>
  rename_with(\(x) sub("flag_rate", "Share flagged", sub("capture_all", "Most vulnerable fifth flagged", x)))

policy_table <- function(gate_name) {
  review$policies |>
    filter(gate == gate_name) |>
    mutate(Method = METHODS[as.character(arm)],
           Contact = paste0(sub("c", "", setting), "%")) |>
    arrange(match(as.character(arm), names(METHODS)), setting) |>
    transmute(Method, Contact,
              `Inclusion error` = 100 * cl10_misalloc,
              `True top decile selected` = 100 * cl10_precision,
              `High-need people reached` = 100 * cl10_hn_reached,
              `Residents in selected clusters` = 100 * residents_in_selected_clusters,
              `Individual slots used` = 100 * individual_budget_used)
}
big <- \(x) format(x, big.mark = ",", trim = TRUE)
DISEASE_LABELS <- c(auc_index_scored = "AUC of the index", auc_model_scored = "AUC of the model",
                    log_or_top = "Top-tier odds ratio, development",
                    log_or_top_target = "Top-tier odds ratio, target",
                    beta_within = "Within-cluster slope", beta_between = "Between-cluster slope")
TARGET_LABELS <- c(cluster_rho_all = "Cluster rank correlation", cl10_misalloc = "Inclusion error",
                   b10_hn_capture = "High-need share selected")
FLAG_LABELS <- c(separation_target = "Separation in the target", separation_dev = "Separation in development")
s6 <- review$exclusions |>
  transmute(`Disease measure` = DISEASE_LABELS[metric], `Flag applied` = FLAG_LABELS[validity_flag],
            `Rows excluded` = big(flagged_rows), `Finite values excluded` = big(finite_flagged_rows))
stopifnot(!anyNA(s6$`Disease measure`), !anyNA(s6$`Flag applied`))
POOL_LABELS <- c(Q1 = "Structural pool Q1", Q2 = "Structural pool Q2", Q3 = "Structural pool Q3",
                 Q4 = "Structural pool Q4", bd_allage = "Bangladesh, 18 to 95", donor_india = "India",
                 donor_nepal = "Nepal", donor_albania = "Albania", donor_bangladesh = "Bangladesh, 18 to 49")
COUNTRY <- c(india = "India", nepal = "Nepal", albania = "Albania", bangladesh = "Bangladesh")
s7 <- review$method_mean_sensitivity |>
  transmute(`Disease measure` = DISEASE_LABELS[a], `Targeting measure` = TARGET_LABELS[t],
            `Scenarios` = as.character(cells), Median = median_rho,
            Minimum = min_rho, Maximum = max_rho)
s8a <- review$donor_counts |>
  transmute(Pool = POOL_LABELS[population], People = big(people),
            Clusters = big(clusters),
            Ages = paste(min_age, max_age, sep = " to "),
            Women = big(women), Men = big(men))
s8b <- review$donor_flow |>
  transmute(Survey = COUNTRY[donor], `Age-eligible records, 18 to 49` = big(age_eligible),
            `Complete on six inputs` = big(six_inputs_complete))
s8c <- review$pool_composition |>
  filter(population %in% c("Q1", "Q2", "Q3", "Q4")) |>
  transmute(Pool = POOL_LABELS[population], Survey = COUNTRY[donor], People = big(people))

# S9: rank of each method on cluster rank correlation over the 45 combinations of
# population, filter and contact level at the reference scenario.
pools <- bind_rows(
  ds$cellwise |> filter(r2 == 0.6, q1_age_slope == 0) |>
    select(population, gate, setting, arm, cluster_rho_all),
  ds$variants |> filter(block %in% c("quadrant", "donor")) |>
    select(population, gate, setting, arm, cluster_rho_all)) |>
  mutate(arm = as.character(arm))
stopifnot(nrow(distinct(pools, population, gate, setting)) == 45)
rank_table <- function(d) {
  d |> group_by(population, gate, setting) |>
    mutate(rank = rank(-cluster_rho_all, ties.method = "min")) |> ungroup()
}
s9_ranks <- rank_table(pools)
s9 <- s9_ranks |> group_by(arm) |>
  summarise(Best = min(rank), Median = median(rank), Worst = max(rank),
            c20 = paste(unique(range(rank[gate == "contact" & setting == "c20"])), collapse = " to "),
            .groups = "drop") |>
  arrange(match(arm, names(METHODS))) |>
  transmute(Method = METHODS[arm], Best = as.character(Best), Median = as.character(Median),
            Worst = as.character(Worst), `Range, contact filter at 20%` = c20)
albania <- pools |> filter(population == "donor_albania", gate == "contact", setting == "c20")
stopifnot(nrow(theory$pools_fa) == 45)
fa_ranks <- bind_rows(pools, theory$pools_fa |> mutate(arm = "fuzzy_average")) |>
  rank_table() |>
  filter(arm == "fuzzy_average")
fa_gap <- theory$pools_fa |>
  left_join(pools |> filter(arm == "drop_items") |> select(population, gate, setting, drop = cluster_rho_all),
            by = c("population", "gate", "setting")) |>
  mutate(gap = cluster_rho_all - drop)

# S10: missing fuzzy average against every production method under the frozen
# rule and margins, 72 core scenarios, 1,000 replicates each (post hoc).
VERDICTS <- c(robust_superiority = "consistently ahead", robust_inferiority = "consistently behind",
              scenario_dependent = "scenario dependent", equivalent = "equivalent",
              unresolved = "undecided")
s10 <- theory$claims |>
  mutate(other = as.character(worse)) |>
  arrange(match(metric, names(METRICS)), match(other, names(METHODS))) |>
  transmute(Metric = METRICS[metric], `Compared with` = METHODS[other],
            `Missing fuzzy average is` = VERDICTS[claim],
            `Smallest difference` = min_diff, `Largest difference` = max_diff,
            `Largest MCSE` = worst_mcse)
stopifnot(!anyNA(s10$`Missing fuzzy average is`), nrow(s10) == 30)

# S11 to S14: the post hoc extensions at the Q3 reference configuration
# (JCE/extensions/extensions_summary.R on outputs/methodology/extensions_v1).
ext <- readRDS("JCE/extensions/extensions_summary.rds")
stopifnot(all(ext$replicates$replicates == 1000), ext$failed == 18)
EXT_METHODS <- c(METHODS[c("oracle", "simple_score", "drop_items", "impute_adaptive")],
                 impute_flic_adaptive = "Imputation, target FLIC models",
                 impute_frozen = METHODS[["impute_frozen"]],
                 impute_flic_frozen = "Imputation, development FLIC models",
                 rimca = "Regularized iterative MCA",
                 fuzzy_average = "Missing fuzzy average",
                 METHODS[c("raw_items", "zero", "passive_mca", "complete_case", "category")])
GATE_LABEL <- c(contact = "Contact filter", diagnosis = "Diagnosis filter")
ext_num <- \(x, k = 3) sub("^(-?)0\\.", "\\1.", formatC(round(x, k) + 0, format = "f", digits = k))

s11 <- ext$sweep |>
  filter(arm %in% names(EXT_METHODS)) |>
  select(gate, loading, arm, mean) |>
  pivot_wider(names_from = loading, values_from = mean, names_prefix = "Loading ") |>
  left_join(ext$crossover |> select(gate, arm, crossover, ahead_everywhere), by = c("gate", "arm")) |>
  left_join(ext$sweep |> group_by(gate, arm) |> summarise(gap = max(abs(diff)), .groups = "drop"),
            by = c("gate", "arm")) |>
  mutate(crossover = case_when(arm == "drop_items" ~ "",
                               gap < .005 ~ paste("within", ext_num(gap), "at every loading"),
                               ahead_everywhere ~ "ahead at every loading",
                               is.na(crossover) ~ "behind at every loading",
                               TRUE ~ ext_num(crossover, 2))) |>
  arrange(gate, match(arm, names(EXT_METHODS))) |>
  transmute(Filter = GATE_LABEL[gate], Method = EXT_METHODS[arm],
            across(starts_with("Loading")), `Overtakes dropping at` = crossover)
names(s11) <- sub("Loading 0\\.", "Loading .", names(s11))
stopifnot(nrow(s11) == 26)
sweep_mcse <- max(ext$sweep$mcse)

s12 <- ext$devsize |>
  left_join(ext$devsize_models, by = "n_dev") |>
  arrange(n_dev, match(arm, names(EXT_METHODS))) |>
  transmute(`Development sample` = format(n_dev, big.mark = ","),
            `Diagnosed in development` = ext_num(dev_eligible, 0),
            `Events per parameter` = ext_num(events_per_parameter, 1),
            Method = EXT_METHODS[arm], `Rank correlation` = mean, `Difference from dropping` = diff,
            `Median difference` = median_diff, `Share trailing by more than .10` = behind_10,
            MCSE = ext_num(mcse, 4))
s12_drop <- ext$devsize |> filter(arm == "fuzzy_average") |> transmute(n_dev, drop = mean - diff)

CODING_ARM <- c(category = "category", subset = "passive_mca", zero = "zero", fuzzy_average = "fuzzy_average")
s13 <- ext$devcontact_theory |>
  mutate(arm = CODING_ARM[coding]) |>
  left_join(ext$devcontact |> select(gate, dev_setting, setting, arm, rho) |>
              pivot_wider(names_from = setting, values_from = rho), by = c("gate", "dev_setting", "arm")) |>
  arrange(gate, match(arm, CODING_ARM), dev_setting) |>
  transmute(Filter = GATE_LABEL[gate], Coding = EXT_METHODS[arm],
            `Development contact` = paste0(sub("c", "", dev_setting), "%"),
            `Share asked, p` = p, `Inertia ratio` = ext_num(inertia_ratio, 2),
            `R², block` = if_else(is.na(r2_block), "", ext_num(r2_block)), `R², eligibility` = r2_eligibility,
            `Target 20%` = c20, `Target 50%` = c50, `Target 80%` = c80)
stopifnot(nrow(s13) == 24, !anyNA(s13$`Target 50%`))
worst_at_50 <- sum(ext$worst$worst_target == "c50")

PCA_LABEL <- c(no = "Skip coded 0", fuzzy = "Mean imputation", category = "Own indicators (yes, no, skipped)",
               drop = "Items dropped")
PCA_ARM <- c(no = "pca_no", fuzzy = "pca_fuzzy", category = "pca_category", drop = "pca_drop")
MCA_ARM <- c(no = "zero", fuzzy = "fuzzy_average", category = "category", drop = "drop_items")
pca_rho <- ext$pca_vs_mca |>
  pivot_longer(-c(gate, loading), names_to = "arm", values_to = "rho") |>
  group_by(gate, arm) |> summarise(lo = min(rho), hi = max(rho), .groups = "drop") |>
  mutate(range = if_else(abs(hi - lo) < .0005, ext_num(lo), paste(ext_num(lo), "to", ext_num(hi))))
s14 <-tidyr::expand_grid(gate = names(GATE_LABEL), coding = names(PCA_LABEL)) |>
  left_join(ext$pca_theory, by = c("gate", "coding")) |>
  mutate(pca_arm = PCA_ARM[coding], mca_arm = MCA_ARM[coding]) |>
  left_join(pca_rho |> select(gate, pca_arm = arm, pca = range), by = c("gate", "pca_arm")) |>
  left_join(pca_rho |> select(gate, mca_arm = arm, mca = range), by = c("gate", "mca_arm")) |>
  transmute(Filter = GATE_LABEL[gate], Coding = PCA_LABEL[coding],
            `Gap from closed form` = if_else(is.na(max_formula_diff), "",
                                                         formatC(max_formula_diff, format = "e", digits = 0)),
            `Squared correlation with eligibility` = if_else(is.na(eta2), "", ext_num(eta2)),
            `Inertia ratio` = if_else(is.na(inertia_ratio), "", ext_num(inertia_ratio, 2)),
            `R², block` = if_else(is.na(r2_block), "", ext_num(r2_block)),
            `R², eligibility` = if_else(is.na(r2_eligibility), "", ext_num(r2_eligibility)),
            `Rank correlation, principal components` = pca, `Rank correlation, MCA counterpart` = mca)
stopifnot(nrow(s14) == 8, !anyNA(s14$`Rank correlation, principal components`))

lines <- c(
  "# Supplementary material",
  "",
  "*Coding questionnaire skips in composite health indices: a plasmode simulation study.*",
  "A. F. M. Tahsin Shahriar and Nasrin Lipi. Code and the summary results behind every table and figure: https://doi.org/10.5281/zenodo.23303798 (all versions).",
  "",
  "## Table S1. Design grid, 52 configurations",
  "",
  "| Block | Configurations | What varies |",
  "|---|---:|---|",
  "| Core | 24 | two populations (pooled stratum Q3, with high deprivation and weakly associated inputs, ages 18 to 49; Bangladesh 2017-18, ages 18 to 95) x *r²* .4, .6, .8 x age slope of diagnosis 0, .5 x contact and diagnosis filters; each evaluated at 20, 50 and 80% contact, 72 scenarios |",
  "| One-at-a-time variations | 13 | around the reference scenario (pooled stratum Q3, *r²* .6, no age slope, contact filter): diagnosis given disease .20 and .80 (also under the diagnosis filter), false-positive diagnosis .15, access unrelated to deprivation, access correlation -.5, no between-cluster spread in access, two and eight advice items, equal loadings, a development sample of 500, twenty imputations |",
  "| Mechanism checks | 4 | diagnosis unrelated to disease; four always-observed proxies in place of the advice items; advice strongly informative about need; a monotone BMI truth |",
  "| Other structural strata | 3 | the three remaining cells of deprivation level against input association (coherence, A1), at the reference scenario |",
  "| Country pools | 8 | India, Nepal, Albania and Bangladesh at 18 to 49 x both filters, at the reference scenario |",
  "",
  "Every configuration: 1,000 replicates, nominal development size 1,000 at 50% contact (500 in one variation), nominal target size 20,000. Production (23 to 25 September 2026) completed all 1,000 replicates of every configuration with no failed replicate; realized contact and prevalence were within .002 of target, and 31,200 metrics recomputed in base R from the saved row-level artifacts of 104 replicates matched the stored values exactly. Whole-cluster sampling usually produces slightly fewer people; the actual target N sets every budget and denominator. After production, every replicate of all 52 configurations was regenerated from its seeds to compute the theory diagnostics and to score missing fuzzy average (A6.9); these results are post hoc.",
  "",
  "## Table S2. Contrasts the decision rule did not order",
  "",
  "The ten undecided contrasts and the one equivalence. No undecided contrast changed sign: in each, every scenario that clears three MCSE favours the same method.",
  "",
  md_table(s2),
  "",
  "## Table S3a. Cluster rank correlation in the variations and the mechanism checks, 20% contact",
  "",
  "Pooled stratum Q3 at the reference scenario, contact filter unless stated. Oracle; Count, deprivation count; Imp-T and Imp-D, imputation with target or development models; Dropped, items dropped; Superv, supervised model; No, skip coded as \"no\"; Subset, subset MCA; CC, complete-case analysis; Single, missing single (skip as its own category).",
  "",
  md_table(s3a),
  "",
  "## Table S3b. Change in cluster rank correlation from 20% to 80% contact",
  "",
  "Methods that read no filtered item (Oracle, Count, Dropped) cannot change, so their change is exactly zero.",
  "",
  md_table(s3b),
  "",
  "## Table S4. Fixed threshold carried from development",
  "",
  "Development quintile cut-points applied to the target. The share of the target in the top tier moves away from a fifth when contact in the target differs from the 50% of the development sample, and the share of the most vulnerable fifth flagged has to be read beside it. Means over the core configurations.",
  "",
  md_table(s4),
  "",
  "## Table S5a. Allocation outcomes under the contact filter",
  "",
  "Percentages, averaged across the 12 core configurations per filter. Cluster selection uses a 10% cluster budget. True top decile selected is the share of selected clusters in the true top decile; high-need people reached is the share of all high-need people living in selected clusters. Individual slots used refers to the separate 10% individual budget. Every method filled the cluster budget in every core replicate.",
  "", md_table(policy_table("contact"), digits = 1), "",
  "## Table S5b. Allocation outcomes under the diagnosis filter",
  "", "Definitions and percentage units as in Table S5a.",
  "", md_table(policy_table("diagnosis"), digits = 1), "",
  "## Table S6. Separation exclusions in the correlation analysis",
  "",
  "Counts among 720,000 method-by-replicate-by-setting rows in 72 core scenarios. Development values repeat across three contact settings. The 6,453 development exclusions represent 2,151 distinct complete-case development fits: 2,143 under the diagnosis filter and eight under the contact filter. The four target exclusions are method-setting rows. Other targeting metrics on these rows remain included.",
  "", md_table(s6), "",
  "## Table S7. Post-review sensitivity using method means",
  "",
  "Spearman correlation across the ten method means within each scenario, after metric-specific validity filtering. The table gives the median and range over 72 scenarios. Inclusion error is sign-reversed. This is a descriptive sensitivity added after production; it does not change the frozen paired-contrast rule.",
  "", md_table(s7), "",
  "## Table S8a. Complete-input analysis pools used in production",
  "",
  "Counts reconstructed from the manifest-matched donor cache and frozen quadrant membership. Pools overlap. The structural pools Q1 to Q4 are not exclusively Indian.",
  "", md_table(s8a), "",
  "## Table S8b. Age eligibility and completeness in the country pools",
  "",
  "The first count includes all harmonized PR records aged 18 to 49, before requiring the six inputs. Attrition includes the BMI and individual-questionnaire sampling universes, not only item nonresponse. For Bangladesh ages 18 and over, Table S8a gives the separate complete-input pool.",
  "", md_table(s8b), "",
  "## Table S8c. Country composition of the structural pools",
  "", md_table(s8c), "",
  "## Table S9. Rank of each method on cluster rank correlation across the nine populations",
  "",
  "Ranks from 1 (best) to 10 among the ten prespecified methods over the 45 combinations of population, filter and contact level at the reference scenario (*r²* .6, no age slope); the last column is the range over the nine populations at 20% contact under the contact filter.",
  "", md_table(s9), "",
  sprintf("In the Albanian pool at 20%% contact under the contact filter, dropping the items reached %s, below both imputers (%s and %s), while the deprivation count held at %s.",
          sub("^0", "", sprintf("%.3f", albania$cluster_rho_all[albania$arm == "drop_items"])),
          sub("^0", "", sprintf("%.3f", albania$cluster_rho_all[albania$arm == "impute_adaptive"])),
          sub("^0", "", sprintf("%.3f", albania$cluster_rho_all[albania$arm == "impute_frozen"])),
          sub("^0", "", sprintf("%.3f", albania$cluster_rho_all[albania$arm == "simple_score"]))),
  "",
  sprintf("Missing fuzzy average (post hoc), ranked among eleven methods in the same %d combinations: best %d, median %g, worst %d. Its cluster rank correlation minus that of dropping the items ranged from %s to %s.",
          nrow(fa_ranks), min(fa_ranks$rank), median(fa_ranks$rank), max(fa_ranks$rank),
          sub("^0", "", sprintf("%.3f", min(fa_gap$gap))), sub("^0", "", sprintf("%.3f", max(fa_gap$gap)))),
  "",
  "## Table S10. Missing fuzzy average against each prespecified method (post hoc)",
  "",
  "The frozen decision rule and margins applied to missing fuzzy average paired with each production method in the 72 core scenarios, 1,000 replicates each. Differences are oriented so that a positive value favours missing fuzzy average (inclusion error is sign-reversed); smallest and largest are over the 72 scenarios. These contrasts were not prespecified and do not change the 132 production verdicts.",
  "", md_table(s10), "",
  "## Table S11. Cluster rank correlation by advice loading (post hoc)",
  "",
  sprintf("Pooled stratum Q3 at the reference scenario (*r²* .6, no age slope), development sample of 1,000 at 50%% contact; means over 1,000 replicates and the three target contact levels. The advice loading is the coefficient of *V* in the advice model (A3), .25 in production, where the ten production methods reproduce production exactly. The last column interpolates linearly where the mean paired difference from dropping the items turns positive. FLIC item models were fitted under the diagnosis filter only. Monte Carlo standard errors of the paired differences are at most %s. Figure 4A plots the differences.",
          ext_num(sweep_mcse, 4)),
  "", md_table(s11), "",
  "## Table S12. Development sample size under the diagnosis filter (post hoc)",
  "",
  sprintf("Advice loading .25, development at 50%% contact; 1,000 replicates per size, 3,000 paired values per row over the three target contact levels. Differences are from dropping the items, which scored %s, %s and %s at 1,000, 2,000 and 4,000. Events per parameter is the smaller outcome class over the number of coefficients, averaged over the four item models; under the diagnosis filter the people asked are the diagnosed. Samples of 2,000 and 4,000 draw their own seed streams. Figure 4B plots the differences.",
          ext_num(s12_drop$drop[1]), ext_num(s12_drop$drop[2]), ext_num(s12_drop$drop[3])),
  "", md_table(s12), "",
  "## Table S13. Development contact against target contact (post hoc)",
  "",
  sprintf("Advice loading .25, development sample of 1,000 at 20, 50 or 80%% contact, target evaluated at 20, 50 and 80%%. The share asked, the inertia ratio and the R² of the first axis with the filtered block's own axis and with eligibility come from development data. Means over 1,000 replicates; Monte Carlo standard errors at most %s. Dropping the items scored %s in every cell. Missing single, subset MCA and the skip coded \"no\" ranked worst at 50%% target contact in %d of the 18 combinations of filter, coding and development contact. Figure 4C plots the three codings.",
          ext_num(max(ext$devcontact$mcse), 4),
          ext_num(mean(ext$devcontact$rho[ext$devcontact$arm == "drop_items"])), worst_at_50),
  "", md_table(s13), "",
  "## Table S14. Principal component indices on indicator variables (post hoc)",
  "",
  sprintf("Every input category, and each advice item as coded, enters as a 0/1 indicator; the score is the first principal component of their correlation matrix, as for asset indices. The squared correlation with eligibility is the mean over the item's indicator columns. A filtered item enters as its yes indicator with the skip coded 0, as the same indicator with the skip replaced by the share of those asked answering yes, as yes, no and skipped indicators, or not at all. The closed forms of A6.8 were checked in every development sample (%s under the contact filter, %s under the diagnosis filter); the other diagnostic columns are means over the reference development (1,000 people at 50%% contact, all five advice loadings). Rank correlations are ranges over the five loadings, beside the MCA coding that treats the skip the same way.",
          format(ext$pca_theory$samples[ext$pca_theory$gate == "contact"][1], big.mark = ","),
          format(ext$pca_theory$samples[ext$pca_theory$gate == "diagnosis"][1], big.mark = ",")),
  "", md_table(s14), "",
  "## Figure S1. Disease validation measures against targeting measures",
  "",
  "Spearman correlations pooling ten methods and 1,000 replicates within each scenario, after metric-specific separation exclusions (A5). Disease measures use each method's scoreable population. Inclusion error and shortfall are sign-reversed. Dots are medians across 72 core scenarios and bars their range, not confidence intervals. The high-need share selected and the shortfall at the 10 percent budget are algebraically linked and give identical sign-adjusted correlations. All 54 medians are negative, from -.015 to -.197, with individual scenario correlations from -.638 to .570; all 18 medians using method means and the three principal metrics are also negative (Table S7). File: `JCE/figures/figS1.pdf`.",
  "",
  readLines("JCE/supplement/supplement_methods.md"), ""
)
writeLines(unlist(lines), "JCE/supplement/supplement_draft.md")
