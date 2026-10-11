# Donor harmonization for the v2 plasmode simulation.
# Maps four DHS surveys onto one set of ordered ladders.
#
# Mapping is by value LABEL, not by numeric code. Codes are not stable across
# DHS phases or countries: hv226 code 8 is wood in DHS-7 and kerosene in DHS-8,
# and v717 code 5 is an agricultural employee in Nepal and Albania but services
# in India. Only hv109 is verified standardized across all four donors and is
# the one variable mapped by code.

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tibble)
})

DONOR_ROOT <- "/Users/tahsin/Desktop/BS_project/Datasets"

donor_registry <- function() {
  tribble(
    ~donor,        ~label,               ~phase,  ~pr,                      ~ir,                      ~mr,                      ~occupation,
    "india",       "India 2019-21",      "DHS-7", "IAPR7EDT/IAPR7EFL.DTA",  "IAIR7EDT/IAIR7EFL.DTA",  "IAMR7EDT/IAMR7EFL.DTA",  "v717",
    "nepal",       "Nepal 2022",         "DHS-8", "NPPR82DT/NPPR82FL.DTA",  "NPIR82DT/NPIR82FL.DTA",  "NPMR82DT/NPMR82FL.DTA",  "v717",
    "albania",     "Albania 2017-18",    "DHS-7", "ALPR71DT/ALPR71FL.DTA",  "ALIR71DT/ALIR71FL.DTA",  "ALMR71DT/ALMR71FL.DTA",  "v717",
    "bangladesh",  "Bangladesh 2017-18", "DHS-7", "BDPR7RDT/BDPR7RFL.DTA",  NA,                       NA,                       "sb308"
  )
}

# ---- Label machinery -----------------------------------------------------

value_label <- function(x) {
  lv <- attr(x, "labels")
  names(lv)[match(as.numeric(x), as.numeric(lv))]
}

normalize_label <- function(s) {
  s <- tolower(as.character(s))
  s <- gsub("[^a-z0-9]+", " ", s)
  trimws(gsub(" +", " ", s))
}

# Ordered rules, first match wins. Returns NA for anything unmatched so the
# caller can report it rather than silently binning it.
map_by_label <- function(labels, rules) {
  out <- rep(NA_character_, length(labels))
  norm <- normalize_label(labels)
  for (rung in names(rules)) {
    hit <- is.na(out) & !is.na(norm) & grepl(rules[[rung]], norm)
    out[hit] <- rung
  }
  out
}

# ---- Ladders -------------------------------------------------------------

WATER_RULES <- c(
  surface     = "river|dam|lake|pond|stream|canal|irrigation",
  unimproved  = "unprotected|^other$|dug well|open well",
  on_premises = "piped into dwelling|piped to yard|piped to plot",
  improved    = "piped|tube well|borehole|protected|rainwater|bottled|tanker|cart with|public tap|standpipe|spring|well|ro plant|reverse osmosis"
)
WATER_LEVELS <- c("on_premises", "improved", "unimproved", "surface")

SANITATION_RULES <- c(
  open       = "no facility|bush|field",
  unimproved = "without slab|open pit|hanging|bucket|somewhere else|don t know where|dry toilet|^other$",
  improved   = "flush|ventilated improved|with slab|composting|biogas attached|pit toilet|latrine"
)
SANITATION_LEVELS <- c("improved_private", "improved_shared", "unimproved", "open")

FUEL_RULES <- c(
  none      = "no food cooked|^other$",
  clean     = "electricity|solar|lpg|liquefied petroleum|natural gas|biogas|ethanol|alcohol",
  kero_coal = "kerosene|paraffin|coal|lignite|gasoline|diesel",
  charcoal  = "charcoal",
  biomass   = "wood|straw|shrub|grass|crop|dung|sawdust|pellet|garbage|plastic|biomass"
)
FUEL_LEVELS <- c("clean", "kero_coal", "charcoal", "biomass")

WORK_RULES <- c(
  not_working = "not working|unemployed|student|retired",
  agriculture = "agricultur|farmer|cultivator|fisher|forestry|plantation|livestock|poultry|cattle|crop grower",
  services    = "professional|technical|managerial|clerical|sales|services|business|trader|land owner",
  manual      = "manual|household and domestic|domestic|construction|rickshaw|brick|road building|boatman|carpenter|mason|driver|tailor|seamstress|factory|blue collar|handicraft|labour|laborer|non agricultural worker|worker"
)
WORK_LEVELS <- c("not_working", "services", "manual", "agriculture")

EDUCATION_LEVELS <- c("none", "incomplete_primary", "primary",
                      "incomplete_secondary", "secondary", "higher")
BMI_LEVELS <- c("underweight", "normal", "overweight", "obese")

bmi_class <- function(bmi100) {
  x <- ifelse(!is.na(bmi100) & bmi100 < 6000, bmi100 / 100, NA_real_)
  cut(x, c(-Inf, 18.5, 25, 30, Inf), labels = BMI_LEVELS, right = FALSE)
}

# ---- Occupation ----------------------------------------------------------

# Bangladesh sb308 is country-specific and has no v717 analogue in the PR file.
# Codes follow R/analysis_main.R:45-47, which placed land owner (11) with the
# non-manual occupations rather than with farming. The v1 labor-taxing class is
# split here into manual and agriculture.
SB308_MAP <- c(
  "0" = "not_working", "61" = "not_working", "62" = "not_working",
  "11" = "services", "41" = "services", "51" = "services", "52" = "services",
  "16" = "manual", "21" = "manual", "22" = "manual", "23" = "manual", "31" = "manual",
  "12" = "agriculture", "13" = "agriculture", "14" = "agriculture", "15" = "agriculture"
)

occupation_from_recodes <- function(row) {
  read_one <- function(path, pre) {
    v <- paste0(pre, c("001", "002", "003", "717"))
    d <- read_dta(file.path(DONOR_ROOT, path), col_select = all_of(v))
    tibble(hv001 = as.numeric(d[[v[1]]]), hv002 = as.numeric(d[[v[2]]]),
           hvidx = as.numeric(d[[v[3]]]), work_label = value_label(d[[v[4]]]))
  }
  out <- read_one(row$ir, "v")
  if (!is.na(row$mr)) out <- bind_rows(out, read_one(row$mr, "mv"))
  out
}

# ---- Harmonization -------------------------------------------------------

harmonize_donor <- function(donor, age_min = 18, age_max = 49) {
  row <- donor_registry() |> filter(donor == !!donor) |> as.list()

  base <- c("hv001", "hv002", "hvidx", "hv021", "hv022", "hv005",
            "hv024", "hv025", "hv104", "hv105", "hv109", "hv201", "hv205",
            "hv225", "hv226", "hv009", "hv216", "hv270")
  sch <- read_dta(file.path(DONOR_ROOT, row$pr), n_max = 0)
  extra <- intersect(c("ha40", "hb40", "sbbm", "sb308"), names(sch))
  pr <- read_dta(file.path(DONOR_ROOT, row$pr), col_select = all_of(c(base, extra)))

  d <- tibble(
    donor      = donor,
    hv001      = as.numeric(pr$hv001),
    hv002      = as.numeric(pr$hv002),
    hvidx      = as.numeric(pr$hvidx),
    cluster    = paste0(donor, "_", as.numeric(pr$hv021)),
    stratum    = paste0(donor, "_", as.numeric(pr$hv022)),
    region     = value_label(pr$hv024),
    urban      = as.numeric(pr$hv025),
    weight     = as.numeric(pr$hv005) / 1e6,
    sex        = factor(as.numeric(pr$hv104), levels = c(1, 2), labels = c("male", "female")),
    age        = as.numeric(pr$hv105),
    hh_size    = as.numeric(pr$hv009),
    rooms      = as.numeric(pr$hv216),
    wealth     = as.numeric(pr$hv270),
    water_lab  = value_label(pr$hv201),
    sanit_lab  = value_label(pr$hv205),
    fuel_lab   = value_label(pr$hv226),
    shared     = as.numeric(pr$hv225),
    educ_code  = as.numeric(pr$hv109)
  )

  d$bmi100 <- if ("sbbm" %in% extra) as.numeric(pr$sbbm) else
    coalesce(as.numeric(pr$ha40), if ("hb40" %in% extra) as.numeric(pr$hb40) else NA_real_)

  if (row$occupation == "sb308") {
    d$work <- unname(SB308_MAP[as.character(as.numeric(pr$sb308))])
    d$work_label <- value_label(pr$sb308)
  } else {
    occ <- occupation_from_recodes(row)
    d <- left_join(d, occ, by = c("hv001", "hv002", "hvidx"))
    d$work <- map_by_label(d$work_label, WORK_RULES)
  }

  d <- d |>
    filter(age >= age_min, age <= age_max) |>
    mutate(
      water = map_by_label(water_lab, WATER_RULES),
      sanitation_base = map_by_label(sanit_lab, SANITATION_RULES),
      sanitation = case_when(
        sanitation_base == "improved" & shared == 1 ~ "improved_shared",
        sanitation_base == "improved" ~ "improved_private",
        .default = sanitation_base),
      fuel_raw = map_by_label(fuel_lab, FUEL_RULES),
      fuel = na_if(fuel_raw, "none"),
      education = if_else(educ_code %in% 0:5, EDUCATION_LEVELS[educ_code + 1], NA_character_),
      bmi = bmi_class(bmi100),
      water = factor(water, levels = WATER_LEVELS),
      sanitation = factor(sanitation, levels = SANITATION_LEVELS),
      fuel = factor(fuel, levels = FUEL_LEVELS),
      education = factor(education, levels = EDUCATION_LEVELS),
      work = factor(work, levels = WORK_LEVELS))

  prefix_levels(d)
}

INDEX_INPUTS <- c("water", "sanitation", "fuel", "bmi", "education", "work")

# predict.MCA fails when two active variables share a level label, and water and
# sanitation both carry "unimproved". Prefix every level with its variable name.
# v1 hit the same failure in the Phase 10 leakage check on the advice items.
prefix_levels <- function(d, vars = INDEX_INPUTS) {
  for (v in vars) {
    lv <- levels(d[[v]])
    if (!all(startsWith(lv, paste0(v, "_")))) levels(d[[v]]) <- paste0(v, "_", lv)
  }
  d
}

unmapped_labels <- function(d) {
  bind_rows(
    d |> filter(is.na(water), !is.na(water_lab)) |> count(variable = "water", label = water_lab),
    d |> filter(is.na(sanitation), !is.na(sanit_lab)) |> count(variable = "sanitation", label = sanit_lab),
    d |> filter(is.na(fuel_raw), !is.na(fuel_lab)) |> count(variable = "fuel", label = fuel_lab),
    d |> filter(is.na(work), !is.na(work_label)) |> count(variable = "work", label = work_label)
  )
}

complete_donor <- function(d) {
  d |> filter(if_all(all_of(INDEX_INPUTS), ~ !is.na(.x)))
}
