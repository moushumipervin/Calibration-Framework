###############################################################################
# FINAL_REAL_TRANSPORT_SIMULATION_MATCHED_IMPLEMENTATION.R
#
# Put this file and
# FINAL_REAL_TRANSPORT_SIMULATION_MATCHED_FUNCTIONS.R
# in the same working directory.
#
# Run ONLY this file.
###############################################################################

rm(list = ls())

source(
  "FINAL_REAL_TRANSPORT_SIMULATION_MATCHED_FUNCTIONS (1).R"
)

library(dplyr)
library(nnet)
library(mgcv)
library(DAAG)
library(haven)
library(numDeriv)


###############################################################################
# 1. Settings
###############################################################################

k <- 4
K <- 4
seed <- 2024202

xvars <- c(
  "age",
  "education",
  "black",
  "hispanic",
  "married",
  "nodegree",
  "re74",
  "re75"
)


###############################################################################
# 2. Data loading -- EXACTLY the previous real-data construction
#
# DO NOT:
#   - recreate/convert exp_dat$data_id
#   - convert PSID/CPS data_id to character
#   - use bind_rows() here
#   - add a complete-case filter
#   - reorder the rows
#
# The old analysis used base rbind(), which automatically coerced data_id
# to a common type and preserved the historical row order used by the seeded
# cross-fitting routine.
###############################################################################

exp_dat <- haven::read_dta(
  "http://www.nber.org/~rdehejia/data/nsw_dw.dta"
)

data("psid1", package = "DAAG")
data("cps1", package = "DAAG")


# -----------------------------
# CPS: EXACT old preparation
# -----------------------------
cps1_renamed <- cps1

names(cps1_renamed)[
  names(cps1_renamed) == "trt"
] <- "treat"

names(cps1_renamed)[
  names(cps1_renamed) == "educ"
] <- "education"

names(cps1_renamed)[
  names(cps1_renamed) == "hisp"
] <- "hispanic"

names(cps1_renamed)[
  names(cps1_renamed) == "marr"
] <- "married"

names(cps1_renamed)[
  names(cps1_renamed) == "nodeg"
] <- "nodegree"

# EXACT old code: integer sequence; do not force character.
cps1_renamed$data_id <-
  seq_len(
    nrow(cps1_renamed)
  )

cps1_renamed <-
  cps1_renamed[
    ,
    c(
      "data_id",
      "treat",
      "age",
      "education",
      "black",
      "hispanic",
      "married",
      "nodegree",
      "re74",
      "re75",
      "re78"
    )
  ]


# -----------------------------
# PSID: EXACT old preparation
# -----------------------------
psid1_renamed <- psid1

names(psid1_renamed)[
  names(psid1_renamed) == "trt"
] <- "treat"

names(psid1_renamed)[
  names(psid1_renamed) == "educ"
] <- "education"

names(psid1_renamed)[
  names(psid1_renamed) == "hisp"
] <- "hispanic"

names(psid1_renamed)[
  names(psid1_renamed) == "marr"
] <- "married"

names(psid1_renamed)[
  names(psid1_renamed) == "nodeg"
] <- "nodegree"

# EXACT old code.
psid1_renamed$data_id <-
  seq_len(
    nrow(psid1_renamed)
  )

psid1_renamed <-
  psid1_renamed[
    ,
    c(
      "data_id",
      "treat",
      "age",
      "education",
      "black",
      "hispanic",
      "married",
      "nodegree",
      "re74",
      "re75",
      "re78"
    )
  ]


# -----------------------------
# NSW: EXACT old preparation
# -----------------------------
nsw_treat <-
  exp_dat |>
  dplyr::filter(
    treat == 1
  ) |>
  dplyr::mutate(
    G = 1L
  )

nsw_ctrl <-
  exp_dat |>
  dplyr::filter(
    treat == 0
  ) |>
  dplyr::mutate(
    G = 2L
  )


# EXACT old group labels.
psid1_renamed$G <- 3
cps1_renamed$G  <- 4


# CRITICAL:
# The original script used base rbind(), not dplyr::bind_rows().
# Keep this exactly because it also preserves the old row ordering.
dat_all <- rbind(
  nsw_treat,
  nsw_ctrl,
  psid1_renamed,
  cps1_renamed
)

dat_all2 <- dat_all


# EXACT old d4 construction.
d4 <-
  dat_all2 |>
  dplyr::filter(
    G %in% c(
      1,
      2,
      3,
      4
    )
  ) |>
  dplyr::mutate(
    G = factor(G)
  )


# The old cross-fitting functions create/recreate ID internally.
# We set it once here only after reproducing the exact old row order.
d4$ID <- seq_len(
  nrow(d4)
)


###############################################################################
# 2A. Sanity checks against the historical analysis
###############################################################################

cat(
  "\nHistorical-data dimensions:\n"
)

print(
  table(d4$G)
)

cat(
  "\nFirst 10 G values (row-order check):\n"
)

print(
  head(
    as.character(d4$G),
    10
  )
)

###############################################################################
# 3. Group counts and experimental benchmark
###############################################################################

print(
  d4 |>
    dplyr::count(
      G
    )
)


benchmark <-
  mean(
    d4$re78[
      d4$G == "1"
    ]
  ) -
  mean(
    d4$re78[
      d4$G == "2"
    ]
  )


cat(
  "\nNSW randomized benchmark = ",
  benchmark,
  "\n\n",
  sep = ""
)


###############################################################################
# 4. ONE four-category multinomial source-membership fit
#
# Point-estimation fit is EXACTLY the old real-data multinomial fit.
###############################################################################

fml_multi <-
  as.formula(
    paste(
      "G ~",
      paste(
        xvars,
        collapse = "+"
      )
    )
  )

fit_multi <-
  nnet::multinom(
    fml_multi,
    data = dat_all2,
    trace = FALSE
  )

P <-
  predict(
    fit_multi,
    type = "probs"
  )

P <- as.matrix(P)

p1 <- P[, "1"]
p2 <- P[, "2"]
p3 <- P[, "3"]
p4 <- P[, "4"]

pNSW <- p1 + p2

# EXACT old real-data propensity ratios.
w1 <- pNSW / p1
w2 <- pNSW / p2
w3 <- pNSW / p3
w4 <- pNSW / p4


# Build the object needed by the NEW multinomial sandwich,
# without refitting or changing the point-estimation propensity model.
X_multi <-
  model.matrix(
    fml_multi,
    data = dat_all2
  )

B_multi <-
  as.matrix(
    coef(fit_multi)
  )

B_multi <-
  B_multi[
    c(
      "2",
      "3",
      "4"
    ),
    ,
    drop = FALSE
  ]

ps_obj <- list(
  fit = fit_multi,
  formula = fml_multi,
  X = X_multi,
  B = B_multi,
  phi_hat =
    as.vector(
      t(B_multi)
    ),
  P = P,
  xvars = xvars
)


###############################################################################
# 4A. Verify that the original point-estimation weights are being used
###############################################################################

cat(
  "\nOriginal multinomial transport-weight summaries:\n"
)

print(
  summary(
    w2[
      d4$G == "2"
    ]
  )
)

print(
  summary(
    w3[
      d4$G == "3"
    ]
  )
)

print(
  summary(
    w4[
      d4$G == "4"
    ]
  )
)

###############################################################################
# 5. Unweighted
###############################################################################

uw_nsw <- unweighted_function(
  d4,
  cat1 = 1,
  cat0 = 2
)

uw_psid <- unweighted_function(
  d4,
  cat1 = 1,
  cat0 = 3
)

uw_cps <- unweighted_function(
  d4,
  cat1 = 1,
  cat0 = 4
)


###############################################################################
# 6. IPW
#
# Point estimator = previous code.
# SE = new joint multinomial sandwich.
###############################################################################

ipw_nsw <- IPW_multinom_joint_final(
  d4 = d4,
  donor_g = 2,
  ps_obj = ps_obj,
  benchmark = benchmark
)

ipw_psid <- IPW_multinom_joint_final(
  d4 = d4,
  donor_g = 3,
  ps_obj = ps_obj,
  benchmark = benchmark
)

ipw_cps <- IPW_multinom_joint_final(
  d4 = d4,
  donor_g = 4,
  ps_obj = ps_obj,
  benchmark = benchmark
)


###############################################################################
# 7. AIPW
#
# STEP 1:
# Run the EXACT ORIGINAL point estimator first.
#
# STEP 2:
# Separately compute the NEW multinomial-phi sandwich SE.
#
# This prevents the variance revision from redefining the AIPW point estimate.
###############################################################################

# -----------------------------
# 7A. ORIGINAL point estimates
# -----------------------------

res_nsw_lm_old <-
  run_aipw_cv(
    cat0 = 2,
    w_use = w2,
    t0_model = "lm"
  )

res_nsw_gam_old <-
  run_aipw_cv(
    cat0 = 2,
    w_use = w2,
    t0_model = "gam"
  )


res_psid_lm_old <-
  run_aipw_cv(
    cat0 = 3,
    w_use = w3,
    t0_model = "lm"
  )

res_psid_gam_old <-
  run_aipw_cv(
    cat0 = 3,
    w_use = w3,
    t0_model = "gam"
  )


res_cps_lm_old <-
  run_aipw_cv(
    cat0 = 4,
    w_use = w4,
    t0_model = "lm"
  )

res_cps_gam_old <-
  run_aipw_cv(
    cat0 = 4,
    w_use = w4,
    t0_model = "gam"
  )


cat(
  "\n============== ORIGINAL AIPW POINT ESTIMATES ==============\n"
)

original_AIPW_points <- data.frame(
  Donor = rep(
    c(
      "NSW",
      "PSID",
      "CPS"
    ),
    each = 2
  ),
  Model = rep(
    c(
      "LM",
      "GAM"
    ),
    3
  ),
  Estimate = c(
    res_nsw_lm_old$ATE,
    res_nsw_gam_old$ATE,
    res_psid_lm_old$ATE,
    res_psid_gam_old$ATE,
    res_cps_lm_old$ATE,
    res_cps_gam_old$ATE
  ),
  Old_SE = c(
    res_nsw_lm_old$se,
    res_nsw_gam_old$se,
    res_psid_lm_old$se,
    res_psid_gam_old$se,
    res_cps_lm_old$se,
    res_cps_gam_old$se
  )
)

print(
  original_AIPW_points
)


# -----------------------------
# 7B. NEW variance calculations
# -----------------------------

aipw_nsw_lm <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 2,
    ps_obj = ps_obj,
    model = "lm",
    K = K,
    seed = seed,
    benchmark = benchmark
  )

aipw_nsw_gam <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 2,
    ps_obj = ps_obj,
    model = "gam",
    K = K,
    seed = seed,
    benchmark = benchmark
  )


aipw_psid_lm <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 3,
    ps_obj = ps_obj,
    model = "lm",
    K = K,
    seed = seed,
    benchmark = benchmark
  )

aipw_psid_gam <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 3,
    ps_obj = ps_obj,
    model = "gam",
    K = K,
    seed = seed,
    benchmark = benchmark
  )


aipw_cps_lm <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 4,
    ps_obj = ps_obj,
    model = "lm",
    K = K,
    seed = seed,
    benchmark = benchmark
  )

aipw_cps_gam <-
  AIPW_multinom_joint_final(
    d4 = d4,
    donor_g = 4,
    ps_obj = ps_obj,
    model = "gam",
    K = K,
    seed = seed,
    benchmark = benchmark
  )


# FORCE the final reported ATE to be the EXACT original run_aipw_cv() result.
# Only the SE comes from the new multinomial sandwich.
aipw_nsw_lm$ATE   <- res_nsw_lm_old$ATE
aipw_nsw_gam$ATE  <- res_nsw_gam_old$ATE

aipw_psid_lm$ATE  <- res_psid_lm_old$ATE
aipw_psid_gam$ATE <- res_psid_gam_old$ATE

aipw_cps_lm$ATE   <- res_cps_lm_old$ATE
aipw_cps_gam$ATE  <- res_cps_gam_old$ATE


###############################################################################
# 7C. Historical-table comparison
#
# These are ONLY diagnostics from the table you supplied.
###############################################################################

historical_table_AIPW <- data.frame(
  Donor = rep(
    c(
      "NSW",
      "PSID",
      "CPS"
    ),
    each = 2
  ),
  Model = rep(
    c(
      "LM",
      "GAM"
    ),
    3
  ),
  HistoricalRounded = c(
    1779,
    1780,
    1283,
    1426,
    1291,
    1683
  )
)

historical_table_AIPW$CurrentRounded <- round(
  original_AIPW_points$Estimate
)

historical_table_AIPW$Difference <- 
  historical_table_AIPW$CurrentRounded -
  historical_table_AIPW$HistoricalRounded


cat(
  "\n============== HISTORICAL TABLE AIPW CHECK ==============\n"
)

print(
  historical_table_AIPW
)


###############################################################################
# PACKAGE BENCHMARK METHODS
# EBPS / CBPS / EBCW
#
# oCBPS is intentionally omitted from the real-data analysis because the
# CBPSOptimal implementation does not support the ATT=1 transport formulation
# required for the external PSID/CPS donor comparisons.
###############################################################################

# Experimental NSW benchmark (kept on the full numerical scale for estimation)
benchmark <-
  mean(
    d4$re78[
      as.character(d4$G) == "1"
    ]
  ) -
  mean(
    d4$re78[
      as.character(d4$G) == "2"
    ]
  )


###############################################################################
# 7D. Construct source objects for package benchmark methods
###############################################################################

source_nsw <-
  make_transport_source_data(
    d4 = d4,
    donor_g = 2,
    xvars = xvars
  )

source_psid <-
  make_transport_source_data(
    d4 = d4,
    donor_g = 3,
    xvars = xvars
  )

source_cps <-
  make_transport_source_data(
    d4 = d4,
    donor_g = 4,
    xvars = xvars
  )


###############################################################################
# 7E. EBPS / CBPS / EBCW
###############################################################################

# -----------------------------
# NSW
# -----------------------------
ebps_nsw <-
  transport_ebps(
    source_nsw,
    benchmark
  )

cbps_nsw <-
  transport_cbps(
    source_nsw,
    benchmark
  )

ebcw_nsw <-
  transport_ebcw(
    source_nsw,
    benchmark
  )


# -----------------------------
# PSID
# -----------------------------
ebps_psid <-
  transport_ebps(
    source_psid,
    benchmark
  )

cbps_psid <-
  transport_cbps(
    source_psid,
    benchmark
  )

ebcw_psid <-
  transport_ebcw(
    source_psid,
    benchmark
  )


# -----------------------------
# CPS
# -----------------------------
ebps_cps <-
  transport_ebps(
    source_cps,
    benchmark
  )

cbps_cps <-
  transport_cbps(
    source_cps,
    benchmark
  )

ebcw_cps <-
  transport_ebcw(
    source_cps,
    benchmark
  )


###############################################################################
# 8. ET / HD / CE
#
# All three use:
#   - exact old GAM predictions;
#   - simulation entropy-specific g direction;
#   - simulation inverse-dual weight map;
#   - Newton + Armijo backtracking;
#   - theta estimated AFTER lambda;
#   - joint sandwich with multinomial phi.
###############################################################################

run_three_gec <- function(
    donor_g
) {

  ans <- lapply(
    c(
      "ET",
      "HD",
      "CE"
    ),
    function(ent) {

      tryCatch(
        run_GEC_transport_simulation_style(
          d4 = d4,
          donor_g = donor_g,
          ps_obj = ps_obj,
          entropy = ent,
          K = K,
          seed = seed
        ),
        error = function(e) {

          list(
            ATE = NA_real_,
            SE = NA_real_,
            variance = NA_real_,
            CI = c(
              NA_real_,
              NA_real_
            ),
            success = FALSE,
            ESS = NA_real_,
            MaxW = NA_real_,
            CalibrationMax = NA_real_,
            variance_type =
              paste0(
                "failed: ",
                conditionMessage(e)
              )
          )
        }
      )
    }
  )

  names(ans) <- c(
    "ET",
    "HD",
    "CE"
  )

  ans
}


gec_nsw <- run_three_gec(
  2
)

gec_psid <- run_three_gec(
  3
)

gec_cps <- run_three_gec(
  4
)


###############################################################################
# 9. AIPW verification
###############################################################################

aipw_identity_check <- historical_table_AIPW

cat(
  "\n============== AIPW HISTORICAL CHECK ==============\n"
)

print(
  aipw_identity_check
)


###############################################################################
# 10. Main long-form results object
###############################################################################

rows <- list()


add_uw <- function(
    donor,
    x
) {

  make_final_real_row(
    donor = donor,
    method = "Unweighted",
    est = as.numeric(
      x[1]
    ),
    se = as.numeric(
      x[2]
    ),
    benchmark = benchmark,
    variance_type =
      "ordinary unweighted variance"
  )
}


add_fit <- function(
    donor,
    method,
    x
) {

  make_final_real_row(
    donor = donor,
    method = method,
    est = x$ATE,
    se = x$SE,
    benchmark = benchmark,

    ess = if (
      !is.null(x$ESS)
    ) {
      x$ESS
    } else if (
      !is.null(x$ESS0)
    ) {
      x$ESS0
    } else {
      NA_real_
    },

    maxw = if (
      !is.null(x$MaxW)
    ) {
      x$MaxW
    } else if (
      !is.null(x$MaxW0)
    ) {
      x$MaxW0
    } else {
      NA_real_
    },

    cal = if (
      !is.null(x$CalibrationMax)
    ) {
      x$CalibrationMax
    } else {
      NA_real_
    },

    success = if (
      !is.null(x$success)
    ) {
      isTRUE(x$success)
    } else {
      is.finite(x$ATE) &&
        is.finite(x$SE)
    },

    variance_type = if (
      !is.null(x$variance_type)
    ) {
      x$variance_type
    } else {
      ""
    }
  )
}


###############################################################################
# NSW rows
###############################################################################

rows[[length(rows) + 1L]] <-
  add_uw(
    "NSW",
    uw_nsw
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "IPW",
    ipw_nsw
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "EBPS",
    ebps_nsw
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "CBPS",
    cbps_nsw
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "EBCW",
    ebcw_nsw
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "AIPW (LM)",
    aipw_nsw_lm
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "AIPW (GAM)",
    aipw_nsw_gam
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "ET",
    gec_nsw[["ET"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "HD",
    gec_nsw[["HD"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "NSW",
    "CE",
    gec_nsw[["CE"]]
  )


###############################################################################
# PSID rows
###############################################################################

rows[[length(rows) + 1L]] <-
  add_uw(
    "PSID",
    uw_psid
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "IPW",
    ipw_psid
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "EBPS",
    ebps_psid
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "CBPS",
    cbps_psid
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "EBCW",
    ebcw_psid
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "AIPW (LM)",
    aipw_psid_lm
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "AIPW (GAM)",
    aipw_psid_gam
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "ET",
    gec_psid[["ET"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "HD",
    gec_psid[["HD"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "PSID",
    "CE",
    gec_psid[["CE"]]
  )


###############################################################################
# CPS rows
###############################################################################

rows[[length(rows) + 1L]] <-
  add_uw(
    "CPS",
    uw_cps
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "IPW",
    ipw_cps
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "EBPS",
    ebps_cps
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "CBPS",
    cbps_cps
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "EBCW",
    ebcw_cps
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "AIPW (LM)",
    aipw_cps_lm
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "AIPW (GAM)",
    aipw_cps_gam
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "ET",
    gec_cps[["ET"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "HD",
    gec_cps[["HD"]]
  )

rows[[length(rows) + 1L]] <-
  add_fit(
    "CPS",
    "CE",
    gec_cps[["CE"]]
  )


###############################################################################
# Bind and order long-form results
###############################################################################

final_long <-
  dplyr::bind_rows(
    rows
  )

method_order <- c(
  "Unweighted",
  "IPW",
  "EBPS",
  "CBPS",
  "EBCW",
  "AIPW (LM)",
  "AIPW (GAM)",
  "ET",
  "HD",
  "CE"
)

final_long$Method <-
  factor(
    final_long$Method,
    levels = method_order
  )

final_long <-
  final_long |>
  dplyr::mutate(
    Donor = factor(
      Donor,
      levels = c(
        "NSW",
        "PSID",
        "CPS"
      )
    )
  ) |>
  dplyr::arrange(
    Method,
    Donor
  )


###############################################################################
# 10A. Detailed diagnostic print table
###############################################################################

final_long_print <-
  final_long |>
  dplyr::mutate(
    Estimate = round(
      Estimate,
      1
    ),
    EvalDiff = round(
      EvalDiff,
      1
    ),
    SE = round(
      SE,
      1
    ),
    CI_L = round(
      CI_L,
      1
    ),
    CI_U = round(
      CI_U,
      1
    ),
    ESS0 = round(
      ESS0,
      1
    ),
    MaxW0 = round(
      MaxW0,
      3
    ),
    CalibrationMax = signif(
      CalibrationMax,
      5
    )
  )

cat(
  "\n================ DETAILED FINAL RESULTS ================\n"
)

print(
  final_long_print
)


###############################################################################
# 10B. ORIGINAL-PAPER STYLE TABLE
#
# Original paper layout:
#   NSW:  Estimate, SE
#   PSID: Estimate, EB, SE
#   CPS:  Estimate, EB, SE
#
# The published table defines EB relative to the rounded NSW benchmark 1794.
###############################################################################

paper_benchmark <- 1794

get_result <- function(
    donor_name,
    method_name
) {

  z <-
    final_long |>
    dplyr::filter(
      as.character(Donor) == donor_name,
      as.character(Method) == method_name
    )

  if (nrow(z) != 1L) {
    return(
      c(
        Estimate = NA_real_,
        SE = NA_real_
      )
    )
  }

  c(
    Estimate = z$Estimate[1],
    SE = z$SE[1]
  )
}

paper_rows <-
  lapply(
    method_order,
    function(m) {

      nsw <- get_result(
        "NSW",
        m
      )

      psid <- get_result(
        "PSID",
        m
      )

      cps <- get_result(
        "CPS",
        m
      )

      data.frame(
        Estimators = m,
        NSW_Estimates = round(
          nsw["Estimate"]
        ),
        NSW_SE = round(
          nsw["SE"]
        ),
        PSID_Estimates = round(
          psid["Estimate"]
        ),
        PSID_EB = round(
          psid["Estimate"] -
            paper_benchmark
        ),
        PSID_SE = round(
          psid["SE"]
        ),
        CPS_Estimates = round(
          cps["Estimate"]
        ),
        CPS_EB = round(
          cps["Estimate"] -
            paper_benchmark
        ),
        CPS_SE = round(
          cps["SE"]
        ),
        stringsAsFactors = FALSE
      )
    }
  )

final_paper_table <-
  dplyr::bind_rows(
    paper_rows
  )

cat(
  "\n============================================================\n"
)
cat(
  "FINAL LALONDE TABLE -- ORIGINAL PAPER FORMAT\n"
)
cat(
  "============================================================\n\n"
)

print(
  final_paper_table,
  row.names = FALSE
)


###############################################################################
# 10C. Historical seven-row reproduction check
###############################################################################

historical_original_table <- data.frame(
  Estimators = c(
    "Unweighted",
    "IPW",
    "CBPS",
    "EBCW",
    "AIPW (LM)",
    "AIPW (GAM)",
    "ET"
  ),

  NSW_Estimates = c(
    1794,
    1796,
    1636,
    1792,
    1779,
    1780,
    1787
  ),

  NSW_SE = c(
    671,
    673,
    687,
    666,
    673,
    672,
    670
  ),

  PSID_Estimates = c(
    -15205,
    1474,
    1389,
    2316,
    1283,
    1426,
    1655
  ),

  PSID_EB = c(
    -17000,
    -320,
    -405,
    522,
    -511,
    -368,
    -139
  ),

  PSID_SE = c(
    657,
    901,
    886,
    799,
    909,
    893,
    912
  ),

  CPS_Estimates = c(
    -8498,
    1064,
    1288,
    1284,
    1291,
    1683,
    1615
  ),

  CPS_EB = c(
    -10292,
    -730,
    -506,
    -510,
    -503,
    -111,
    -179
  ),

  CPS_SE = c(
    583,
    644,
    641,
    633,
    652,
    656,
    635
  )
)

current_old_methods <-
  final_paper_table |>
  dplyr::filter(
    Estimators %in%
      historical_original_table$Estimators
  )

reproduction_check <-
  historical_original_table |>
  dplyr::left_join(
    current_old_methods,
    by = "Estimators",
    suffix = c(
      "_Paper",
      "_Current"
    )
  )

cat(
  "\n============================================================\n"
)
cat(
  "ORIGINAL PAPER REPRODUCTION CHECK\n"
)
cat(
  "============================================================\n\n"
)

print(
  reproduction_check,
  row.names = FALSE
)


###############################################################################
# 11. GEC diagnostic table
###############################################################################

gec_diagnostics <-
  final_long |>
  dplyr::filter(
    as.character(Method) %in%
      c(
        "ET",
        "HD",
        "CE"
      )
  ) |>
  dplyr::select(
    Donor,
    Method,
    Estimate,
    SE,
    ESS0,
    MaxW0,
    CalibrationMax,
    Success,
    VarianceType
  )

cat(
  "\n================ GEC DIAGNOSTICS ================\n"
)

print(
  gec_diagnostics
)


###############################################################################
# 12. Save outputs
###############################################################################

write.csv(
  final_paper_table,
  "FINAL_LALONDE_TABLE_ORIGINAL_PAPER_FORMAT.csv",
  row.names = FALSE
)

write.csv(
  final_long_print,
  "FINAL_REAL_TRANSPORT_SIMULATION_MATCHED_RESULTS.csv",
  row.names = FALSE
)

write.csv(
  reproduction_check,
  "FINAL_ORIGINAL_PAPER_REPRODUCTION_CHECK.csv",
  row.names = FALSE
)

write.csv(
  aipw_identity_check,
  "FINAL_AIPW_POINT_IDENTITY_CHECK.csv",
  row.names = FALSE
)

write.csv(
  gec_diagnostics,
  "FINAL_GEC_SIMULATION_STYLE_DIAGNOSTICS.csv",
  row.names = FALSE
)


###############################################################################
# Main objects to send back for checking:
#
# final_paper_table
# reproduction_check
# final_long_print
# aipw_identity_check
# gec_diagnostics
###############################################################################
