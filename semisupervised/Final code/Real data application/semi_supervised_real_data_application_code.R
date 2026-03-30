###############################################################################
# SEMI-SUPERVISED REGRESSION UNDER MAR (NHANES 2017–2018)
# - X fully observed (age, sex, race, BMI, SBP, DBP)
# - Y (fasting plasma glucose LBXGLU) missing for most → unlabeled >> labeled
# - Compare:
#       (1) Supervised OLS using labeled only
#       (2) PSSE estimator
#       (3) ET estimator (EM + CVXR + K-fold)
###############################################################################

rm(list = ls())
setwd("C:/Users/moushumi/Desktop/Codes_Reproducibility/Final semi-supervised code")

library(dplyr)
library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(nhanesA)

#source("semi_supervised_methods.R")
source("semi_supervised_real_data_functions.R")

###############################################################################
# 1. DOWNLOAD + MERGE NHANES DATA
###############################################################################
demo <- nhanes("DEMO_J")   # demographics
bmx  <- nhanes("BMX_J")    # body measures
bpx  <- nhanes("BPX_J")    # blood pressure
glu  <- nhanes("GLU_J")    # fasting glucose (Y)

dat <- demo %>%
  transmute(SEQN,
            age  = RIDAGEYR,
            sex  = factor(RIAGENDR, labels = c("Male", "Female")),
            race = factor(RIDRETH3)) %>%
  left_join(bmx %>% dplyr::select(SEQN, BMI = BMXBMI), by = "SEQN") %>%
  left_join(bpx %>% dplyr::select(SEQN, SBP = BPXSY1, DBP = BPXDI1), by = "SEQN") %>%
  left_join(glu %>% dplyr::select(SEQN, Y = LBXGLU), by = "SEQN")


###############################################################################
# 2. KEEP ONLY COMPLETE X (MAR STRUCTURE FOR Y)
###############################################################################
datX <- dat %>%
  filter(!if_any(c(age, sex, race, BMI, SBP, DBP), is.na))

cat("Total N:", nrow(datX), "\n")
cat("Number labeled (Y observed):", sum(!is.na(datX$Y)), "\n")
cat("Number unlabeled (Y missing):", sum(is.na(datX$Y)), "\n")
cat("Fraction unlabeled:", mean(is.na(datX$Y)), "\n")

###############################################################################
# 3. CREATE LABELED / UNLABELED SETS (MAR)
###############################################################################
#labeled   <- datX %>% filter(!is.na(Y))
#unlabeled <- datX %>% filter( is.na(Y))
##############################################################################
#Subsampling of labeled data
#############################################################################
set.seed(123)
labeled   <- datX %>% filter(!is.na(Y)) %>% slice_sample(prop = 0.5)
unlabeled <- datX %>% filter(is.na(Y))
# ------------------------------------------------------------------
# OPTIONAL: SUBSAMPLING OF LABELED DATA (KEEPING THIS “ALIVE”)
# Uncomment to work with a smaller labeled subset:
# set.seed(123)
# labeled <- labeled %>% slice_sample(prop = 0.5)   # e.g., 50% of labeled
# ------------------------------------------------------------------

###############################################################################
# 4. ONE-HOT ENCODING FOR LABELED / UNLABELED
###############################################################################
encode_data <- function(df) {
  Xmm <- model.matrix(~ age + sex + race + BMI + SBP + DBP, data = df)[, -1]
  Xmm <- as.data.frame(Xmm)
  # Clean column names
  colnames(Xmm) <- gsub("[^A-Za-z0-9_]", "_", colnames(Xmm))
  colnames(Xmm) <- gsub("_+", "_", colnames(Xmm))
  colnames(Xmm) <- gsub("_$", "", colnames(Xmm))
  Xmm
}

# Labeled: Y + encoded X
labeled_final <- cbind(
  Y = labeled$Y,
  encode_data(labeled)
)

# Unlabeled: only X (Y is missing)
unlabeled_final <- encode_data(unlabeled)

###############################################################################
# 5. SUPERVISED OLS USING LABELED ONLY (FULL TABLE WITH CI, WIDTH, VARIANCE)
###############################################################################
model.formula <- as.formula(
  paste("Y ~", paste(sprintf("`%s`", colnames(labeled_final)[-1]), collapse = " + "))
)

model.fit <- lm(model.formula, data = as.data.frame(labeled_final))

# Coefficient summary (Estimate, Std. Error, etc.)
S <- summary(model.fit)$coefficients

# Compute variance and 95% CI
Variance <- (S[, "Std. Error"])^2
CI.lower <- S[, "Estimate"] - 1.96 * S[, "Std. Error"]
CI.upper <- S[, "Estimate"] + 1.96 * S[, "Std. Error"]
CI.width <- CI.upper - CI.lower

# ✅ Full supervised table (same structure you wrote)
coef_table <- data.frame(
  Variable = rownames(S),
  Estimate = S[, "Estimate"],
  Variance = Variance,
  StdError = S[, "Std. Error"],
  CI_Lower = CI.lower,
  CI_Upper = CI.upper,
  CI_Width = CI.width,
  row.names = NULL
)

cat("\n--- Supervised OLS (labeled only) ---\n")
print(coef_table)

###############################################################################
# 6. PSSE ESTIMATION (FULL TABLE WITH CI, WIDTH, VARIANCE)
###############################################################################
PSSE_estimate <- PSSE1(
  labelled_data   = as.matrix(labeled_final),
  unlabelled_data = as.matrix(unlabeled_final),
  c1     = NULL,
  type   = "linear",
  tau    = 0,
  alpha  = 1,
  gamma  = 10,
  sd     = TRUE,
  Kfolds = 5
)

Estimate <- as.numeric(PSSE_estimate$Hattheta)
StdError <- as.numeric(PSSE_estimate$sd.of.hattheta)
Variance <- StdError^2

CI_Lower <- Estimate - 1.96 * StdError
CI_Upper <- Estimate + 1.96 * StdError
CI_Width <- CI_Upper - CI_Lower
ARE<-round((as.numeric(S[,"Std. Error"])^2)/(as.numeric(PSSE_estimate$sd.of.hattheta)^2),digits=3)
# Variable names aligned with the supervised model
Variable <- names(coef(lm(model.formula, data = as.data.frame(labeled_final))))

# ✅ Full PSSE table
PSSE_table <- data.frame(
  Variable = Variable,
  Estimate = Estimate,
  StdError = StdError,
  Variance = Variance,
  CI_Lower = CI_Lower,
  CI_Upper = CI_Upper,
  CI_Width = CI_Width,
  ARE=ARE,
  row.names = NULL
)
print(PSSE_table)


################################################################################
#DRESS
################################################################################

DRESS1<-DRESS(labelled_data   = as.matrix(labeled_final),
              unlabelled_data = as.matrix(unlabeled_final),L=1,Kfolds=5)
cat("\n--- PSSE Estimates ---\n")

round((as.numeric(S[,"Std. Error"])^2)/(as.numeric(DRESS1$sd.of.hattheta)^2),digits=3)



Estimate <- as.numeric(DRESS1$Hattheta)
StdError <- as.numeric(DRESS1$sd.of.hattheta)
Variance <- StdError^2

CI_Lower <- Estimate - 1.96 * StdError
CI_Upper <- Estimate + 1.96 * StdError
CI_Width <- CI_Upper - CI_Lower
ARE<-round((as.numeric(S[,"Std. Error"])^2)/(Variance),digits=3)
# Variable names aligned with the supervised model
Variable <- names(coef(lm(model.formula, data = as.data.frame(labeled_final))))

# ✅ Full PSSE table
DRESS_table <- data.frame(
  Variable = Variable,
  Estimate = Estimate,
  StdError = StdError,
  Variance = Variance,
  CI_Lower = CI_Lower,
  CI_Upper = CI_Upper,
  CI_Width = CI_Width,
  ARE=ARE,
  row.names = NULL
)
print(DRESS_table)


# ---- DRESS -> Overleaf block (2 decimals) ----
v <- DRESS_table$Variable
to2 <- function(x) sprintf("%.2f", x)

est <- paste(to2(DRESS_table$Estimate),  collapse = " & ")
se  <- paste(to2(DRESS_table$StdError),  collapse = " & ")
ciw <- paste(to2(DRESS_table$CI_Width),  collapse = " & ")
are <- paste(to2(DRESS_table$ARE),       collapse = " & ")

cat(
  "% ================= DRESS ==================\n",
  "\\multirow{4}{*}{\\shortstack{\\textbf{DRESS}}}\n",
  "& Est  & ", est, " \\\\\n",
  "& SE   & ", se,  " \\\\\n",
  "& CIW  & ", ciw, " \\\\\n",
  "& ARE  & ", are, " \\\\\n",
  "\\hline\n",
  sep = ""
)



################PI########################################




###############################################################################
# 7. BUILD SEMI-SUPERVISED FULL DATA WITH D, pi.hat, ID
###############################################################################
# Remove SEQN from both parts
labeled_noseqn   <- labeled   %>% dplyr::select(-SEQN)
unlabeled_noseqn <- unlabeled %>% dplyr::select(-SEQN)

# Stack labeled and unlabeled
data_full_real <- as.data.frame(rbind(
  cbind(labeled_noseqn,   D = 1),
  cbind(unlabeled_noseqn, D = 0)
))
# Reorder columns so Y is first
data_full_real <- data_full_real %>%
  dplyr::select(Y, everything())

glm.model.formula <- as.formula(
  paste("D ~", paste(colnames(data_full_real)[!colnames(data_full_real) %in% c("Y", "D")],
                     collapse = " + "))
)

pi.hat <- fitted(glm(glm.model.formula,
                     family = binomial(link = "logit"),
                     data   = as.data.frame(data_full_real[, -c(1)])))

data_full_real$pi.hat <- pi.hat
data_full_real$ID     <- seq_len(nrow(data_full_real))

###############################################################################
# 8. ET ESTIMATOR (EM + CVXR + K-FOLD), note that ET and HD may give you slightly different results than in the original table due to cross fitting 
###############################################################################
# (uses your existing function definitions:
#  variance_theta_diag, build_SU_folds, k_fold_function,
#  estimate_theta_EM_kfold_CVXR_ET, etc.)
# Make sure those functions are defined above or in sourced files.

ET <- estimate_theta_EM_kfold_CVXR_ET(
  th        = as.numeric(coef(model.fit)),
  data_full = data_full_real,
  K         = 2,
  seed      = 2025,
  max.iter  = 50,
  eps       = 1e-4
)

cat("\n--- ET (EM + CVXR) Estimated Table ---\n")
print(ET$estimate.table)

ARE<-round((as.numeric(S[,"Std. Error"])^2)/as.numeric(ET$estimate.table$Variance),digits=3)
cbind(ET$estimate.table,ARE)



###############################################################################
# 8. HD ESTIMATOR (EM + CVXR + K-FOLD)
###############################################################################
# (uses your existing function definitions:
#  variance_theta_diag, build_SU_folds, k_fold_function,
#  estimate_theta_EM_kfold_CVXR_ET, etc.)
# Make sure those functions are defined above or in sourced files.

HD <- estimate_theta_EM_kfold_CVXR_HD(
  th        = as.numeric(coef(model.fit)),
  data_full = data_full_real,
  K         = 2,
  seed      = 2025,
  max.iter  = 500,
  eps       = 1e-4
)

cat("\n--- HD (EM + CVXR) Estimated Table ---\n")
print(HD$estimate.table)

ARE<-round((as.numeric(S[,"Std. Error"])^2)/as.numeric(HD$estimate.table$Variance),digits=3)
cbind(HD$estimate.table,ARE)









