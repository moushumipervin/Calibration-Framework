################################################################################
# A general M-estimation theory in the semi-supervised framework
# Part 4: Reproducing the subtable of Table 6 in Section 4.1
#
# This script:
#   1. Loads required packages and source files
#   2. Sets simulation parameters
#   3. Runs Monte Carlo replications
#   4. Computes summary measures (bias, SE, ARE)
#   5. Produces boxplots and a ggplot figure
################################################################################

rm(list = ls())

################################################################################
# 1. Required packages
################################################################################
library(MASS)
library(mgcv)
library(CVXR)
library(caret)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)

################################################################################
# 2. Working directory and source files
################################################################################

source("semi_supervised_methods.R")
source("Data_generation.R")
source("SupervisedEstimation.R")
source("semi_supervised_real_data_functions (1).R")
source("variance_estimation.R")

################################################################################
# 3. Global parameters
################################################################################
n <- 1000                 # labeled sample size
N <- 1500                 # unlabeled sample size
p <- 4                    # number of predictors
rep <- 1000               # number of Monte Carlo replications
polyOrder <- 4
OR <- 2                   # OR =1 indicates OR1 model, otherwise OR2 model
MAR <- 1                  #MAR =1 indicates MAR mechanism, otherwise it indicates MCAR mechanism
q <- 5                    # number of folds

###############################################################################
# TRUE TARGET PARAMETER
###############################################################################

set.seed(1234)

n_target <- 10^7

alpha0 <- 1
alpha1 <- rep(1, p)
alpha2 <- rep(1, p)

X_target <- MASS::mvrnorm(
  n = n_target,
  mu = rep(0, p),
  Sigma = diag(p)
)

if (OR == 1) {
  
  Y_target <-
    alpha0 +
    X_target %*% alpha1 +
    rnorm(n_target, 0, 1)
  
} else {
  
  Y_target <-
    alpha0 +
    X_target %*% alpha1 +
    (
      X_target^3 -
        X_target^2 +
        exp(X_target)
    ) %*%
    alpha2 +
    rnorm(n_target, 0, 2)
}

target_parameter <-
  as.numeric(
    coef(
      lm(
        Y_target ~ X_target
      )
    )
  )

parameter_names <-
  c(
    "(Intercept)",
    paste0("X", seq_len(p))
  )
###############################################################################
# Storage
###############################################################################

simulation_results <- vector(
  "list",
  rep
)
################################################################################
# 6. Monte Carlo replications
################################################################################
###############################################################################
# MONTE CARLO SIMULATION
###############################################################################

for (k in seq_len(rep)) {
  
  message(
    "Replication ",
    k,
    "/",
    rep
  )
  
  set.seed(
    k + 20220122
  )
  
  one_rep <- list()
  
  ###########################################################################
  # 1. Generate data
  ###########################################################################
  
  DesiredData <-
    GenerateData(
      n = n,
      N = N,
      p = p,
      OR = OR,
      MAR = MAR
    )
  
  data_labelled <-
    DesiredData$Data.labelled
  
  data_unlabelled <-
    DesiredData$Data.unlabelled
  
  data_full <-
    DesiredData$data_full
  
  
  ###########################################################################
  # 2. Initial OLS
  ###########################################################################
  
  formula_lm <-
    as.formula(
      paste0(
        "Y ~ ",
        paste0(
          colnames(
            data_labelled[, -1]
          ),
          collapse = " + "
        )
      )
    )
  
  fit_initial <-
    lm(
      formula_lm,
      data =
        as.data.frame(
          data_labelled
        )
    )
  
  theta_init <-
    as.numeric(
      coef(fit_initial)
    )
  
  
  ###########################################################################
  # 3. SUPERVISED
  ###########################################################################
  
  fit_sup <-
    tryCatch(
      SupervisedEst(
        data_labelled
      ),
      error =
        function(e) NULL
    )
  
  if (!is.null(fit_sup)) {
    
    sup_est <-
      as.numeric(
        fit_sup$Est.coef
      )
    
    # If SupervisedEst does not return SE,
    # use standard OLS SE from the labeled fit.
    sup_se <-
      as.numeric(
        summary(fit_initial)$coefficients[
          ,
          "Std. Error"
        ]
      )
    
    one_rep[["Supervised"]] <-
      make_simulation_rows(
        rep_id = k,
        method = "Supervised",
        estimate = sup_est,
        truth = target_parameter,
        se = sup_se,
        parameter_names =
          parameter_names
      )
  }
  ###########################################################################
  # 4. PSSE
  ###########################################################################
  
  fit_PSSE <-
    tryCatch(
      PSSE1(
        data_labelled,
        data_unlabelled,
        type = "linear",
        sd = TRUE,
        alpha = polyOrder
      ),
      error = function(e) {
        message(
          "PSSE failed in replication ",
          k,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
  
  if (!is.null(fit_PSSE)) {
    
    PSSE_est <-
      as.numeric(
        fit_PSSE$Hattheta
      )
    
    PSSE_se <-
      if (!is.null(fit_PSSE$sd.of.hattheta)) {
        as.numeric(
          fit_PSSE$sd.of.hattheta
        )
      } else {
        rep(
          NA_real_,
          length(target_parameter)
        )
      }
    
    one_rep[["PSSE"]] <-
      make_simulation_rows(
        rep_id = k,
        method = "PSSE",
        estimate = PSSE_est,
        truth = target_parameter,
        se = PSSE_se,
        parameter_names = parameter_names
      )
  }
  
  
  
  
  ###########################################################################
  # 5. PI
  ###########################################################################
  
  fit_PI <-
    tryCatch(
      PI(
        data_labelled,
        data_unlabelled
      ),
      error =
        function(e) NULL
    )
  
  if (!is.null(fit_PI)) {
    
    PI_est <-
      as.numeric(
        fit_PI$Hattheta
      )
    
    PI_se <-
      if (!is.null(
        fit_PI$sd.of.hattheta
      )) {
        as.numeric(
          fit_PI$sd.of.hattheta
        )
      } else {
        rep(
          NA_real_,
          length(
            target_parameter
          )
        )
      }
    
    one_rep[["PI"]] <-
      make_simulation_rows(
        rep_id = k,
        method = "PI",
        estimate = PI_est,
        truth = target_parameter,
        se = PI_se,
        parameter_names =
          parameter_names
      )
  }
  
  
  ###########################################################################
  # 6. EASE
  ###########################################################################
  
  fit_EASE <-
    tryCatch(
      EASE(
        data_labelled,
        data_unlabelled,
        K = 2,
        H = 2,
        r = 2
      ),
      error =
        function(e) NULL
    )
  
  if (!is.null(fit_EASE)) {
    
    EASE_est <-
      as.numeric(
        fit_EASE$Hattheta
      )
    
    EASE_se <-
      if (!is.null(
        fit_EASE$sd.of.hattheta
      )) {
        as.numeric(
          fit_EASE$sd.of.hattheta
        )
      } else {
        rep(
          NA_real_,
          length(
            target_parameter
          )
        )
      }
    
    one_rep[["EASE"]] <-
      make_simulation_rows(
        rep_id = k,
        method = "EASE",
        estimate = EASE_est,
        truth = target_parameter,
        se = EASE_se,
        parameter_names =
          parameter_names
      )
  }
  
  
  ###########################################################################
  # 7. DRESS
  ###########################################################################
  
  fit_DRESS <-
    tryCatch(
      DRESS1(
        data_labelled,
        data_unlabelled,
        L = polyOrder,
        sd = TRUE,
        Kfolds = q
      ),
      error =
        function(e) NULL
    )
  
  if (!is.null(fit_DRESS)) {
    
    DRESS_est <-
      as.numeric(
        fit_DRESS$Hattheta
      )
    
    DRESS_se <-
      if (!is.null(
        fit_DRESS$sd.of.hattheta
      )) {
        as.numeric(
          fit_DRESS$sd.of.hattheta
        )
      } else {
        rep(
          NA_real_,
          length(
            target_parameter
          )
        )
      }
    
    one_rep[["DRESS"]] <-
      make_simulation_rows(
        rep_id = k,
        method = "DRESS",
        estimate = DRESS_est,
        truth = target_parameter,
        se = DRESS_se,
        parameter_names =
          parameter_names
      )
  }
  
  
  ###########################################################################
  # 8. ET
  ###########################################################################
  ET_fit <- tryCatch(
    
    estimate_theta_EM_kfold_dual_ET(
      th = theta_init,
      data_full = data_full,
      K = 5,
      seed = k + 20220122,
      max.iter = 500,
      eps = 1e-4,
      damping = 0.1
    ),
    
    error = function(e) {
      message(
        "ET estimator failed in replication ",
        k,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  
  one_rep[["ET"]] <-
    process_gec_simulation(
      fit_gec = ET_fit,
      entropy = "ET",
      rep_id = k,
      target_parameter =
        target_parameter,
      parameter_names =
        parameter_names
    )
  
  
    
  
  ###########################################################################
  # 9. HD
  ###########################################################################
  
  HD_fit <- tryCatch(
    
    estimate_theta_EM_kfold_dual_HD(
      th = theta_init,
      data_full = data_full,
      K = 5,
      seed = k + 20220122,
      max.iter = 500,
      eps = 1e-4,
      damping = 0.1
    ),
    
    error = function(e) {
      message(
        "HD estimator failed in replication ",
        k,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  
  one_rep[["HD"]] <-
    process_gec_simulation(
      fit_gec = HD_fit,
      entropy = "HD",
      rep_id = k,
      target_parameter =
        target_parameter,
      parameter_names =
        parameter_names
    )
  
  
  ###########################################################################
  # 10. CE
  ###########################################################################
  
  CE_fit <- tryCatch(
    
    estimate_theta_EM_kfold_CE(
      th = theta_init,
      data_full = data_full,
      K = 5,
      seed = k + 20220122,
      max.iter = 500,
      eps = 1e-4,
      damping = 0.1
    ),
    
    error = function(e) {
      message(
        "CE estimator failed in replication ",
        k,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  
  one_rep[["CE"]] <-
    process_gec_simulation(
      fit_gec = CE_fit,
      entropy = "CE",
      rep_id = k,
      target_parameter =
        target_parameter,
      parameter_names =
        parameter_names
    )
  
  ###########################################################################
  # Save replication
  ###########################################################################
  simulation_results[[k]] <-
    dplyr::bind_rows(
      one_rep
    )
  
}

###############################################################################
# MONTE CARLO SUMMARY
###############################################################################

simulation_results_df <-
  dplyr::bind_rows(
    simulation_results
  )

simulation_summary <-
  simulation_results_df %>%
  
  dplyr::group_by(
    method,
    parameter
  ) %>%
  
  dplyr::summarise(
    
    N_success =
      sum(
        is.finite(estimate)
      ),
    
    Mean_Estimate =
      mean(
        estimate,
        na.rm = TRUE
      ),
    
    Bias =
      mean(
        estimate - truth,
        na.rm = TRUE
      ),
    
    MC_SD =
      sd(
        estimate,
        na.rm = TRUE
      ),
    
    RMSE =
      sqrt(
        mean(
          (estimate - truth)^2,
          na.rm = TRUE
        )
      ),
    
    Avg_SE =
      mean(
        analytic_se,
        na.rm = TRUE
      ),
    
    SE_MC_Ratio =
      Avg_SE / MC_SD,
    
    Coverage =
      mean(
        covered,
        na.rm = TRUE
      ),
    
    Avg_CI_Width =
      mean(
        ci_width,
        na.rm = TRUE
      ),
    
    Mean_ESS =
      if (
        all(is.na(ESS))
      ) {
        NA_real_
      } else {
        mean(
          ESS,
          na.rm = TRUE
        )
      },
    
    Mean_Max_Weight =
      if (
        all(is.na(max_weight))
      ) {
        NA_real_
      } else {
        mean(
          max_weight,
          na.rm = TRUE
        )
      },
    
    Max_Weight_Overall =
      if (
        all(is.na(max_weight))
      ) {
        NA_real_
      } else {
        max(
          max_weight,
          na.rm = TRUE
        )
      },
    
    .groups = "drop"
  ) %>%
  
  dplyr::group_by(parameter) %>%
  
  dplyr::mutate(
    
    Supervised_MC_SD =
      MC_SD[method == "Supervised"][1],
    
    ARE =
      (Supervised_MC_SD^2) /
      (MC_SD^2)
    
  ) %>%
  
  dplyr::ungroup()

simulation_summary <-
  simulation_summary %>%
  dplyr::mutate(
    method = factor(
      method,
      levels = c(
        "Supervised",
        "PI",
        "EASE",
        "DRESS",
        "PSSE",
        "ET",
        "HD",
        "CE"
      )
    ),
    parameter = factor(
      parameter,
      levels = c(
        "(Intercept)",
        "X1",
        "X2",
        "X3",
        "X4"
      )
    )
  ) %>%
  dplyr::arrange(
    method,
    parameter
  )
View(simulation_summary)

failure_summary <-
  simulation_results_df  %>%
  dplyr::group_by(
    method,
    parameter
  ) %>%
  dplyr::summarise(
    N_total = rep,
    
    N_success =
      sum(
        is.finite(
          estimate
        )
      ),
    
    Failure_Rate =
      1 -
      N_success /
      N_total,
    
    .groups = "drop"
  )

simulation_summary <-
  simulation_summary %>%
  dplyr::left_join(
    failure_summary,
    by = c(
      "method",
      "parameter"
    )
  )

###############################################################################
#Boxplot
###############################################################################
simulation_results <-
  dplyr::bind_rows(simulation_results)

simulation_results <- simulation_results %>%
  dplyr::mutate(
    beta_idx = match(parameter, parameter_names) - 1,
    
    beta_lab = factor(
      paste0("beta[", beta_idx, "]"),
      levels = paste0("beta[", 0:(length(parameter_names) - 1), "]")
    )
  )

simulation_results$method <- factor(
  simulation_results$method,
  levels = c(
    "Supervised",
    "PI",
    "EASE",
    "DRESS",
    "PSSE",
    "ET",
    "HD",
    "CE"
  )
)

target_df <- data.frame(
  beta_lab = factor(
    paste0(
      "beta[",
      0:(length(target_parameter) - 1),
      "]"
    ),
    levels = paste0(
      "beta[",
      0:(length(target_parameter) - 1),
      "]"
    )
  ),
  
  target = target_parameter
)
library(ggplot2)
library(dplyr)

p <- ggplot(
  simulation_results,
  aes(
    x = method,
    y = estimate
  )
) +
  
  geom_boxplot(
    outlier.size = 0.5,
    fill = "grey80",
    na.rm = TRUE
  ) +
  
  geom_hline(
    data = target_df,
    aes(
      yintercept = target
    ),
    color = "red",
    linewidth = 0.6
  ) +
  
  facet_wrap(
    ~ beta_lab,
    nrow = 1,
    labeller = label_parsed,
    scales = "free_y"
  ) +
  
  labs(
    x = "Method",
    y = "Point estimate"
  ) +
  
  theme_bw() +
  
  theme(
    text = element_text(size = 12),
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    
    strip.text =
      element_text(size = 11)
  )

ggsave(
  filename = "linreg_estimates_gg.pdf",
  plot = p,
  width = 15,
  height = 4
)
plot_data <- simulation_results %>%
  filter(
    method %in%
      c(
        "Supervised",
        "PSSE",
        "DRESS",
        "ET",
        "HD",
        "CE"
      )
  )


ggplot(
  plot_data,
  aes(
    x = method,
    y = estimate
  )
) +
  geom_boxplot(
    fill = "grey80",
    outlier.size = 0.5
  ) +
  geom_hline(
    data = target_df,
    aes(yintercept = target),
    color = "red"
  ) +
  facet_wrap(
    ~ beta_lab,
    nrow = 1,
    labeller = label_parsed,
    scales = "free_y"
  ) +
  theme_bw()

################################################################################
# 11. Final summary table
################################################################################
results_total <- round(
  cbind(
    target_parameter,
    
    colMeans(results_supervised, na.rm = TRUE) - target_parameter,
    apply(results_supervised, 2, sd, na.rm = TRUE),
    
    colMeans(results_PI, na.rm = TRUE) - target_parameter,
    apply(results_PI, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_PI, 2, var, na.rm = TRUE),
    
    colMeans(results_EASE, na.rm = TRUE) - target_parameter,
    apply(results_EASE, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_EASE, 2, var, na.rm = TRUE),
    
    colMeans(results_DRESS, na.rm = TRUE) - target_parameter,
    apply(results_DRESS, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_DRESS, 2, var, na.rm = TRUE),
    
    colMeans(results_PSSE, na.rm = TRUE) - target_parameter,
    apply(results_PSSE, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_PSSE, 2, var, na.rm = TRUE),
    
    colMeans(results_ET, na.rm = TRUE) - target_parameter,
    apply(results_ET, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_ET, 2, var, na.rm = TRUE),
    
    colMeans(results_HD, na.rm = TRUE) - target_parameter,
    apply(results_HD, 2, sd, na.rm = TRUE),
    apply(results_supervised, 2, var, na.rm = TRUE) / apply(results_HD, 2, var, na.rm = TRUE)
  ),
  2
)

colnames(results_total) <- c(
  "real_value",
  "bias_sup", "SE_sup",
  "bias_PI", "SE_PI", "ARE_PI",
  "bias_EASE", "SE_EASE", "ARE_EASE",
  "bias_DRESS", "SE_DRESS", "ARE_DRESS",
  "bias_PSSE", "SE_PSSE", "ARE_PSSE",
  "bias_ET", "SE_ET", "ARE_ET",
  "bias_HD", "SE_HD", "ARE_HD"
)

results_total <- as.data.frame(results_total)
print(results_total)

################################################################################
# 12. Save results
################################################################################
saveRDS(results_EASE, "results_EASE_n1000N1000_OM1PS1.RDS")
# saveRDS(results_total, "Final_summary_table.RDS")