##50yrs ssp245 model set 1
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._4_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##50yrs ssp245 model set 2
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._6_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._6_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##50yrs ssp245 model set 3
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._0_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._0_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##50yrs ssp245 model set 4
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._8_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._8_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##50yrs ssp585 model set 1
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._4_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##50yrs ssp585 model set 2
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._2_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._2_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")



##50yrs ssp585 model set 3
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._7_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._7_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")



##50yrs ssp585 model set 4
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._1_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._1_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##Alltime ssp245 model set 1
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._4_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##Alltime ssp245 model set 2
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._8_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._8_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##Alltime ssp245 model set 3
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._5_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._5_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##Alltime ssp245 model set 4
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._0_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._0_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##Alltime ssp585 model set 1
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._8_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._8_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")


##Alltime ssp585 model set 2
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._1_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._1_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")



##Alltime ssp585 model set 3
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._2_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._2_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")



##Alltime ssp585 model set 4
# ==== Packages ====
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(readr)
})

# ==== 0) INPUTS: edit these ====
presence_csv   <- "Plectranthus_esculentus_N.E.Br._5_samplePredictions.csv"
background_csv <- "Plectranthus_esculentus_N.E.Br._5_backgroundPredictions.csv"

# If you know the prediction column name, set it here; otherwise leave NULL to auto-detect
prediction_col <- NULL   # e.g., "logistic" or "cloglog" or "pred"

# ==== 1) Helpers ====
pick_prediction_col <- function(df, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(df)) return(user_col)
  # Heuristic: choose the last numeric column (typical in MaxEnt outputs)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) == 0) stop("No numeric prediction column found. Set `prediction_col` explicitly.")
  tail(num_cols, 1)
}

confusion_from_threshold <- function(scores, labels01, thr) {
  pred_bin <- ifelse(scores >= thr, 1, 0)
  TP <- sum(pred_bin == 1 & labels01 == 1)
  FP <- sum(pred_bin == 1 & labels01 == 0)
  FN <- sum(pred_bin == 0 & labels01 == 1)
  TN <- sum(pred_bin == 0 & labels01 == 0)
  sens <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  tss  <- sens + spec - 1
  omission <- if ((TP + FN) == 0) NA else FN / (TP + FN)
  list(TP=TP, FP=FP, FN=FN, TN=TN, sensitivity=sens, specificity=spec,
       TSS=tss, omission=omission, thr=thr, pred_bin=pred_bin)
}

print_conf_matrix <- function(TP, FP, FN, TN) {
  m <- matrix(c(TN, FN, FP, TP), nrow = 2, byrow = TRUE,
              dimnames = list(Predicted = c("0","1"), Actual = c("0","1")))
  print(m)
}

# ==== 2) Load data ====
pres <- read_csv(presence_csv, show_col_types = FALSE)
back <- read_csv(background_csv, show_col_types = FALSE)

col_name <- pick_prediction_col(pres, prediction_col)
if (!(col_name %in% names(back))) {
  # try to auto-pick for background too
  col_name_bg <- pick_prediction_col(back, NULL)
  if (col_name_bg != col_name) {
    message(sprintf("Note: using '%s' for presence and '%s' for background.", col_name, col_name_bg))
    col_name <- col_name  # presence column fixed
    # We'll align by using separate picks
  }
}

scores_presence   <- pres[[pick_prediction_col(pres, prediction_col)]]
scores_background <- back[[pick_prediction_col(back, prediction_col)]]

# Labels: 1 = presence, 0 = background
labels <- c(rep(1, length(scores_presence)),
            rep(0, length(scores_background)))
scores <- c(scores_presence, scores_background)

# ==== 3) ROC & AUC (test-style on pooled set here; for CV, do per-fold) ====
roc_obj <- roc(response = labels, predictor = scores, levels = c(0,1), direction = ">", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# ---- Best threshold = Youden's J ----
youden <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_youden <- youden$threshold

# ---- Alternative thresholds ----
# Equal Sensitivity & Specificity
eq_ss <- coords(roc_obj, x = "best", best.method = "closest.topleft", ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
thr_eqss <- eq_ss$threshold

# 10th percentile training presence (commonly reported in ENM)
thr_p10 <- as.numeric(quantile(scores_presence, probs = 0.10, na.rm = TRUE))

# ==== 4) Metrics at each threshold ====
met_youden <- confusion_from_threshold(scores, labels, thr_youden)
met_eqss   <- confusion_from_threshold(scores, labels, thr_eqss)
met_p10    <- confusion_from_threshold(scores, labels, thr_p10)

# ==== 5) Print summary ====
cat("=== ROC Summary ===\n")
cat(sprintf("AUC = %.4f\n\n", auc_val))

cat("=== Thresholds Considered ===\n")
cat(sprintf("Youden's J threshold            = %.6f\n", met_youden$thr))
cat(sprintf("Equal Sens & Spec threshold     = %.6f\n", met_eqss$thr))
cat(sprintf("10th percentile presence (P10)  = %.6f\n\n", met_p10$thr))

report_block <- function(name, met) {
  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("Sensitivity = %.4f | Specificity = %.4f | TSS = %.4f\n",
              met$sensitivity, met$specificity, met$TSS))
  cat(sprintf("Omission rate = %.4f\n", met$omission))
  cat("Confusion matrix (rows = Predicted, cols = Actual):\n")
  print_conf_matrix(TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN)
  cat("\n")
}

report_block("At Youden's J", met_youden)
report_block("At Equal Sens & Spec", met_eqss)
report_block("At 10th percentile presence", met_p10)

# ==== 6) (Optional) If you need a single 'best' model threshold → usually Youden's J ====
cat(">>> Recommended (commonly used): Youden’s J threshold.\n")
