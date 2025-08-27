options(grf.num.threads = parallel::detectCores())
n_cores <- parallel::detectCores() - 1

library(grf)
library(hstats)
library(tidyverse)
library(patchwork)
library(kernelshap)
library(shapviz)
library(pdp)
library(parallel)
library(scico)


#---LOAD AND FILTER---#
df <- read.csv("/Users/bobkohler/Desktop/Manuscripts/Machine Learning/Prebaseline Effect Manuscript Materials/vars_all_recode.csv") %>%
  mutate(
    Y_pre = Average_PreBaseline,
    Y_trt = Average_Treatment,
    Y_diff = Y_trt - Y_pre,
    type_treatment = factor(type_treatment)) %>%
  drop_na(type_treatment, Y_pre, Y_trt)

clusters_all <- as.integer(factor(df$study))

cov_df <- df %>%
  select(-subject_id, -X,
         -ends_with(".Final"), 
         -starts_with("Treatment_"), -starts_with("PreBaseline_"),
         -Average_PreBaseline, -Average_Treatment, -Average_Baseline, -pretx_SDU, 
         -SITE, -study, -YearlyIncome,
         -Y_pre, -Y_trt, -Y_diff, -type_treatment, -treatment,
         -Baseline_Change_Group, -Baseline_Change)

cov_df_cc <- drop_na(cov_df)

X_full <- model.matrix(~ . -1, data = cov_df_cc, na.action = na.pass)
kept_idx <- as.integer(rownames(X_full))

Y_trt <- df$Y_trt[kept_idx]
W_pre <- df$Y_pre[kept_idx]
clusters <- clusters_all[kept_idx]


# RENAME Variables
rename_vars <- c(
  "GGT Levels"                 = "baseline_GGT",
  "Male"                           = "GenderMale",
  "Female"                         = "GenderFemale")

colnames(X_full) <- recode(colnames(X_full), !!!rename_vars)



# CAUSAL FOREST FIT
set.seed(42)
cf_med <- causal_forest(
  X = X_full,
  Y = Y_trt,
  W = W_pre,
  clusters = clusters,
  num.trees = 5000,
  equalize.cluster.weights = TRUE,
  tune.parameters = "all",
  tune.num.reps = 100,
  tune.num.trees = 100,
  num.threads = n_cores)

saveRDS(cf_med, "/Users/bobkohler/Desktop/Manuscripts/Machine Learning/Prebaseline Effect Manuscript Materials/cf_output.rds")

test_calibration(cf_med, vcov.type = "HC3")


# SHAP 
pred_fun <- function(object, newdata) {
  predict(object, newdata)$predictions
}

set.seed(42)
ks_train <- kernelshap(cf_med, bg_n = 100, X = X_full, pred_fun = pred_fun)
saveRDS(ks_train, "/Users/bobkohler/Desktop/Manuscripts/Machine Learning/Prebaseline Effect Manuscript Materials/shap_output_train.rds")

# CATE
pred <- predict(cf_med, estimate.variance = TRUE)
cates <- pred$predictions
cate_se <- sqrt(pred$variance.estimates)
z_scores <- cates / cate_se
p_values <- 2 * (1 - pnorm(abs(z_scores)))

cate_results <- data.frame(
  subject_id = 1:length(cates),
  cate_value = cates,
  standard_error = cate_se,
  zscore = z_scores,
  pvalue = p_values)
cate_significant_results <- cate_results %>% filter(pvalue < 0.05)


# CF Variable Importance
imp_cf <- sort(setNames(variable_importance(cf_med), colnames(X_full)))
imp_df <- data.frame(Feature = names(imp_cf), Value = as.numeric(imp_cf)) %>%
  mutate(FeatureLabel = recode(Feature, !!!rename_vars))


p_cf <- ggplot(imp_df[11:25,], aes(x = reorder(FeatureLabel, Value), y = Value)) +
  geom_bar(stat = "identity", fill = "dodgerblue", color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, y = "Causal Forest Importance", title = "Causal Forest Feature Importance")
p_cf

# SHAP Variable Importance
sv <- shapviz(ks_train)
p_shap <- sv_importance(sv) +
  geom_bar(stat = "identity", fill = "orange", color = "black") +
  labs(x = "Average SHAP", y = NULL, title = "SHAP Feature Importance") +
  theme_minimal()

# SHAP dependence plot
sv_dependence(sv, v = "baseline_ALT")


# hstats interactions test 
interaction_effects <- hstats(cf_med, X = X_full, pred_fun = pred_fun)
saveRDS(interaction_effects, "/Users/bobkohler/Desktop/Manuscripts/Machine Learning/Prebaseline Effect Manuscript Materials/interaction_effects.rds")


interaction_summary <- summary(interaction_effects)
# Overall interactions
interaction_overall <- data.frame(
  Feature = names(interaction_summary$h2_overall$M[,1]),
  Overall = as.numeric(interaction_summary$h2_overall$M[,1]))

ggplot(interaction_overall[1:15,], aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "#C4A77D", color = "black") +
  coord_flip() +
  labs(title = "Overall Interaction Strength") +
  theme_minimal()

# Pairwise interactions
interaction_pairwise <- data.frame(
  Feature = names(interaction_summary$h2_pairwise$M[,1]),
  Overall = as.numeric(interaction_summary$h2_pairwise$M[,1]))

ggplot(interaction_pairwise, aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "#C4A77D", color = "black") +
  coord_flip() +
  labs(title = "Pairwise Interaction Strength", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        legend.position = "none")


# Partial Dependence Plots
cf_med <- readRDS("/Users/bobkohler/Desktop/Manuscripts/Machine Learning/Prebaseline Effect Manuscript Materials/cf_output.rds")
grf::test_calibration(cf_med)
tau.hat <- predict(cf_med)$predictions


X_full <- as.data.frame(X_full)
summary(lm(Y_trt ~ W_pre + tau.hat))
plot(pred$predictions, resid(lm(Y_trt ~ W_pre)), 
     pch = 19, col = "blue", 
     xlab = "CATEs (tau.hat)", ylab = "Residual Treatment Drinking")
abline(lm(resid(lm(Y_trt ~ W_pre)) ~ pred$predictions), col = "red")


pd1 <- partial_dep(cf_med, v = "baseline_CREATINE", X = X_full, BY = X_full[, "baseline_dep_anx_composite"], pred_fun = pred_fun)
plot(pd1) + labs(title = "Creatine Levels vs. Depression & Anxiety Score", y = "Average Treatment Effect", x = "Creatine Levels", color = "Depression & Anxiety \nScore Quantiles")


