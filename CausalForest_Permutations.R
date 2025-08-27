library(grf)
library(hstats)

#~~~~~~Must run CausalForest_SHAP.R FIRST~~~~~~#

hstat_overall_original <- data.frame(Features = rownames(interaction_effects$h2_overall$num), Score=interaction_effects$h2_overall$num/interaction_effects$h2_overall$denom)
hstat_pairwise_original <- data.frame(Features = rownames(interaction_effects$h2_pairwise$num), Score=interaction_effects$h2_pairwise$num/interaction_effects$h2_pairwise$denom)

# Define number of permutations
n_perm <- 100

# empty 1-way and 2-way permutations arrays 
oneway_perm <- matrix(NA, nrow = length(interaction_effects_overall$Overall), ncol = n_perm)
rownames(oneway_perm) <- interaction_effects_overall$Feature

twoway_perm <- matrix(NA, nrow = length(interaction_effects_pairwise$Pairwise), ncol = n_perm)
rownames(twoway_perm) <- interaction_effects_pairwise$Pairwise_Features

for (i in 1:n_perm) {
  # Permute outcome
  Y_perm <- sample(cf_med$Y.orig)
  
  # Re-fit causal forest
  cf_perm <- causal_forest(cf_med$X.orig, Y_perm, cf_med$W.orig)
  
  # Calculate hstats
  hstats_perm <- hstats(cf_perm, cf_med$X.orig)
  
  # Fill in 1-way
  perm_oneway <- hstats_perm$h2_overall[1]/ hstats_perm$h2_overall[2]
  oneway_perm[rownames(oneway_perm) %in% names(perm_oneway), i] <- perm_oneway
  
  # Fill in 2-way
  perm_twoway <- hstats_perm$h2_pairwise[1]/hstats_perm$h2_pairwise[2]
  twoway_perm[rownames(twoway_perm) %in% names(perm_twoway), i] <- perm_twoway
}


# Calculate empirical p-values for both
p_values_oneway <- apply(oneway_perm, 1, function(null_dist) {
  (sum(null_dist >= hstat_overall_original$Score) + 1) / (n_perm + 1)
})

p_values_twoway <- apply(twoway_perm, 1, function(null_dist) {
  (sum(null_dist >= hstat_pairwise_original$Score) + 1) / (n_perm + 1)
})


