library(tidyr)
library(stringr)

extract_all_terms <- function(formula) {
  tt <- terms(formula)
  term_labels <- attr(tt, "term.labels")

  # Split random effect terms like (1|f_yr) into variable names
  random_terms <- grep("\\|", term_labels, value = TRUE)
  random_vars <- unlist(
    lapply(random_terms, function(x) {
      x <- gsub("^.*\\|", "", x)  # keep only RHS of |
      strsplit(x, ":|/")[[1]]     # split nested or crossed
    })
  )

  fixed_terms <- setdiff(term_labels, random_terms)
  all_terms <- unique(c(fixed_terms, random_vars))
  all_terms
}

split_group_levels <- function(proc, formulas, group_factor, estimates, remove_original = TRUE) {
  # Get the formula and terms
  form <- formulas[[proc]]
  terms <- attr(terms(form), "term.labels")

  # Get the levels of the group factor (e.g., "2013.Female", "2013:Female", or both)
  levels_grp <- levels(group_factor[[proc]])

  # Split by "." and ":" to capture fixed and random effects
  split_vals <- strsplit(levels_grp, split = "[.:]")

  # Get unique variable names from formula (use all.vars to support nested interactions)
  vars <- all.vars(formula(paste("~", paste(terms, collapse = "+"))))

  # Build data frame
  df <- as.data.frame(do.call(rbind, split_vals), stringsAsFactors = FALSE)
  names(df) <- make.unique(vars, sep = "_")  # Avoid duplicate column names

  df$estimate <- estimates
  if (!remove_original) {
    df$.grp <- levels_grp
  }
  df
}

split_group_levels("Nsuper", formulas, data$group_factor, rep$Nsuper_vec) %>%
  ggplot(aes(x = as.numeric(as.character(f_yr)), y = estimate, color = f_sex)) +
  geom_line()


library(dplyr)
library(ggplot2)

# assume you’ve already computed:
#   uTot_g     = observed total detections per group (length K)
#   lambda_k   = fitted Poisson means per group (length K)

df1 <- tibble(
  group = factor(1:rep$K),
  observed = rep$uTot_g,
  fitted   = rep$lambda_k
)

ggplot(df1, aes(x = observed, y = fitted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title    = "Poisson fit: observed vs fitted total detections",
    x        = "Observed uTot_g",
    y        = "Fitted λ_k"
  ) +
  theme_minimal()

# u_k_t is K×s matrix of observed counts
# multP_k_t is K×s matrix of fitted probs

df2 <- bind_rows(lapply(1:rep$K, function(k) {
  tibble(
    group    = factor(k),
    week     = 1:rep$s,
    observed = rep$u_k_t2[k, ],
    fitted   = rep$multP_k_t[k, ]
  )
}), .id = "dummy") %>% select(-dummy)

# plot observed vs fitted proportions for each week & group:
df2 <- df2 %>%
  group_by(group) %>%
  mutate(
    obs_prop = observed / sum(observed)
  )

ggplot(df2, aes(x = week)) +
  geom_col(aes(y = obs_prop), fill = "steelblue", alpha = 0.6) +
  geom_line(aes(y = fitted), color = "firebrick", size = 1) +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    title = "Multinomial fit by group: observed proportions vs fitted",
    y     = "proportion",
    x     = "week (t)"
  ) +
  theme_minimal()
