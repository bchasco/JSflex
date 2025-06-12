# -----------------------------------------------------------
# 1.  helper -- collect all grouping vars except the period
# -----------------------------------------------------------
get_joint_group_vars <- function(formulas, period = "f_tk") {
  vars <- character()

  for (fml in formulas) {
    ## fixed-effect RHS
    rhs_vars  <- if (length(fml) == 2) all.vars(fml[[2]])
    else                    all.vars(fml[[3]])

    ## random-effect bars
    bar_vars  <- unlist(lapply(lme4::findbars(fml),
                               function(b) all.vars(b[[3]])))

    vars <- union(vars, c(rhs_vars, bar_vars))
  }

  setdiff(vars, period)          # drop the period factor
}

# -----------------------------------------------------------
# 2.  build the joint grouping factor *once*
#     (call this right after you read the data)
# -----------------------------------------------------------
build_joint_group <- function(data, formulas, period = "f_tk") {

  joint_vars <- get_joint_group_vars(formulas, period)

  if (length(joint_vars) == 0L) {
    joint_fac <- factor(rep(1L, nrow(data)))
  } else {
    joint_fac <- interaction(data[ joint_vars ], drop = TRUE)
  }

  list(joint      = joint_fac,
       joint_vars = joint_vars)   # handy to keep the names!
}


jg_out <- build_joint_group(d, formulas, period = "f_tk")
joint_group <- jg_out$joint          # factor with 15 levels
vars_used   <- jg_out$joint_vars     # c("f_sex", "f_tl")

levels(joint_group)
## [1] "Female.1" "Male.1" "Jack.1" "Female.2" … "Jack.5"

# -----------------------------------------------------------
# 4.  tidy label data-frame for plotting --------------------
# -----------------------------------------------------------
library(tidyr)
library(dplyr)

label_df <- tibble(level = levels(joint_group)) %>%
  separate(level,
           into = vars_used,          # c("f_sex","f_tl")
           sep  = "\\.", convert = TRUE) %>%
  mutate(row = row_number())          # 1 … 15  (row index in pent_mat etc.)

label_df
## # A tibble: 15 × 3
##    f_sex  f_tl   row
##    <chr> <int> <int>
##  1 Female     1     1
##  2 Male       1     2
##  3 Jack       1     3
##  … …

# -----------------------------------------------------------
# 5.  give names to the matrices and plot -------------------
# -----------------------------------------------------------
pent_mat <- rep$pent_mat
colnames(pent_mat) <- 1:nlevels(data$period_fac)
# S2 <- t(as.matrix(rep$S2[,1,1,]))
S2 <- rep$S2[,1,1,]
colnames(S2) <- 1:nlevels(data$period_fac)
# D2 <- t(as.matrix(rep$D2[,1,1,]))
D2 <- rep$D2[,1,1,]
colnames(D2) <- 1:nlevels(data$period_fac)

rownames(pent_mat) <- label_df$level                # optional
rownames(S2)       <- label_df$level
rownames(D2)       <- label_df$level

## long data for pent_mat
pent_long <- as.data.frame(pent_mat) |>
  mutate(row = row_number()) |>
  mutate(proc = "Arrival/pent_mat") |>
  pivot_longer(-c(row,proc),
               names_to  = "period",
               values_to = "val") |>
  left_join(label_df, by = "row")

S2_long <- as.data.frame(S2) |>
  mutate(row = row_number()) |>
  mutate(proc = "S/phi") |>
  pivot_longer(-c(row,proc),
               names_to  = "period",
               values_to = "val") |>
  left_join(label_df, by = "row")

D2_long <- as.data.frame(D2) |>
  mutate(row = row_number()) |>
  mutate(proc = "D/p") |>
  pivot_longer(-c(row,proc),
               names_to  = "period",
               values_to = "val") |>
  left_join(label_df, by = "row")

df <- rbind(pent_long,S2_long,D2_long)

library(ggplot2)
ggplot(df ,
       aes(x = as.numeric(period),
           y = val,
           colour  = f_yr,
           # linetype = as.factor(f_tl),
           group = row)) +
  geom_line(size = 1.4) +
  # ylim(0,1) +
  labs(x = "Period",
       y = "Probability"
       ,colour   = "Year"
       # ,linetype = "Location (f_tl)"
       ) +
  theme_classic() +
  facet_wrap(~proc, scales = "free_y", ncol = 1) +
  theme(text = element_text(size = 20))
