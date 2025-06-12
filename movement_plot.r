# Assume rep$out$M is a 3D array: [G, n, n] — G groups of n × n matrices
library(Matrix)

# Number of matrix slices (G), rows and cols (n)
G <- dim(rep$out$M)[1]
n <- dim(rep$out$M)[2]

# Initialize output array
if(G>1){
  M_exp <- array(NA_real_,
                 dim = c(G, n, n),
                 dimnames = list(sex = levels(d$f_sex), from = sort(unique(d$t_l)), to = sort(unique(d$t_l))))
}else{
  M_exp <- array(NA_real_, dim = c(G, n, n), dimnames = list(sex = "Combined", from = sort(unique(d$t_l)),
                                                             to = sort(unique(d$t_l))))
}

for (g in seq_len(G)) {
  mat <- rep$out$M[g, , ]
  M_exp[g, , ] <- as.matrix(expm(mat))
}

M_exp <- as.data.frame(reshape2::melt(M_exp))
M_exp <- M_exp[M_exp$value>0,]
M_exp <- M_exp[M_exp$from<5,]

ggplot(M_exp, aes(y = from, x = to, fill = value)) +
  facet_wrap(~sex, ncol = 2) +
  geom_tile() +
  scale_y_reverse() +
  scale_fill_gradient2(
    low = "red",
    mid = "yellow",
    high = "green",
    midpoint = 0.5
  ) +
  theme_classic() +
  theme(text = element_text(size = 16))
