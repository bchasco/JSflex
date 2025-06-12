
Nv_e <- as.list(sdr, "Estimate",report=TRUE)$Nsuper_vec
Nv_sd <- as.list(sdr, "Std. Error",report=TRUE)$Nsuper_vec
N_e <- as.list(sdr, "Estimate",report=TRUE)$Ntotal
N_sd <- as.list(sdr, "Std. Error",report=TRUE)$Ntotal


g1 <- data.frame(N = Nv_e, sd = Nv_sd, level = levels(data$group_factor$Nsuper)) %>%
  ggplot(aes(x = level, y = N)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = N - 1.96 * sd, ymax = N + 1.96 * sd)) +
  theme_bw() +
  theme(text= element_text(size = 16))


g2 <- data.frame(N = N_e, sd = N_sd, level = "Total") %>%
  ggplot(aes(x = level, y = N)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = N - 1.96 * sd, ymax = N + 1.96 * sd)) +
  theme_bw() +
  theme(text= element_text(size = 16))


ggpubr::ggarrange(g1,g2,ncol=2)
