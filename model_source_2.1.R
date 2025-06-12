u <- rep(0,max(d$t_k))
u_i <- aggregate(list(u = d$n),by = list(s = d$t_k),sum)
data$u <- u[u_i$s] <- u_i$u

sp_count <- d %>%
  group_by(t_l) %>%
  summarise(n = sum(n))
data$sp_pr <- as.vector(sp_count$n/sum(sp_count$n))


# make_param_list('~1', d, )
params <- make_all_params(formulas
                          , d
                          , s = max(na.omit(c(d$t_k,d$r_k)))
                          , state_col = data$state)


# make_param_list('~1', d, )
struct <- extract_all_structures(formulas
                                 , d
                                 , s = max(na.omit(c(d$t_k,d$r_k)))
                                 , state = data$state)

data$struct <- struct
data$sp_count <- sp_count

Nsuper_g <- data$struct$groups$Nsuper[[1]]
if(length(Nsuper_g)>1){
  data$uTot_g <- as.numeric( with(data, tapply(n, Nsuper_g, sum)) )
}else{
  data$uTot_g <- sum(data$uTot)
}


n_state <- if(!is.null(data$state)) length(unique(data[[data$state]])) else 1L

map <- list()
random <- c( "logit_v")
if(length(data$struct$Z_list$Nsuper)>0){
  random <- c(random, "u_Nsuper")
}
if(length(data$struct$Z_list$phi)>0){
  random <- c(random, "u_phi")
}
if(length(data$struct$Z_list$p)>0){
  random <- c(random, "u_p")
}
if(length(data$struct$Z_list$w)>0){
  random <- c(random, "u_w")
  map$log_sd_w <- as.factor(NA)
}
if(length(params$t_var)>1){
  random <- c(random, "t_var")
}


environment(model) <- .GlobalEnv
map <- list(M_sigma = as.factor(NA)
            # logit_v = as.factor(rep(NA,length(params$logit_v)))
)

obj <- RTMB::MakeADFun(func = model,
                       parameters = params,
                       data = data,
                       random = random,
                       map = map,
                       silent = FALSE)

# Optimize Model
opt <- nlminb(obj$par, obj$fn, obj$gr)
# sdr <- sdreport(obj)
rep <- obj$report()


