design <- make_design_list_2(formulas, list(state = input$state, input = input$time), d)
lk <- make_process_lookup(names(design$formulas), design, d, state_var = design$state, period_var = "f_tk")
params    <- make_param_vectors(design, names(design$formulas))
random <- make_random(params$u)

flatten_params <- function(params) {
  out <- list()
  for (group in names(params)) {
    sublist <- params[[group]]
    for (nm in names(sublist)) {
      mat <- as.vector(sublist[[nm]])
      colnames(mat) <- NULL
      out_name <- paste(group, nm, sep = "_")
      out[[out_name]] <- mat
    }
  }
  out
}

# example
params <- flatten_params(params)
names(params)

data <- make_RTMB_data_list(design, d, state = design$state, period = "f_tk")
data$lk <- lk


source("model_2.3.r")
obj <- RTMB::MakeADFun(func = model,
                       data = data,
                       random = random,
                       map = list(
                         beta_t_var = as.factor(NA)),
                       parameters = params)
# # # Optimize Model
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep_ <- obj$report()
# sdr <- sdreport(obj)
# sd_val <- as.list(sdr, "Estimate", report = TRUE)
# sd_sd <- as.list(sdr, "Std. Error", report = TRUE)

out <- list(opt = opt,
            obj = obj,
            data = data,
            params = params,
            rep = rep_
            # ,sd = list(est = sd_val,
            #           sd = sd_sd)
            )

