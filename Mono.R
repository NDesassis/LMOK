param_ini <- list(nodes   = nodes, margin  = margin,
                  type = "MonoF", 
                  nu = nu,
                  sigma_1 = 10,
                  sigma   = 10,
                  range_1 = 10000
)

#' ------------------------------------------------------------
#'  Auxiliary function: 1 factor in mono-variable
#' x[1] = sigma_1 (Standard deviation of Z_1)
#' x[2] = range_1 (Practical range of Z_1)
#' x[3] = sigma   (Standard deviation of the measurement error)
#' ------------------------------------------------------------
modelMulti<- function(x, verbose = TRUE){
  nu = param$nu
  model <- model.create(8, range = abs(x[2]), sill = x[1], param = param$nu)
  if (verbose){
    print(model)
    print(paste0("sigma = ", x[3]))
  }
  list(model = model, sigma = x[3])
}

x_ini= c(param_ini$sigma_1, param_ini$range_1, param_ini$sigma)

cont = list(fnscale = 1,parscale = c(1,100000,1), ndeps=rep(1,3)*1e-4,trace = 1, maxit=10,REPORT=1)


namesParam = c("sigma_1", "range_1", "sigma")
