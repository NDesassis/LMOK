if(type == "flat")
  param_ini <- list(nodes   = nodes, margin  = margin,
                    type = "IntFac", 
                    nu = nu,
                    sigma_0 = 18.878809,
                    sigma_1 = 5.477902,
                    tau     = -1.519629,
                    sigma   = 11.587706,
                    range   = 302.334218
  )

if(type =="house")
param_ini <- list(nodes   = nodes, margin  = margin,
                  type = "IntFac", 
                  nu = nu,
                  sigma_0 = 58.61,
                  sigma_1 = 119.83,
                  tau     = 0.2,
                  sigma   = 10,
                  range   = 639805.5
)

#' ------------------------------------------------------------
#'  Auxiliary function: 2 factors in intrinsic correlation
#' x[1] = sigma_0 (Standard deviation of Z_0)
#' x[2] = range_1 (Standard deviation of Z_1)
#' x[3] = atanh(rho) (Coefficient of correlation between the two factors)
#' x[4] = range   (Practical range of Z_0 and Z_1)
#' x[5] = sigma   (Standard deviation of the measurement error)
#' ------------------------------------------------------------
modelMulti<- function(x, verbose = TRUE){
  nu = param$nu
  c0 <- matrix(c(x[1]^2, rep(x[1]*x[2]*tanh(x[3]),2), x[2]^2),nrow=2)
  model <- model.create(8, range = abs(x[4]), sill = c0)
  if (verbose){
    print(model)
    print(paste0("sigma = ", x[5]))
  }
  list(model = model, sigma = x[5])
}

x_ini= c(param_ini$sigma_0, param_ini$sigma_1, atanh(param_ini$tau), param_ini$range, param_ini$sigma)

cont = list(fnscale = 1,parscale = c(1,1,1,100000,1), ndeps=rep(1,5)*1e-4,trace = 1, maxit=300,REPORT=1)


namesParam = c("sigma_0", "sigma_1", "tau", "range", "sigma")
