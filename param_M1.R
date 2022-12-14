if(type == "flat")
  param_ini <- list(nodes   = nodes, margin  = margin,
                    type = "IndFac", 
                    nu = nu,
                    sigma_0 = 6.888328,
                    range_0 = 184.315752,
                    sigma_1 = 8.770781,
                    range_1 = 47581.246797,
                    sigma   = 12.273570
                  
 )


if(type == "house")
  param_ini <- list(nodes   = nodes, margin  = margin,
                    type = "IndFac", 
                    nu = nu,
                    sigma_0 = 10,
                    range_0 = 10000,
                    sigma_1 = 10,
                    range_1 = 10000,
                    sigma   = 31
                    
  )

#' ------------------------------------------------------------
#'  Auxiliary function: BiInt (independence of the 4 factors)
#' x[1] = sigma_1 (Standard deviation of Z_1)
#' x[2] = range_1 (Practical range of Z_1)
#' x[3] = sigma_2 (Standard deviation of Z_2)
#' x[4] = range_1 (Practical range of Z_2)
#' x[5] = sigma   (Standard deviation of the measurement error)
#' ------------------------------------------------------------
modelMulti<- function(x, verbose = TRUE){
  nu = param$nu
  c0 <- matrix(c(x[1]^2, rep(0,3)),nrow=2)
  model <- model.create(8, range = abs(x[2]), sill = c0)
  c1 <- matrix(c(rep(0,3),x[3]^2),nrow=2)
  model <- model.create(8, range = abs(x[4]), sill = c1,model=model)
  if (verbose){
    print(model)
    print(paste0("sigma = ", x[5]))
  }
  list(model = model, sigma = x[5])
}

x_ini= c(param_ini$sigma_0,param_ini$range_0,param_ini$sigma_1,param_ini$range_1,param_ini$sigma)

cont = list(fnscale = 1,parscale = c(1,100000,1,100000,1), ndeps=rep(1,5)*1e-4,trace = 1, maxit=300,REPORT=1)

x_ini=c(param_ini$sigma_0,
        param_ini$range_0,
        param_ini$sigma_1,
        param_ini$range_1,
        param_ini$sigma
)



namesParam = c("sigma_0", "range_0", "sigma_1", "range_1", "sigma")
