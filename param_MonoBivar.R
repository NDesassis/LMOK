param_ini <- list(nodes   = nodes, margin  = margin,
                  type = "MonoBivar", 
                  nu = nu,
                  sigma_1a = 10,
                  sigma_1b = 15,
                  tau_1  = 0.01,
                  range_1 = 100000,
                  sigma_a = 20,
                  sigma_b = 15
)
#' ------------------------------------------------------------
#' Auxiliary function: bivar (intrinsic coregionalisation of the factors)
#' Factors: 
#'          F1 = "N x indicator of type == a"
#'          F2 = "N x indicator of type == b"
#'          
#' x[1] = sigma_1a (Standard deviation of F1)
#' x[2] = sigma_1b (Standard deviation of F2) 
#' x[3] = tau_1    correlation coefficient between F1 and F2
#' x[4] = range_1 range of spatial structure of F1 and F2
#' x[5] = sigma_a Standard deviation of the measurement error for type == a
#' x[6] = sigma_b Standard deviation of the measurement error for type == b
#' ------------------------------------------------------------

modelMulti <- function(x, verbose = TRUE){
  nu = 1
  # F1 and F2
  c1 <- matrix(c(x[1]^2, rep(x[1]*x[2]*tanh(x[3]), 2), x[2]^2), 2, 2)
  model <- model.create(8, range = abs(x[4]), sill = c1)
  if (verbose){
    print(model)
    print(paste0("sigma_a = ", x[5]))
    print(paste0("sigma_b = ", x[6]))
  }
  list(model = model, sigma = c(x[5], x[6]))
}


x_ini=c(param_ini$sigma_1a,
        param_ini$sigma_1b,
        param_ini$tau_1   ,
        param_ini$range_1 ,
        param_ini$sigma_a ,
        param_ini$sigma_b
)

cont = list(fnscale = 1,parscale = c(1,1,1,10000,1,1), ndeps=rep(1,6)*1e-4,trace = 1, maxit=300,REPORT=1)

namesParam =  c("sigma_1a", 
                "sigma_1b",
                "tau_1"   ,  
                "range_1" ,
                "sigma_a" ,
                "sigma_b"
)
