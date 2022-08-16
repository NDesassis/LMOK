param_ini <- list(nodes   = nodes, margin  = margin,
                  type = "IndFac", 
                  nu = nu,
                  sigma_0a = 10,
                  sigma_0b = 10,
                  sigma_1a = 10,
                  sigma_1b = 10,
                  tau_0  = 0.01,
                  tau_1  = 0.01,
                  range_1 = 100000,
                  range_2 = 100000,
                  sigma_a = 20,
                  sigma_b = 20
)
#' ------------------------------------------------------------
#' Auxiliary function: bivar (intrinsic coregionalisation of the factors)
#' Factors: F1 = "indicator of type == a"
#'          F2 = "indicator of type == b"
#'          F3 = "N x indicator of type == a"
#'          F4 = "N x indicator of type == b"
#'          
#' x[1] = sigma_0a Standard deviation of Z_0) F1
#' x[2] = sigma_0b Standard deviation of Z_1) F2
#' x[3] = sigma_1a Standard deviation of Z_2) F3
#' x[4] = sigma_1b Standard deviation of Z_3) F4
#' x[5] = tau_0   correlation coefficient between F1 and F2
#' x[6] = tau_1   correlation coefficient between F3 and F4
#' x[7] = range_1 range of spatial structure of F1 and F2
#' x[8] = range_2 range of spatial structure of F3 and F4
#' x[9] = sigma_a Standard deviation of the measurement error for type == a
#' x[10] = sigma_b Standard deviation of the measurement error for type == b
#' ------------------------------------------------------------

modelMulti <- function(x, verbose = TRUE){
  nu = 1
  # F1 and F3
  c1 <- as.matrix(rbind(
    cbind(matrix(c(x[1]^2, rep(x[1]*x[2]*tanh(x[5]), 2), x[2]^2), 2, 2),
          matrix(0.0,2,2)),
    matrix(0.0, 2, 4)))
  
  model <- model.create(8, range = abs(x[7]), sill = c1)
  
  c2 <- as.matrix(rbind(
    matrix(0.0,2,4),
    cbind(matrix(0.0, 2,2), 
          matrix(c(x[3]^2, rep(x[3]*x[4]*tanh(x[6]), 2), x[4]^2), 2, 2))))
  model <- model.create(8, range = abs(x[8]), sill = c2, model = model)
  
  if (verbose){
    print(model)
    print(paste0("sigma_a = ", x[9]))
    print(paste0("sigma_b = ", x[10]))
  }
  list(model = model, sigma = c(x[9], x[10]))
}


x_ini=c(param_ini$sigma_0a,
        param_ini$sigma_0b,
        param_ini$sigma_1a, 
        param_ini$sigma_1b,
        param_ini$tau_0   ,
        param_ini$tau_1   ,
        param_ini$range_1 ,
        param_ini$range_2 ,
        param_ini$sigma_a ,
        param_ini$sigma_b
)

cont = list(fnscale = 1,parscale = c(1,1,1,1,1,1,10000,10000,1,1), ndeps=rep(1,10)*1e-4,trace = 1, maxit=200,REPORT=1)

namesParam =  c("sigma_0a",
                "sigma_0b",
                "sigma_1a", 
                "sigma_1b",
                "tau_0"   ,
                "tau_1"   ,  
                "range_1" ,
                "range_2" , 
                "sigma_a" ,
                "sigma_b"
)
