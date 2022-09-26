# functions to implement LMOK with the ordinary kriging approach

#' compute the estimation of the locally varying parameters
#' 
#' The function computes the kriging of using a unique neighborhood. 
#' Hence, the number of the input data points should be small enough
#' to make the computations tractable (e.g. dbin$nactive <= 100).
#' 
#'@param dbout The db-class structure containing the target file
#'The estimation is performed only for the selected samples of the target file.
#'@param dbin The db-class structure containing the data file
#'The following variables should be defined in the data file:
#' - the observed values (locator "z")
#' - the factors associated to the latent fields to be estimated (locator "f")
#' - the standard deviation of the modeling error (locator "v")
#' The number of latent fields/factors should be equal to the number of variables 
#' of the model (length(db.getcols(dbin, "f")) == model$nvar).
#'@param model The model-class defining the covariance model of the latent fields
#'Only stationary covariance functions should be used as basic structures
#'of the linear model of coregionalization. 
#'@param idx The integer vector containing the indices of the latent fields to be estimated.
#'If idx == NA, the estimation is computed for all the latent fields.
#'@param flag.estim A Boolean value controlling the computation of the kriging value.
#'@param flag.stdev A Boolean value controlling the computation of the kriging standard deviation.
#'@param radix A string used as prefix of the variables created in the output file *dbout*.
#'@param verbose A Boolean variable to control the printed messages
#'@return The target Db where the following variables have been added:
#' - the estimation of the latent fields
#' - the standard deviation of the kriging error.
#' These variables are multiplied for each one of the latent fields (numbered from 1).
LMOK_kriging <- function(dbout, dbin, model, 
                         idx = NA, flag.estim = TRUE, flag.stdev = FALSE, 
                         radix = "LMOK", verbose = FALSE){
  
  idx_in_z <- db.getcols(dbin, "z")
  idx_in_f <- db.getcols(dbin, "f")
  idx_in_v <- db.getcols(dbin, "v")
  
  L <- model$nvar
  stopifnot((idx_in_z > 0)&(idx_in_v > 0))
  stopifnot(L == length(idx_in_f))  
  stopifnot(flag.estim|flag.stdev)
  
  if(min(is.na(idx))){idx = 1:L}
  stopifnot(is.element(idx, 1:L))
  
  if(verbose){ 
    print(paste0("LMOK::Kriging: Input data"))
    print(paste0(" - number of data:    ", dbin$nactive, "/", dbin$nech))
    print(paste0(" - number of targets: ", dbout$nactive, "/", dbout$nech))    
    print(paste0(" - number of factors: ", L))
  }
  
  if(verbose){ print(paste0("LMOK::Kriging: computing LHS"))}
  
  Factors <- as.matrix(db.extract(dbin,  names = idx_in_f, flag.compress = TRUE),
                       nrow = dbin$nactive, ncol = length(idx_in_f))
  Values  <- db.extract(dbin,  names = idx_in_z, flag.compress = TRUE)
  Sigma2  <- db.extract(dbin,  names = idx_in_v, flag.compress = TRUE)
  C <- diag(Sigma2)
  for (i in 1:L) {
    for (j in 1:L) {
      C <- C + model.covmat(db1 = dbin, db2 = dbin, model = model,
                            ivar = i, jvar = j, as.cov = TRUE)*
        outer(X = Factors[,i], Y = Factors[,j], FUN = "*")
    }
  }
  LHS <- as.matrix(rbind(
    as.matrix(cbind(C, Factors)),
    as.matrix(cbind(t(Factors), matrix(0,L,L)))
  ))
  invLHS <- solve(LHS)
  
  for (i0 in idx) {
    
    if(verbose){ print(paste0("LMOK::Kriging: factor ", i0))}
    # Computing RHS for i0
    c0 <- matrix(0.0, nrow = dbin$nactive, ncol = dbout$nactive)
    for (j in 1:L){
      c0 <- c0 + 
        model.covmat(
          db1 = dbin, db2 = dbout, model = model,
          ivar = j, jvar = i0, as.cov = TRUE)*
        outer(X = Factors[,j], Y = rep(1,dbout$nactive), FUN = "*")
    }
    m_i0 <- matrix(0.0, nrow = L, ncol = dbout$nactive)
    m_i0[i0,] <- 1.0
    RHS <- as.matrix(rbind(c0, m_i0))
    
    lambdaNu <- invLHS %*% RHS
    
    if(flag.estim){
      nm_estim <- paste(radix, paste0("Z", i0), "estim", sep = ".")
      Z_estim  <- as.numeric(t(c(Values, rep(0.0, L))) %*% lambdaNu)
      dbout <- db.replace(dbout, name = nm_estim, tab = Z_estim, flag.sel = TRUE)
    }
    if(flag.stdev){
      nm_stdev <- paste(radix, paste0("Z", i0), "stdev", sep = ".")
      Z_stdev  <- rep(model$silltot[i0, i0], dbout$nactive) -
        vapply(
          X = 1:dbout$nactive, 
          FUN = function(i){
            sum(lambdaNu[,i]*RHS[,i])}, 
          FUN.VALUE = numeric(1))
      Z_stdev <- sqrt(Z_stdev)
      dbout <- db.replace(dbout, name = nm_stdev, tab = Z_stdev, flag.sel = TRUE)
    }
  } # Loop over the factors to be estimated
  dbout
}

