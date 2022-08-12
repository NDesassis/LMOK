#' New version for the kriging in a linear model (can handle both mono and multi-variable case for CCR)
#' Input: dbin input data base with variables selected with the roles "z" observed value, "f" factors (), and "v" measurement error
#'        dbout output data base to store the results (a selection can be active)
#'        model multivariate stationary model of the covariance between the factors
#'        idx  index of the factor to estimate (by default all the factors are estimated)
#'        flag.estim  computation of the estimation if TRUE
#'        flag.stdev  computation of the standard deviation of the estimation error
#'        radix prefix used to rename the created variables (estimation and standard deviation of the estimation error)
#' Value: the updated dbout data base        
#'
#' 2022-08-10        

#' --------------------------------------
#' Kriging the spatially varying coefficients of a linear model
#' - a single variable ("z")
#' - a measurement error ("v") and 
#' - N factors ("f1", ..., "fN")
#' --------------------------------------

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

#' --------------------------------------
#' Mono-variable and two factors (1, N)
#'  param$type == "IndFac" or "IntFac"
#' --------------------------------------
OK_kriging <- function(dbout, dbin, param, 
                flag.estim = TRUE, flag.stdev = TRUE, 
                radix = "OK", verbose = TRUE){
  
  idx_in_z <- db.getcols(dbin, loctype = "z")
  idx_in_f <- db.getcols(dbin, loctype = "f")

  # Adding Factors and variance of error
  dbin <- db.add(dbin, Un = 1.0, auto.locate = FALSE)
  dbin <- db.add(dbin, sigma2 = param$sigma^2)

  # reseting the locators
  nm_factors <- c("Un", names(dbin@items)[idx_in_f])
  dbin <- db.locerase(dbin, "z")
  dbin <- db.locerase(dbin, "f")
  dbin <- db.locate(dbin, names = nm_factors, loctype = "f", flag.locnew = TRUE)
  dbin <- db.locate(dbin, names = idx_in_z, loctype = "z", flag.locnew = TRUE)
  dbin <- db.locate(dbin, names = "sigma2", loctype = "v", flag.locnew = TRUE)
  
  # model
  if(param$type == "IntFac"){
    cov <- matrix(c(
      param$sigma_1^2, 
      rep(param$tau*param$sigma_1*param$sigma_2,2),
      param$sigma_2^2), 2,2)
    model <- model.create(8,range = param$range, param= param$nu, ndim = 2, nvar = 2, sill = cov)
  }
  if(param$type == "IndFac"){
    cov <- matrix(c(param$sigma_1^2, rep(0.0, 3)), 2,2)
    model <- model.create(8,range = param$range_1, param= param$nu, ndim = 2, nvar = 2, sill = cov)
    cov <- matrix(c(rep(0.0, 3), param$sigma_2^2), 2,2)
    model <- model.create(8,range = param$range_2, param= param$nu, ndim = 2, nvar = 2, sill = cov,
                      model = model)
  }
  
  # Computing the kriging
  dbout <- LMOK_kriging(dbout = dbout, dbin = dbin, model = model, 
               idx = 1:2, radix = radix, flag.estim = flag.estim, flag.stdev = flag.stdev,
               verbose = verbose)
  
  # Renaming according to the naming convenetion
  # renaming the estimated variable
  nm_out_old <- paste(radix, paste0("Z", 1:2), sep = ".")
  nm_out_new <- paste(radix, c("Z0", "Z1"), sep = ".")
  if(flag.estim) {
    dbout = db.rename(dbout, 
                      names =    paste0(nm_out_old, ".estim"), 
                      newnames = paste0(nm_out_new, ".estim") )
  }
  if(flag.stdev) {
    dbout = db.rename(dbout, 
                      names =    paste0(nm_out_old, ".stdev"), 
                      newnames = paste0(nm_out_new, ".stdev") )
  }
  
  
  dbout
}

#' --------------------------------------
#' Bi-variable defined from "code" and 4 factors (1, N)
#'  param$type == "bivar"
#' --------------------------------------
OK_cokriging <- function(dbout, dbin, param, 
                     flag.estim = TRUE, flag.stdev = TRUE,
                     radix = "COK", verbose = TRUE){
  
  idx_in_z <- db.getcols(dbin, "z")
  idx_in_f <- db.getcols(dbin, "f")
  idx_in_c <- db.getcols(dbin, "code")
  stopifnot((idx_in_z > 0)&(idx_in_f > 0)&(idx_in_c > 0))
  stopifnot(param$type == "bivar")

  # Factors
  Type <- dbin[,idx_in_c]
  N    <- dbin[,idx_in_f]
  Un.a <- as.numeric(Type == 1)
  Un.b <- as.numeric(Type == 2)
  N.a  <- N; N.a[Type != 1] <- 0.0
  N.b  <- N; N.b[Type != 2] <- 0.0
  Sigma2 <- rep(0.0, length(N))
  Sigma2[Type == 1] <- param$sigma[1]^2
  Sigma2[Type == 2] <- param$sigma[2]^2
  
  dbin <- db.add(dbin, Un.a, Un.b, N.a, N.b, Sigma2)
  dbin <- db.locerase(dbin, "z")
  dbin <- db.locerase(dbin, "f")
  nm_factors <- c("Un.a", "Un.b", "N.a", "N.b")
  dbin <- db.locate(dbin, names = nm_factors, loctype = "f", flag.locnew = TRUE)
  dbin <- db.locate(dbin, names = idx_in_z, loctype = "z", flag.locnew = TRUE)
  dbin <- db.locate(dbin, names = "Sigma2", loctype = "v", flag.locnew = TRUE)
  
  # models
  cov0 <- matrix(c(param$coef_0[1]^2, rep(prod(param$coef_0),2), param$coef_0[2]^2), 2,2)
  cov0 <- as.matrix(rbind(
    cbind(cov0, matrix(0.0, nrow = 2, ncol = 2)),
    matrix(0.0, nrow = 2, ncol = 4)))
  model <- model.create(8,range = param$range[1], param= param$nu, ndim = 2, nvar = 2, sill = cov0)
  cov1 <- matrix(c(param$coef_1[1]^2, rep(prod(param$coef_1),2), param$coef_1[2]^2), 2,2)
  cov1 <- as.matrix(rbind(
    matrix(0.0, nrow = 2, ncol = 4),
    cbind(matrix(0.0, nrow = 2, ncol = 2), cov1)
  ))
  model <- model.create(8,range = param$range[2], param= param$nu, ndim = 2, nvar = 2, sill = cov1, model = model)
  
  # do the estimation
  dbout <- LMOK_kriging(dbout = dbout, dbin = dbin, model = model, idx = 1:4, 
               flag.estim = flag.estim, flag.stdev = flag.stdev, 
               radix = radix, verbose = TRUE)
  
  # renaming the estimated variable
  nm_out_old <- paste(radix, paste0("Z", 1:4), sep = ".")
  nm_out_new <- paste(radix, c("a.Z0", "b.Z0", "a.Z1", "b.Z1"), sep = ".")
  if(flag.estim) {
    dbout = db.rename(dbout, 
                      names =    paste0(nm_out_old, ".estim"), 
                      newnames = paste0(nm_out_new, ".estim") )
  }
  if(flag.stdev) {
    dbout = db.rename(dbout, 
                      names =    paste0(nm_out_old, ".stdev"), 
                      newnames = paste0(nm_out_new, ".stdev") )
  }
  
 dbout 
}

