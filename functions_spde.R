# --------------------------------------------------------------------
# Kriging functions
# --------------------------------------------------------------------
# New version of the functions (combining the two mono-variate model)

# 
# param: a list of parameters used for the kriging of the drift parameters

#  nu:  is the regularity coefficient of the matern covariance function
#  sigma: is the nugget effect of the modeling error

#  type  \in {"IndFact", "IntFact"} defines the type of model
#    "IndFact" the two FA are independent and defined by their respective sill and range
#    In this case, the following parameters should be defined
#     range_1: range of the Matern covariance function of the first (constant) parameter
#     sigma_1: standard deviation of the first parameter
#     range_2: range of the Matern covariance function of the second parameter
#     sigma_2: standard deviation of the second parameter

#    "IntFact" The two FA are in intrinsic correlation defined by their respective sill, 
#     the coefficient of correlation, and the range of the common spatial structure
#    In this case, the following parameters should be defined
#     range:   range of the common Matern covariance function
#     sigma_1: standard deviation of the first parameter
#     range_2: range of the Matern covariance function of the second parameter
#     tau    : coefficient of corstandard deviation of the second parameter

# If the mesh is not provided, the mesh is computed using the bounding box of the
# input data set and the following parameters

#  border: proportion of the range to extend the bounding box (if mesh not provided)
#  nodes:  number of vertices along the main direction (if mesh not provide)


# Performing the kriging or the conditional simulation according to the value of nsim.
# If nsim == 0, only the kriging is computed and stored
# If nsim >  0, only the simulations are performed and according to the value of the flag.ce 
# the simulations are stored or the conditinal expectation (mean and std) are computed and stored.

#' ----------------------------------------
#' Function to compute LMOK using SPDE
#' ----------------------------------------
#' Multi-facteur
#'  function to initialize

SPDE_LMOK_init_geometry <- function(dbin, mesh = NA, verbose = TRUE){
  
  idx_in_z <- db.getcols(dbin, "z")
  idx_in_f <- db.getcols(dbin, "f")
  stopifnot((idx_in_z > 0))

  # Geometry (if the mesh is not available, it is computed from the limits of dbin)
  if (is.na(mesh)){
    bb    <- dbin$limits
    bb    <-  bb + outer(c(-1,1), (bb[2,] - bb[1,])* param$margin, "*") 
    dx    <- (bb[2,] - bb[1,]) / (param$nodes - 1)
    mesh  <- meshing(extendmin = bb[1,], extendmax = bb[2,], cellsize = dx)
  }
  n.mesh <- mesh$npoints
  
  # Input data (active selection & factors & observation & measurement error)
  X  <- as.matrix(db.extract(dbin, names = idx_in_f, flag.compress = TRUE))
  Y  <- db.extract(dbin, names = idx_in_z, flag.compress = TRUE)

  # Computing the projection matrix for the data
  AprojMean  <- mesh.barycenter(dbin, mesh)
  for (i in 1:length(idx_in_f)){
    if(i == 1) {
      Aproj <- Diagonal(x=X[,i])%*%AprojMean
    } else {
      Aproj <- cbind(Aproj, Diagonal(x=X[,i])%*%AprojMean)
    }
  }
  if(verbose){ 
    print(paste0("SPDE_LMOK::geometry: Input data"))
    print(paste0(" - number of data:    ", dbin$nactive, "/", dbin$nech))
    print(paste0(" - size of the mesh : ", n.mesh))
    print(paste0(" - number of factors: ", length(idx_in_f)))
    print(paste0(" - return values    : "))
    print(paste0(" - dim(X)    : ", dim(X)))
    print(paste0(" - dim(Y)    : ", length(Y)))
    print(paste0(" - dim(Aproj): ", dim(Aproj)))
  }
  

  return(list(
    X     = X,
    Y     = Y,
    mesh  = mesh,
    Aproj = Aproj
  ))
}

SPDE_LMOK_init_model <- function(model, sigma, geo, verbose = TRUE){
  
  n.var  <- model$nvar
  n.cova <- model$ncova
  stopifnot(n.var == ncol(geo$X))  

  if(length(sigma) == 1) {sigma <- rep(sigma, length(geo$Y))}
  if(verbose){
    print(paste0("SPDE_LMOK_init_model: length(sigma) =", length(sigma)))
  }
  stopifnot(length(sigma) == length(geo$Y))
  
  # defining the link between the structures and the factors
  nnF <- rep(0, n.var)
  idF <- rep(0, n.var)
  for (i in 1:n.cova){
    idx <- (1:n.var)[diag(model.reduce(model, structs = i)$silltot) > 0]
    nnF[idx] <- nnF[idx] + 1
    idF[idx] <- i
  }
  stopifnot(nnF == 1)
  
  if(verbose){ 
    print(paste0("SPDE_LMOK::init: Model"))
    print(paste0(" - structures nbr.  : ", n.cova))    
  }

  # Computing the precision matrix on the mesh
  expand.idx <- function(idx){
    I = as.numeric(outer(X=idx, Y=idx, FUN = function(i,j){i}))
    J = as.numeric(outer(X=idx, Y=idx, FUN = function(i,j){j}))
    list(I = I, J = J)
  }
  
  # Loop over the structures of the model
  D    <- Diagonal(x= 1/sigma^2)
  prec <- t(geo$Aproj)%*%D%*%geo$Aproj
  items<- list()
  for (i in 1:n.cova){
    ci  <- model@basics[[i]]
    # Matrix Mi
    cMi  <- ci@sill
    idx  <- which(diag(cMi) > 0)
    qMi  <- solve(cMi[idx, idx])
    i.Mi <- expand.idx(idx)  
    x.Mi <- as.numeric(qMi)
    Mi   <- sparseMatrix(i = i.Mi$I, j = i.Mi$J, x = x.Mi, dims = c(n.var, n.var))    
    # Matrix Qi
    mi  <- model.create(
      vartype = ci@vartype, 
      range   = ci@range,
      param   = ci@param
    )
    Qi   <- spde.matrices(model=mi,flag.Q = TRUE, mesh=geo$mesh)$Q
    
    # updating precision matrix
    prec <- prec + kronecker(Mi, Qi)
    
    if(verbose) {
      print(paste0(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"))
      print(paste0("SPDE_LMOK: structure #", i))
      print(Mi)
      print(mi)
      print(paste0(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"))
    }
    items[[1+length(items)]] <- list(M = qMi, Q = Qi)
  }
  return(list(
    items    = items,
    D        = D,
    cholPrec = Cholesky(prec, LDL=FALSE)
  ))
}
  
SPDE_LMOK_init_spde <- function(geo, mod, verbose = TRUE){
  computeInvSigma <- function(X, A, cholPrec, D){
    D %*% (X - A %*% (solve(cholPrec, t(A) %*% D %*% X)))
  }
  # computing the estimate on the mesh
  invSigmaX <- computeInvSigma(geo$X, geo$Aproj, mod$cholPrec, mod$D)
  coeff     <- solve(t(geo$X)%*%invSigmaX, t(invSigmaX))%*%geo$Y
  Sc        <- geo$Y - geo$X %*% coeff
  Z         <- solve(mod$cholPrec, t(geo$Aproj)%*%mod$D%*%Sc)
  invSigmaSc<- computeInvSigma(Sc, geo$Aproj, mod$cholPrec, mod$D)
  if (verbose) {
    for (i in seq_along(coeff)){
      print(paste0("SPDE_LMOK: coeff[",i,"] = ", coeff[i]))
    }
    print(paste0("SPDE_LMOK: length(Z)       = ", length(Z)))
    print(paste0("SPDE_LMOK: length(Sc)      = ", length(Sc)))
    print(paste0("SPDE_LMOK: dim(invSigmaX)  = ", dim(invSigmaX)))
    print(paste0("SPDE_LMOK: dim(invSigmaSc) = ", dim(invSigmaSc)))
  }
  return(list(
    invSigmaX  = invSigmaX,
    invSigmaSc = invSigmaSc,
    Sc         = Sc,
    Z          = Z,
    coeff      = coeff
  ))
}


SPDE_LMOK_krigsim <- function(dbout, dbin, model, mesh = NA,
                              nsim = 0, seed = NA, flag.ce = FALSE, 
                              radix = "SPDE_LMOK", verbose = TRUE){
  
  # measurement error
  sigma <- sqrt(db.extract(dbin, db.getname(dbin, "v", 1), flag.compress = TRUE))

  geo   <- SPDE_LMOK_init_geometry(dbin, mesh, verbose)
  mod   <- SPDE_LMOK_init_model(model, sigma, geo, verbose)
  spde  <- SPDE_LMOK_init_spde(geo, mod, verbose)
  
  # computing simulations on the mesh (if necessary)
  if (nsim > 0){
    if(!is.na(seed)) {set.seed(seed = seed)}
    U      <- matrix(rnorm(length(spde$Z)*nsim), nrow = length(spde$Z), ncol = nsim)
    Res    <- t(expand(mod$cholPrec)$P)%*%solve(mod$cholPrec, U,system="Lt")
  }
  
  # Compute on the output data base dbout
  A_out   <- mesh.barycenter(dbout, geo$mesh)
  n.mesh  <- geo$mesh$npoints
  if(verbose){
    print(paste0("krigsim:"))
    print(paste0(" - length(spde$Z) = ", length(spde$Z)))
    print(paste0(" - length(spde$coeff) = ", length(spde$coeff)))
    print(paste0(" - n.mesh         = ", n.mesh))
    print(paste0(" - dim(A_out)     = ", dim(A_out)))
  }
  for (i in 1:length(spde$coeff)) {
    if (verbose) {print(paste0("Processing factor #",i))}
    # compute the estimate
    Z_out <- (A_out%*%spde$Z[1:n.mesh+(i-1)*n.mesh,1])[,1] + spde$coeff[i]
    # storing the estimate
    if (nsim == 0){
      nm_estim <- paste(radix, paste0("Z",i), "estim", sep = ".")
      dbout    <- db.replace(dbout, name = nm_estim, tab  = Z_out)
    } 
    # processing the simulations
    if (nsim > 0) {
      Z_sim <- as.matrix(Z_out + (A_out%*%Res[1:n.mesh+(i-1)*n.mesh,]))
      
      # storing the conditional expectation  
      if (flag.ce){
        nm_ce_estim <- paste(radix, paste0("Z",i), "ce.estim", sep = ".")
        dbout <- db.replace(dbout, name = nm_ce_estim, 
                            tab = apply(X = Z_sim, MARGIN = 1, FUN = "mean"))
        nm_ce_stdev <- paste(radix, paste0("Z",i), "ce.stdev", sep = ".")
        dbout <- db.replace(dbout, name = nm_ce_stdev, 
                            tab = apply(X = Z_sim, MARGIN = 1, FUN = "sd"))
      }
      # storing the simulations
      if (!flag.ce) {
        for (n in 1:nsim){
          nm_sim <- paste(radix, paste0("Z",i,".S"), n, sep = ".")
          dbout <- db.replace(dbout, name = nm_sim, tab = Z_sim[,i])
        }
      }
    } # end of simulations
  } # loop on the factors
  dbout 
}


SPDE_LMOK_init_MLL <- function(dbin, param2model, mesh = NA, verbose = TRUE){
  
  geo <- SPDE_LMOK_init_geometry(dbin = dbin, mesh = mesh, verbose = verbose)
  
  fn_LL <- function(x, verbose = FALSE){
    # convert the parameters into a RGeostats model
    res   <- param2model(x, verbose = verbose)
    sigma <- res$sigma
    # Multivariate case
    if (length(res$sigma) > 1){
      if (verbose){print(paste0("Multivariate case: nvar = ", length(sigma)))}
      idx_in_c <- db.getcols(dbin, loctype = "code", rank.match = 1)
      stopifnot(!is.null(idx_in_c))
      code <- db.extract(dbin, names = idx_in_c, flag.compress = TRUE)
      l_code <- sort(unique(code))
      stopifnot(is.element(l_code, 1:length(sigma)))
      sigma <- res$sigma[code]
    }
    mod  <- SPDE_LMOK_init_model(res$model, sigma, geo, verbose)
    spde <- SPDE_LMOK_init_spde(geo, mod, verbose = verbose)
    
    # compute the likelihood
    ll_val <- 0.0
    # compute the quadratic form
    ll_val <- -1/2 * sum(spde$Sc*spde$invSigmaSc) 
    # compute the determinant
    detSigma = 2 * determinant(mod$cholPrec,logarithm = TRUE)$modulus[1]
    for (i in seq_along(mod$items)){
      detQi <- 2 * determinant(Cholesky(mod$items[[i]]$Q,LDL=F),logarithm = T)$modulus[1]
      detMi <- log(det(mod$items[[i]]$M))
      detSigma <- detSigma - 2 * detQi - nrow(mod$items[[i]]$Q) * detMi
    }
    detSigma <- detSigma - sum(log(diag(mod$D)))
    ll_val <- ll_val -1/2 * detSigma
    return(ll_val)
  }
  return(fn_LL)
}

#' ---------------------------------------------------------------------
#' Utility functions to implement LMOK with SPDE approach
#' ---------------------------------------------------------------------

SPDE_eval_diff <- function(u,v) {
  tab <- matrix(c(mean(u - v), mean(abs(u-v)), sqrt(mean((u-v)^2))), nrow = 1, ncol = 3)
  colnames(tab) <- c("mean", "MAE", "RMSE")
  tab
}

#' Definition of the parameters of the LMOK
SPDE_param2value <- function(param, nu = 1, margin = 0.2, nodes = c(100, 100), type = "IntFac"){
  ll <- NULL
  if(type == "IntFac"){
    ll <- list(
      nodes   = nodes,
      margin  = margin,
      nu      = nu,
      type    = type,
      range   = param[5],
      sigma_1 = abs(param[1]),
      sigma_2 = abs(sqrt(param[2]^2 + param[3]^2)),
      tau     = param[2] / sqrt(param[2]^2 + param[3]^2),
      sigma   = abs(param[4])
    )
  }
  if(type == "IndFac"){
    ll <- list(
      nodes   = nodes,
      margin  = margin,
      nu      = nu,
      type    = type,
      sigma_1 = abs(param[1]),
      sigma_2 = abs(param[2]),
      sigma   = abs(param[3]),
      range_1 = param[4],
      range_2 = param[5]
    )
  }
  ll
}

SPDE_value2param <- function(value){
  par <- NULL
  if(value$type == "IntFac"){
    par <- c(
      value$sigma_1, 
      value$tau * value$sigma_2, 
      value$sigma_2 * sqrt(1-value$tau^2), 
      value$sigma, 
      value$range
    )
  }
  if(value$type == "IndFac"){
    par <- c(
      value$sigma_1, 
      value$sigma_2, 
      value$sigma, 
      value$range_1,
      value$range_2
    )
  }
  par  
}

# Create a RGeostats model from a parameters structure
SPDE_param2model <- function(param){
  m <- NULL
  if(param$type == "IndFac"){
    m <- model.create(vartype = 8, 
                      param = param$nu, 
                      range = param$range_1,
                      sill  = matrix(c(param$sigma_1^2, rep(0,3)),2,2)
    )
    m <- model.create(vartype = 8, 
                      param = param$nu, 
                      range = param$range_2,
                      sill = matrix(c(rep(0,3),param$sigma_2^2),2,2),
                      model = m
    )
  } else if(param$type == "IntFac") {
    covM <- matrix(c(param$sigma_1^2,
                     rep(param$sigma_1*param$sigma_2*param$tau,2),
                     param$sigma_2^2), 2, 2)
    m <- model.create(vartype = 8, 
                      param = param$nu, 
                      range = param$range,
                      sill  = covM
    )
  } else if(param$type == "bivar") {
    m <- NULL
    covM <- matrix(0, 4, 4)
    covM[1:2, 1:2] <- matrix(
      c(param$coef_1[1]^2, rep(prod(param$coef_1),2), param$coef_1[2]^2),
      2, 2)
    m <- model.create(vartype = 8,
                      param = param$nu,
                      range = param$range[1],
                      sill  = covM)
    covM <- matrix(0, 4, 4)
    covM[3:4, 3:4] <- matrix(
      c(param$coef_2[1]^2, rep(prod(param$coef_2),2), param$coef_2[2]^2),
      2, 2)
    m <- model.create(vartype = 8,
                      param = param$nu,
                      range = param$range[2],
                      sill  = covM, 
                      model = m)
  }
  m
}
