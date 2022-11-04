# functions to implement LMOK with the SPDE approach

#' compute the estimation/simulation of the locally varying parameters
#' 
#' The SPDE approach is used to compute the estimations of the latent fields.
#' 
#'@param dbout The db-class structure containing the target file
#'The estimation is performed only for the selected samples of the target file.
#'@param dbin The db-class structure containing the data file
#'#'The following variables should be defined in the data file:
#' - the observed values (locator "z")
#' - the factors associated to the latent fields to be estimated (locator "f")
#' - the standard deviation of the modeling error (locator "v")
#' The number of latent fields/factors should be equal to the number of variables 
#' of the model (length(db.getcols(dbin, "f")) == model$nvar).
#'@param model The model-class defining the covariance model of the latent fields
#'As the SPDE approach is used, only K-Bessel covariances should be used as basic structures
#'of the linear model of coregionalization. 
#'@param mesh The mesh-class containing the meshing of the domain 
#'If the mesh is not provided (mesh == NA), the mesh is a regular grid computed using 
#'the bounding box of the input data set and the following parameters
#'  - param$nodes The number of vertices along the main directions
#'  - param$margin The inflating coefficients of the bounding box
#'@param nsim The number of conditional simulations of the latent fields 
#'(i.e. the spatially varying parameters). 
#'If nsim == 0, only the kriging value is computed.
#'@param seed The seed used by the generator of the random numbers.
#'If seed == NA, the generator is not initialized. 
#'The random generator should be reset to be able to reproduce simulations. 
#'@param flag.ce A Boolean value controlling the use of the simulations (nsim > 0).
#'If flag.ce == TRUE, the mean and the standard deviation of the simulations are computed.
#'If flag.ce == FALSE, the simulations are stored in output file *dbout*.
#'@param radix A string used as prefix of the variables created in the output file *dbout*.
#'@param verbose A Boolean variable to control the printed messages
#'@return The target Db where the following variables have been added:
#'If nsim == 0,
#' - the estimation of the latent fields (if nsim==0)
#' If nsim > 0 and flag.ce == FALSE,
#' - the conditional simulations of the latent fields (if nsim>=0 and flag.ce == FALSE)
#' If nsim > 0 and flag.ce == TRUE,
#' - the conditional expectation computed as the mean of the simulations
#' - the conditional standard deviation computed as the standard deviation of simulations
#' These variables are multiplied for each one of the latent fields (numbered from 1).
SPDE_LMOK_krigsim <- function(dbout, dbin, model, mesh = NA,
                              nsim = 0, seed = NA, flag.ce = FALSE, 
                              radix = "SPDE_LMOK", verbose = FALSE){
  
  # measurement error
  sigma <- sqrt(db.extract(dbin, db.getname(dbin, "v", 1), flag.compress = TRUE))

  geo  <- SPDE_LMOK_init_geometry(dbin, mesh, verbose)
  mod  <- SPDE_LMOK_init_model(model, sigma, geo, verbose)
  spde <- SPDE_LMOK_init_spde(geo, mod, verbose)
  
  SPDE <- list(
           geometry = geo,
           model    = mod,
           spde     = spde
  )
  
  # computing the results on the output file (if necessary)
  if (!is.na(dbout)){
    res <- SPDE_LMOK_compute(model = SPDE,dbout = dbout,
                             nsim = nsim, seed = seed, flag.ce = flag.ce,
                             radix = radix, verbose = verbose
                            )
  } else {
    res <- SPDE
  }
  return(res)
}

# functions to compute LMOK on a output db from the SPDE model

#' compute the estimation/simulation of the locally varying parameters
#' 
#' The SPDE approach is used to compute the estimations of the latent fields.
#' 
#'@param model The a SPDE model computed by SPDE_LMOK_krigsim with dbout == NA.
#'It is a list with the folloging attributes:
#' - geometry  This is the output of the function SPDE_LMOK_init_geometry.
#' - model     This is the output of the function SPDE_LMOK_init_model
#' - spde      This is the output of the function SPDE_LMOK_init_spde
#'@param dbout The db-class structure containing the target file
#'@param nsim The number of conditional simulations of the latent fields 
#'(i.e. the spatially varying parameters). 
#'If nsim == 0, only the kriging value is computed.
#'@param seed The seed used by the generator of the random numbers.
#'If seed == NA, the generator is not initialized. 
#'The random generator should be reset to be able to reproduce simulations. 
#'@param flag.ce A Boolean value controlling the use of the simulations (nsim > 0).
#'If flag.ce == TRUE, the mean and the standard deviation of the simulations are computed.
#'If flag.ce == FALSE, the simulations are stored in output file *dbout*.
#'@param radix A string used as prefix of the variables created in the output file *dbout*.
#'@param verbose A Boolean variable to control the printed messages
#'@return The target Db where the following variables have been added:
#'If nsim == 0,
#' - the estimation of the latent fields (if nsim==0)
#' If nsim > 0 and flag.ce == FALSE,
#' - the conditional simulations of the latent fields (if nsim>=0 and flag.ce == FALSE)
#' If nsim > 0 and flag.ce == TRUE,
#' - the conditional expectation computed as the mean of the simulations
#' - the conditional standard deviation computed as the standard deviation of simulations
#' These variables are multiplied for each one of the latent fields (numbered from 1).
SPDE_LMOK_compute <- function(model, dbout,
                              nsim = 0, seed = NA, flag.ce = FALSE, 
                              radix = "SPDE_LMOK", verbose = FALSE){

  # initialization from the input model
  geo  <- model$geometry
  mod  <- model$model
  spde <- model$spde
  
  print(paste0("length(spde$Z)     = ", length(spde$Z)))
  print(paste0("length(spde$coeff) = ", length(spde$coeff)))
  print(paste0("geo$mesh$npoints  = ", length(geo$mesh$npoints)))
  
  # check the consistency between the spde model and the mesh
  stopifnot(length(spde$Z)/length(spde$coeff) == geo$mesh$npoints)

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
      
      str(Z_sim)
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
          dbout <- db.replace(dbout, name = nm_sim, tab = Z_sim[,n])
        }
      }
    } # end of simulations
  } # loop on the factors
  dbout 
}

#' function to initialize the log likelihood function
#' 
#' It can be used before the optimization of the log likelihood using *optim*.
#' 
#'@param dbin The db-class structure containing the data file
#'@param param2model The input function converting a vector of parameters into a model. 
#'  The input parameters of this function are (x, verbose)
#'   - x is a vector of real values, 
#'   - verbose is a Boolean variable
#'  the return value list(model, sigma) where
#'   - model is a model-class defining a covariance function
#'   - sigma is a vector of positive values giving the standard deviation of the modeling error
#'@param verbose A Boolean variable to control the printed messages
#'@return a function computing the log likelihood with
#'  - the input parameter (vector of real values x, a Boolean variable verbose)
#'  - the output the log likelihood value for the parameter values x
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
      detSigma <- detSigma - nrow(mod$item[[i]]$M) * detQi - nrow(mod$items[[i]]$Q) * detMi
    }
    detSigma <- detSigma - sum(log(diag(mod$D)))
    ll_val <- ll_val -1/2 * detSigma
    return(-ll_val)
  }
  return(fn_LL)
}

#' function to initialize the geometry structure
#' 
#' This function is used by the functions 
#' *SPDE_LMOK_krigsim* and *SPDE_LMOK_init_MLL*
#' 
#'@param dbin The db-class structure containing the data file
#'@param mesh The mesh-class structure containing the meshing of the domain
#'@param verbose A Boolean variable to control the printed messages
#'@return The geometry structure containing 
#'  - the input data *Y*, 
#'  - the factors *X*, 
#'  - the mesh *mesh* and 
#'  - the projection matrix *Aproj*
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
  
  # Input data (active selection & factors & observation)
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

#' function to initialize the model structure
#' 
#' This function is used by the functions 
#' *SPDE_LMOK_krigsim* and *SPDE_LMOK_init_MLL*
#' 
#'@param model The model-class structure containing the covariance model
#'@param sigma The vector of standard deviation of the modeling error
#'@param geo The geometry structure generated by the function *SPDE_LMOK_init_geometry*
#'@param verbose A Boolean variable to control the printed messages
#'@return The list of items used for the computation of the MLL, kriging and simulation
SPDE_LMOK_init_model <- function(model, sigma, geo, verbose = FALSE){
  
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
  if(verbose){
    print(model)
  }
  cholPrec =  Cholesky(prec, LDL=FALSE)
  return(list(
    items    = items,
    D        = D,
    cholPrec = cholPrec
  ))
}

#' function to initialize the SPDE structure
#' 
#' This function is used by the functions 
#' *SPDE_LMOK_krigsim* and *SPDE_LMOK_init_MLL*
#' 
#'@param geo The geometry structure generated by the function *SPDE_LMOK_init_geometry*
#'@param mod The model geometry generated by the function *SPDE_LMOK_init_model*
#'@param verbose A Boolean variable to control the printed messages
#'@return The list of items used for the computation of the MLL, kriging and simulation
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


