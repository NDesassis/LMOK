# cross validation functions

#' compute the k-fold cross validation using a linear model with locally varying parameters
#' 
#' The estimator of the parameters is provided by a generic function to be defined before calling 
#' the cross validation function. The cross validation is performed on each fold using
#' the other folds as input data. First the parameters are estimated, 
#' then the linear combination of the factors by the estimated coefficients
#' is computed for each selected sample of the input file.
#' 
#'@param dbin The db-class structure containing the data file
#'#'The following variables should be defined in the data file:
#' - the observed values (locator "z")
#' - the factors associated to the latent fields to be estimated (locator "f")
#' - the standard deviation of the modeling error (locator "v")
#' - the variable defining the folds (locator "c")
#' The number of latent fields/spatially varying parameters to be estimated is given by the number of factors
#' specified in the input data file (length(db.getcols(dbin, "f"))).
#' A selection on the input file may be active. 
#' Then, the k-fold cross validation
#' is computed only for the selected samples.
#'@param fn_estim A function performing the estimation of spatially varying parameters
#'on a target file from the input data file.
#'The return value is the target file with the estimated value added.
#'An example of such a function using the SPDE approach is,
#'  fn_estim <- function(dbin, dbout, radix = "estim_SPDE_LMOK"){
#'                     SPDE_LMOK_krigsim(dbout = dbout, dbin = dbin, 
#'                     model = model_TBD, 
#'                     mesh = mesh_TBD,
#'                     nsim = 0, radix = radix)
#'                     }
#' where the model and the mesh should be defined before the function.
#' Beware to the order of the input arguments if they are not named.
#'@param radix A string used as prefix of the variables created in return file.
#'@param verbose A Boolean variable to control the printed messages
#'@return The input db where the estimation of the observed variable (locator "z" in the input file) has been added.
kfold_compute <- function(dbin, fn_estim, radix = "kfold", verbose = TRUE){
  idx_in_z <- db.getcols(dbin, loctype = "z")
  idx_in_f <- db.getcols(dbin, loctype = "f")
  idx_in_c <- db.getcols(dbin, loctype = "code")
  idx_in_s <- db.getcols(dbin, loctype = "sel")
  stopifnot(length(idx_in_z) == 1)
  stopifnot(length(idx_in_f) >= 1)
  stopifnot(length(idx_in_c) == 1)
  
  folds <- db.extract(dbin, names = idx_in_c, flag.compress = FALSE)
  if (is.null(idx_in_s)) {sel <- rep(TRUE, dbin$nech)} 
  else {
    sel   <- (db.extract(dbin, names = idx_in_s) > 0)
  }
  
  l_folds <- sort(unique(folds[sel]))
  K <- length(l_folds)
  
  if(verbose) {
    print(paste0("K-fold: cross validation of ", 
                 names(dbin@items)[idx_in_z], " with factors ",
                 names(dbin@items)[idx_in_f]
      )
    )
    print(paste0("K-fold: K = ", K))
  }
  
  estim <- list()
  for (i in seq_along(idx_in_f)){
    estim[[1+length(estim)]] <- list(Z = rep(NaN, dbin$nech))
  }
  
  for (c in l_folds){
    sel_in  <- sel & (!is.na(folds))&(folds != c)
    sel_out <- sel & (!is.na(folds))&(folds == c)
    
    if(verbose) {
      print(paste0("K-fold: processing fold #",c))
    }
    res <- fn_estim(
      dbin  = db.sel(dbin, sel_in),
      dbout = db.sel(dbin, sel_out),
      radix = radix)
    
    for (i in seq_along(idx_in_f)){
      estim[[i]]$Z[sel_out] <- db.extract(res, paste(radix, paste0("Z", i), "estim", sep = "."), 
                                          flag.compress = TRUE)
    }
  }
  S_estim <- rep(0.0, dbin$nech)
  for (i in seq_along(idx_in_f)){
    S_estim <- S_estim + estim[[i]]$Z * dbin[,idx_in_f[i]]
    dbin <- db.add(dbin, 
                   names = paste(radix, paste0("Z", i), "estim", sep = "."),
                   estim[[i]]$Z)
  }
  db.add(dbin, names = paste(radix, "S", "estim", sep = "."), S_estim)
}

#' compute the scores of a cross validation
#' 
#' 
#'@param dbin The db-class structure containing the data file with the actual values and their estimated values.
#'@param name_real A string giving the name of the variable containing the actual values.
#'@param name_esti A string giving the name of the variable containing the estimated values.
#'@param stats a list of statistics to be computed on the cross validation error. The list actually the name of the functions.
#'Example of statistics are:
#' - number,  number of samples for which the cross validation error is defined
#' - mean,    the mean of error
#' - MAE,     the mean absolute error
#' - RMSE,    the square root of the squared error
#' - sd,      the standard deviation of the error,
#'@param flag.fig A Boolean variable controlling the display of two plots:
#' - actual value vs. estimated value
#' - residual vs. estimated value
#'@param plot_fn A list of functions used A Boolean variable controlling the display of the figures
#'  Used only if flag.fig is equal to TRUE
#' - start() to initialize the output of the graphic
#' - end()   to close the output graphic
#'@return The table containing the computed statistics for each fold and for the total input file.
kfold_results <- function(dbin, name_real, name_esti, 
                          stats = c("number", "mean", "MAE", "RMSE"),
                          flag.fig = TRUE, plot_fn = NA){
  
  res <- as.matrix(db.extract(dbin, names = c(name_real, name_esti), flag.compress = TRUE))
  idx_in_c <- db.getcols(dbin, loctype = "code")
  K <- 0
  if (!is.null(idx_in_c)){
    codes    <- db.extract(dbin, names = idx_in_c, flag.compress = TRUE)
    ll_codes <- sort(unique(codes))
    K <- length(ll_codes)
  }
  
  # functions to be defined
  number <- function(u){sum(!is.na(u))}
  MAE    <- function(u){mean(abs(u), na.rm = TRUE)}
  RMSE   <- function(u){sqrt(mean(u^2, na.rm = TRUE))}
  
  tab <- matrix(NaN, nrow = 1+K, ncol = length(stats))
  colnames(tab) <- stats
  rownames(tab) <- c(paste0(names(dbin@items)[idx_in_c], "==", ll_codes), "Total")
  if (K > 0){
    for (k in seq_along(ll_codes)){
      sel <- (codes == ll_codes[k])
      for (i in seq_along(stats)){
        tab[k,i] <- apply(X = as.matrix(res[sel,1] - res[sel,2]),  MARGIN = 2, FUN = stats[i])
      }
    }
  }
  # All selected samples
  for (i in seq_along(stats)){
    tab[1+K,i] <- apply(X = as.matrix(res[,1] - res[,2]),  MARGIN = 2, FUN = stats[i])
  }
  
  if(flag.fig){
    constant.define("asp",1)
    # figure : actual value vs. estimated value
    if(
      (class(plot_fn) == "list")&
      (plot_fn$start(file_name = paste(name_esti, "real_vs_estim", sep = "_"), radix = "xval"))) {
      correlation(dbin, name2 = name_real, name1 = name_esti,
                  title = paste0("Cross validation of ", name_real, " by ", name_esti),
                  xlab  = "Estimated value", ylab = "Actual value",
                  flag.iso = TRUE, flag.aspoint = ifelse(dbin$nactive > 500, FALSE, TRUE),
                  flag.diag = TRUE, diag.col = "red", diag.lwd = 2,
                  flag.ce = TRUE, ce.col = "blue"
      )
      plot_fn$end()
    }
    # figure : actual value vs. estimated value
    
    # figure : residual vs. estimated value
    if(
      (class(plot_fn) == "list")&
      (plot_fn$start(file_name = paste(name_esti, "residual_vs_estim", sep = "_"), radix = "xval"))) {
      dbin <- db.add(dbin, Residual = dbin[,name_real] - dbin[,name_esti])
      correlation(dbin, name2 = "Residual", name1 = name_esti,
                  title = paste0("Cross validation of ", name_real, " by ", name_esti),
                  xlab  = "Estimated value", ylab = "Residual",
                  flag.iso = FALSE, flag.aspoint = ifelse(dbin$nactive > 500, FALSE, TRUE),
                  flag.diag = FALSE,
                  flag.ce = TRUE, ce.col = "blue"
      )
      abline(h = 0.0, col = "red", lwd = 2)  
      plot_fn$end()
    }
    # figure : residual vs. estimated value
    
  } # flag.fig == TRUE
  tab 
}

#' compute the k-fold cross validation using the SPDE approach for the estimation of 
#' a linear model with locally varying parameters
#' 
#'@param dbin The db-class structure containing the data file
#'#'The following variables should be defined in the data file:
#' - the observed values (locator "z")
#' - the factors associated to the latent fields to be estimated (locator "f")
#' - the standard deviation of the modeling error (locator "v")
#' - the variable defining the folds (locator "c")
#' The number of latent fields/spatially varying parameters to be estimated is given by the number of factors
#' specified in the input data file (length(db.getcols(dbin, "f"))).
#' A selection on the input file may be active. 
#' Then, the k-fold cross validation
#' is computed only for the selected samples.
#'@param model The model-class defining the covariance model of the latent fields
#'As the SPDE approach is used, only K-Bessel covariances should be used as basic structures
#'of the linear model of coregionalization. 
#'#'@param mesh The mesh-class containing the meshing of the domain 
#'If the mesh is not provided (mesh == NA), the mesh is a regular grid computed using 
#'the bounding box of the input data set and the following parameters
#'  - param$nodes The number of vertices along the main directions
#'  - param$margin The inflating coefficients of the bounding box
#'@param radix A string used as prefix of the variables created in return file.
#'@param verbose A Boolean variable to control the printed messages
#'@return The input db where the estimation of the observed variable (locator "z" in the input file) has been added.
SPDE_LMOK_kfold_compute <- function(dbin, model,  mesh = NA, 
                                      radix = "SPDE_LMOK.kfold", verbose = TRUE){
  fn_estim <- function(dbin, dbout, radix = "kfold_SPDE_LMOK"){
    SPDE_LMOK_krigsim(dbout = dbout, dbin = dbin,
                      model = model,
                      mesh  = mesh,
                      nsim  = 0, radix = radix, verbose = verbose)
  }
  kfold_compute(dbin, fn_estim = fn_estim, radix = "kfold",  verbose = verbose)
}
  
