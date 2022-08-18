#'-------------------------------
#' Generic function to compute the K-fold cross validation with factors
#'-------------------------------
#' Two functions compute and results
#'  - "kfold_compute" computes the xval on the input data base and returns it
#'  - "kfold_results" presents different figures and return the score table
#'  
#'  Required in dbin, locators defined "z", "f1", "code" (f0 is the constant 1.0)
#'  Optional in dbin, a selection defined (locator "sel")
#'  
#'  Example of a function to use SPDE_LMOK_krigsim : model_TBD, mesh_TBD
#'  
# fn_estim_SPDE_LMOK <- function(dbin, dbout, radix = "estim_SPDE_LMOK"){
#   SPDE_LMOK_krigsim(dbin, dbout, 
#                     model = model_TBD, 
#                     mesh = mesh_TBD,
#                     nsim = 0, radix = radix)
# }
#'

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

kfold_results <- function(dbin, name_real, name_esti, 
                          stats = c("number", "mean", "MAE", "RMSE"),
                          radix = "K-fold", flag.fig = TRUE){
  
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
    correlation(dbin, name2 = name_real, name1 = name_esti,
                title = paste0("Cross validation of ", name_real, " - Kriging"),
                xlab  = "Estimated value", ylab = "Actual value",
                flag.iso = TRUE, flag.aspoint = ifelse(dbin$nactive > 500, FALSE, TRUE),
                flag.diag = TRUE, diag.col = "red", diag.lwd = 2,
                flag.ce = TRUE, ce.col = "blue"
    )
    # figure : residual vs. estimated value
    dbin <- db.add(dbin, Residual = dbin[,name_real] - dbin[,name_esti])
    correlation(dbin, name2 = "Residual", name1 = name_esti,
                title = paste0("Cross validation of ", name_real, " - Kriging"),
                xlab  = "Estimated value", ylab = "Residual",
                flag.iso = FALSE, flag.aspoint = ifelse(dbin$nactive > 500, FALSE, TRUE),
                flag.diag = FALSE,
                flag.ce = TRUE, ce.col = "blue"
    )
    abline(h = 0.0, col = "red", lwd = 2)  
  }
  tab 
}


#'-------------------------------
#' Compute the K-fold cross validation (with two factors)
#'-------------------------------
#' Two functions compute and results
#'  - "SPDE_kfold_compute" computes the xval on the input data base and returns it
#'  - "SPDE_kfold_results" presents different figures and return the score table
#'  
#'  Required in dbin, locators defined "z", "f1", "code" (f0 is the constant 1.0)
#'  Optional in dbin, a selection defined (locator "sel")

SPDE_LMOK_kfold_compute <- function(dbin, model,  mesh = NA, 
                                    radix = "SPDE_LMOK.kfold", verbose = TRUE){
  idx_in_z <- db.getcols(dbin, loctype = "z")
  idx_in_f <- db.getcols(dbin, loctype = "f")
  idx_in_c <- db.getcols(dbin, loctype = "code")
  idx_in_s <- db.getcols(dbin, loctype = "sel")
  stopifnot(length(idx_in_z) == 1)
  stopifnot(length(idx_in_f) == model$nvar)
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
                 names(dbin@items)[idx_in_z], " with factor ",
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
    res <- SPDE_LMOK_krigsim(
      dbin  = db.sel(dbin, sel_in),
      dbout = db.sel(dbin, sel_out), 
      model = model, mesh = mesh, nsim = 0, radix = radix)
    
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


#'-------------------------------
#' Compute the K-fold cross validation (with two factors)
#'-------------------------------
#' Two functions compute and results
#'  - "SPDE_kfold_compute" computes the xval on the input data base and returns it
#'  - "SPDE_kfold_results" presents different figures and return the score table
#'  
#'  Required in dbin, locators defined "z", "f1", "code" (f0 is the constant 1.0)
#'  Optional in dbin, a selection defined (locator "sel")

#'-------------------------------
#' Compute the K-fold cross validation (with two factors)
#'-------------------------------
#' Two functions compute and results
#'  - "SPDE_kfold_compute" computes the xval on the input data base and returns it
#'  - "SPDE_kfold_results" presents different figures and return the score table
#'  
#'  Required in dbin, locators defined "z", "f1", "code" (f0 is the constant 1.0)
#'  Optional in dbin, a selection defined (locator "sel")

SPDE_LMOK_kfold_compute <- function(dbin, model,  mesh = NA, 
                               radix = "SPDE_LMOK.kfold", verbose = TRUE){
  idx_in_z <- db.getcols(dbin, loctype = "z")
  idx_in_f <- db.getcols(dbin, loctype = "f")
  idx_in_c <- db.getcols(dbin, loctype = "code")
  idx_in_s <- db.getcols(dbin, loctype = "sel")
  stopifnot(length(idx_in_z) == 1)
  stopifnot(length(idx_in_f) == model$nvar)
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
                 names(dbin@items)[idx_in_z], " with factor ",
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
    res <- SPDE_LMOK_krigsim(
      dbin  = db.sel(dbin, sel_in),
      dbout = db.sel(dbin, sel_out), 
      model = model, mesh = mesh, nsim = 0, radix = radix)
    
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

# Old interface with param instead of model
SPDE_kfold_compute <- function(dbin, param,  mesh = NA, 
                               radix = "SPDE.kfold", verbose = TRUE){
  idx_in_z <- db.getcols(dbin, loctype = "z")
  idx_in_f <- db.getcols(dbin, loctype = "f")
  idx_in_c <- db.getcols(dbin, loctype = "code")
  idx_in_s <- db.getcols(dbin, loctype = "sel")
  stopifnot(length(idx_in_z) == 1)
  stopifnot(length(idx_in_f) == 1)
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
                 names(dbin@items)[idx_in_z], " with factor ",
                 names(dbin@items)[idx_in_f]
                 )
          )
    print(paste0("K-fold: K = ", K))
  }
  
  Z0_estim <- rep(NaN, dbin$nech)
  Z1_estim <- rep(NaN, dbin$nech)
  
  for (c in l_folds){
    sel_in  <- sel & (!is.na(folds))&(folds != c)
    sel_out <- sel & (!is.na(folds))&(folds == c)
    
    if(verbose) {
      print(paste0("K-fold: processing fold #",c))
    }
    res <- SPDE_krigsim(
      dbin  = db.sel(dbin, sel_in),
      dbout = db.sel(dbin, sel_out), 
      param = param, mesh = mesh, nsim = 0, radix = "SPDE_kfold")
    Z0_estim[sel_out] <- db.extract(res, "SPDE_kfold.Z0.estim", flag.compress = TRUE)
    Z1_estim[sel_out] <- db.extract(res, "SPDE_kfold.Z1.estim", flag.compress = TRUE)
  }
  
  S_estim <- Z0_estim + Z1_estim * dbin[,idx_in_f]
  db.add(dbin, names = paste(radix, "Z", "estim", sep = "."), S_estim)
}

#' Function to display the results of the cross validation stored in a db
#' stats gives the list of statistics to be computed
#'  number, mean, MAE, RMSE, min, max, median
SPDE_LMOK_kfold_results <- function(dbin, name_real, name_esti, stats = c("number", "mean", "MAE", "RMSE"),
                               radix = "kfold", flag.fig = TRUE){
  SPDE_kfold_results(
    dbin = dbin, 
    name_real = name_real, 
    name_esti = name_esti, 
    stats = stats,
    radix = radix, 
    flag.fig = flag.fig)
}
  

SPDE_kfold_results <- function(dbin, name_real, name_esti, stats = c("number", "mean", "MAE", "RMSE"),
                               radix = "K-fold", flag.fig = TRUE){
res <- as.matrix(db.extract(dbin, names = c(name_real, name_esti), flag.compress = TRUE))
#stopifnot(dim(res)[2] == 2)
# functions to be defined
number <- function(u){sum(!is.na(u))}
MAE    <- function(u){mean(abs(u), na.rm = TRUE)}
RMSE   <- function(u){sqrt(mean(u^2, na.rm = TRUE))}

tab <- matrix(NaN, nrow = 1, ncol = length(stats))
colnames(tab) <- stats
for (i in seq_along(stats)){
  tab[1,i] <- apply(X = as.matrix(res[,1] - res[,2]),  MARGIN = 2, FUN = stats[i])
}

if(flag.fig){
  constant.define("asp",1)
  # figure : actual value vs. estimated value
  correlation(dbin, name2 = name_real, name1 = name_esti,
              title = paste0("Cross validation of ", name_real, " - Kriging"),
              xlab  = "Estimated value", ylab = "Actual value",
              flag.iso = TRUE, flag.aspoint = ifelse(dbin$nactive > 500, FALSE, TRUE),
              flag.diag = TRUE, diag.col = "red", diag.lwd = 2,
              flag.ce = TRUE, ce.col = "blue"
  )
  # figure : residual vs. estimated value
  dbin <- db.add(dbin, Residual = dbin[,name_real] - dbin[,name_esti])
  correlation(dbin, name2 = "Residual", name1 = name_esti,
              title = paste0("Cross validation of ", name_real, " - Kriging"),
              xlab  = "Estimated value", ylab = "Residual",
              flag.iso = FALSE, flag.aspoint = ifelse(dbin$nactive > 500, FALSE, TRUE),
              flag.diag = FALSE,
              flag.ce = TRUE, ce.col = "blue"
  )
  abline(h = 0.0, col = "red", lwd = 2)  
  }
tab 
}

#-------------------------------
# Compute the cross validation (with two factors)
#-------------------------------

OK_Xval <- function(db, sel_test, sel_train, param, radix = "OK.xval", 
                      flag.fig = FALSE, mesh = NA){
  GEN_Xval(method = "OK",
           db = db, 
           sel_test = sel_test, 
           sel_train = sel_train, 
           param = param, 
           mesh = mesh, 
           radix = radix, 
           flag.fig = flag.fig
  )
}


SPDE_Xval <- function(db, sel_test, sel_train, param, radix = "SPDE.xval", 
                      flag.fig = FALSE, mesh = NA){
  GEN_Xval(method = "SPDE",
           db = db, 
           sel_test = sel_test, 
           sel_train = sel_train, 
           param = param, 
           mesh = mesh, 
           radix = radix, 
           flag.fig = flag.fig
           )
}
  
# Compute the Xvalidation
GEN_Xval <- function(method, db, sel_test, sel_train, param, mesh = NA, radix = "GEN.xval", flag.fig = FALSE){
  
  idx_f   <- db.getcols(db, loctype = "f", rank.match = 1)  
  idx_z   <- db.getcols(db, loctype = "z", rank.match = 1)  
  stopifnot((idx_f > 0)&(idx_z > 0))
  
  nm_f <- names(db@items)[idx_f]
  nm_z <- names(db@items)[idx_z]
  
  dbin  <- db.sel(db, selold = sel_train)
  dbout <- db.sel(db, selold = sel_test)
  if (method == "SPDE"){
    dbout <- SPDE_krigsim(dbout, dbin, param = param, mesh = mesh, 
                          nsim = 0, seed = NA, flag.ce = FALSE,
                          radix = radix)
  } else if (method == "OK"){
    dbout <- OK_kriging(dbout, dbin, param = param, 
                        radix = radix)
  } else { 
    print(paste0(method, ", not available !"))
  }
  
  Z0 <- db.extract(dbout, names = paste(radix, "Z0", "estim", sep = "."), flag.compress = TRUE)
  Z1 <- db.extract(dbout, names = paste(radix, "Z1", "estim", sep = "."), flag.compress = TRUE)
  N  <- db.extract(dbout, names = nm_f, flag.compress = TRUE)
  S_real <- db.extract(dbout, names = nm_z, flag.compress = TRUE)
  S_esti <- Z0 + Z1 * N
  tab    <- matrix(unlist(SPDE_eval_diff(S_esti, S_real)), 1, 3)
  colnames(tab) <- c("mean", "MAE", "RMSE")
  rownames(tab) <- c(method)
  
  # Xplot
  if(flag.fig){
    nm_z_estim <-paste(radix, nm_z, "estim", sep = ".")
    dbout <- db.replace(dbout, name = "S_esti", tab = S_esti)
    dbout <- db.rename(dbout, names = "S_esti", newnames = nm_z_estim)
    constant.define("asp",1)
    correlation(dbout, name2 = nm_z, name1 = nm_z_estim,
                title = paste0("Cross validation of ", nm_z, " - Kriging"),
                xlab  = "Estimated value", ylab = "Actual value",
                flag.iso = TRUE, flag.aspoint = TRUE,
                flag.diag = TRUE, diag.col = "red", diag.lwd = 2,
                flag.ce = TRUE, ce.col = "blue"
    )
  }
  tab 
}

# --------------------
# MONO factor and MONO variable
# --------------------


OK_Xval1 <- function(db, sel_test, sel_train, param, radix = "OK.xval", 
                    flag.fig = FALSE, mesh = NA){
  GEN_Xval1(method = "OK",
           db = db, 
           sel_test = sel_test, 
           sel_train = sel_train, 
           param = param, 
           mesh = mesh, 
           radix = radix, 
           flag.fig = flag.fig
  )
}


SPDE_Xval1 <- function(db, sel_test, sel_train, param, radix = "SPDE.xval", 
                      flag.fig = FALSE, mesh = NA){
  GEN_Xval1(method = "SPDE",
           db = db, 
           sel_test = sel_test, 
           sel_train = sel_train, 
           param = param, 
           mesh = mesh, 
           radix = radix, 
           flag.fig = flag.fig
  )
}

# Compute the Xvalidation
GEN_Xval1 <- function(method, db, sel_test, sel_train, param, mesh = NA, radix = "GEN.xval", flag.fig = FALSE){
  
  idx_f   <- db.getcols(db, loctype = "f", rank.match = 1)  
  idx_z   <- db.getcols(db, loctype = "z", rank.match = 1)  
  stopifnot((idx_f > 0)&(idx_z > 0))
  
  nm_f <- names(db@items)[idx_f]
  nm_z <- names(db@items)[idx_z]
  
  dbin  <- db.sel(db, selold = sel_train)
  dbout <- db.sel(db, selold = sel_test)
  if (method == "SPDE"){
    dbout <- SPDE_krigsim1(dbout, dbin, param = param, mesh = mesh, 
                          nsim = 0, seed = NA, flag.ce = FALSE,
                          radix = radix)
  } else if (method == "OK"){
    dbout <- OK_kriging1(dbout, dbin, param = param, 
                        radix = radix)
  } else { 
    print(paste0(method, ", not available !"))
  }
  
  Z1 <- db.extract(dbout, names = paste(radix, "Z1", "estim", sep = "."), flag.compress = TRUE)
  N  <- db.extract(dbout, names = nm_f, flag.compress = TRUE)
  S_real <- db.extract(dbout, names = nm_z, flag.compress = TRUE)
  S_esti <- Z1 * N
  tab    <- matrix(unlist(SPDE_eval_diff(S_esti, S_real)), 1, 3)
  colnames(tab) <- c("mean", "MAE", "RMSE")
  rownames(tab) <- c(method)
  
  # Xplot
  if(flag.fig){
    nm_z_estim <-paste(radix, nm_z, "estim", sep = ".")
    dbout <- db.replace(dbout, name = "S_esti", tab = S_esti)
    dbout <- db.rename(dbout, names = "S_esti", newnames = nm_z_estim)
    constant.define("asp",1)
    correlation(dbout, name2 = nm_z, name1 = nm_z_estim,
                title = paste0("Cross validation of ", nm_z, " - Kriging"),
                xlab  = "Estimated value", ylab = "Actual value",
                flag.iso = TRUE, flag.aspoint = TRUE,
                flag.diag = TRUE, diag.col = "red", diag.lwd = 2,
                flag.ce = TRUE, ce.col = "blue"
    )
  }
  tab 
}


# ---------------------
# Bi-variable
# ---------------------

OK_Xval_bivar <- function(db, sel_test, sel_train, param, radix = "COK.xval", 
                    flag.fig = FALSE, mesh = NA){
  GEN_Xval_bivar(method = "OK",
           db = db, 
           sel_test = sel_test, 
           sel_train = sel_train, 
           param = param, 
           mesh = mesh, 
           radix = radix, 
           flag.fig = flag.fig
  )
}

SPDE_Xval_bivar <- function(db, sel_test, sel_train, param, radix = "SPDE.COK.xval", 
                      flag.fig = FALSE, mesh = NA){
  GEN_Xval_bivar(method = "SPDE",
           db = db, 
           sel_test = sel_test, 
           sel_train = sel_train, 
           param = param, 
           mesh = mesh, 
           radix = radix, 
           flag.fig = flag.fig
  )
}

GEN_Xval_bivar <- function(method, db, sel_test, sel_train, param, mesh = NA, radix = NA, flag.fig = FALSE){
  
  idx_f   <- db.getcols(db, loctype = "f", rank.match = 1)  
  idx_z   <- db.getcols(db, loctype = "z", rank.match = 1)  
  idx_c   <- db.getcols(db, loctype = "code", rank.match = 1)  
  stopifnot((idx_f > 0)&(idx_z > 0)&(idx_c > 0))
  
  nm_f <- names(db@items)[idx_f]
  nm_z <- names(db@items)[idx_z]
  nm_c <- names(db@items)[idx_c]
  
  dbin  <- db.sel(db, selold = sel_train)
  dbout <- db.sel(db, selold = sel_test)
  
  if (method == "SPDE") {
    dbout <- SPDE_cokrigsim(dbout, dbin, param = param, mesh = mesh, 
                            nsim = 0, seed = NA, flag.ce = FALSE,
                            radix = radix)
  } else if(method == "OK") {
    dbout <- OK_cokriging(dbout, dbin, param = param, 
                          radix = radix)
  } else {
    print(paste0(method, "not available!"))
  }
  
  # Actual values  
  Type <- db.extract(dbout, names = nm_c, flag.compress = TRUE)
  N    <- db.extract(dbout, names = nm_f, flag.compress = TRUE)
  S    <- db.extract(dbout, names = nm_z, flag.compress = TRUE)
  # Estimated values
  a.Z0 <- db.extract(dbout, names = paste(radix, "a.Z0", "estim", sep = "."), flag.compress = TRUE)
  a.Z1 <- db.extract(dbout, names = paste(radix, "a.Z1", "estim", sep = "."), flag.compress = TRUE)
  b.Z0 <- db.extract(dbout, names = paste(radix, "b.Z0", "estim", sep = "."), flag.compress = TRUE)
  b.Z1 <- db.extract(dbout, names = paste(radix, "b.Z1", "estim", sep = "."), flag.compress = TRUE)
  
  S_esti <- rep(NaN, length(S))
  sel_a  <- (Type == 1)
  S_esti[sel_a] <- a.Z0[sel_a] + a.Z1[sel_a] * N[sel_a]
  S_esti[!sel_a] <- b.Z0[!sel_a] + b.Z1[!sel_a] * N[!sel_a]
  tab    <- rbind(
    matrix(unlist(SPDE_eval_diff(S_esti[sel_a], S[sel_a])), 1, 3),
    matrix(unlist(SPDE_eval_diff(S_esti[!sel_a], S[!sel_a])), 1, 3),
    matrix(unlist(SPDE_eval_diff(S_esti, S)), 1, 3)
  )
  colnames(tab) <- c("mean", "MAE", "RMSE")
  rownames(tab) <- paste(method, " - ", c("(a)", "(b)", "all"))
  
  # Xplot
  if(flag.fig){
    nm_z_estim <-paste(radix, nm_z, "estim", sep = ".")
    dbout <- db.replace(dbout, name = "S_esti", tab = S_esti)
    dbout <- db.rename(dbout, names = "S_esti", newnames = nm_z_estim)
    constant.define("asp",1)
    # (a) and (b)
    correlation(dbout, name2 = nm_z, name1 = nm_z_estim,
                title = paste0("Cross validation of ", nm_z, " - (a) and (b) - Cokriging"),
                xlab  = "Estimated value", ylab = "Actual value",
                flag.iso = TRUE, flag.aspoint = TRUE,
                flag.diag = TRUE, diag.col = "red", diag.lwd = 2,
                flag.ce = TRUE, ce.col = "blue"
    )
    # (a) Type == 1
    dbout <- db.sel(dbout, selold = sel_test)
    dbout <- db.sel(dbout, dbout@items[,idx_c] == 1, combine = "and", namesel = "xval_a")
    correlation(dbout, name2 = nm_z, name1 = nm_z_estim,
                title = paste0("Cross validation of ", nm_z, " - (a) - Cokriging"),
                xlab  = "Estimated value", ylab = "Actual value",
                flag.iso = TRUE, flag.aspoint = TRUE,
                flag.diag = TRUE, diag.col = "red", diag.lwd = 2,
                flag.ce = TRUE, ce.col = "blue"
    )
    # (b) Type == 2
    dbout <- db.sel(dbout, selold = sel_test)
    dbout <- db.sel(dbout, dbout@items[,idx_c] == 2, combine = "and", namesel = "xval_b")
    correlation(dbout, name2 = nm_z, name1 = nm_z_estim,
                title = paste0("Cross validation of ", nm_z, " - (b) - Cokriging"),
                xlab  = "Estimated value", ylab = "Actual value",
                flag.iso = TRUE, flag.aspoint = TRUE, 
                flag.diag = TRUE, diag.col = "red", diag.lwd = 2,
                flag.ce = TRUE, ce.col = "blue"
    )
  }
  tab 
}

