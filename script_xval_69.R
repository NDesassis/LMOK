# 30/11/2012 : Xavier Freulon - Creation
# Script for the computation of the RMSE Table

knitr::opts_chunk$set(echo = TRUE)
library(RGeostats)
library(Matrix)
rm(list=ls())

# functions
source(file = "./functions_spde.R")
source(file = "./functions_Xval.R")
source(file = "./functions_prepareDB.R")
source(file = "./functions_meshing.R")
source(file = "./filter_examples.R")

# General parameters
set.seed(seed=125875)  # Initialisation of the random generator (to reproduce simulations)
n.fold = 10            # Number of folds to compute the cross-validation

# zone à traiter
area_nm   <- "Rhone"   # Selection de la zone dans "BPL", "Creuse", "Oise", "Rhone"

#filters: a list which contains all the filters to apply. 
# Example of filter functions are given in the file "filter_examples.R"
# filter_list = list(good_cedantes)
filter_list = list(lyon_suburb)
# filter_list = list(good_cedantes, lyon_suburb)
# filter_list = list() # no filter

nm_input_RData <- paste0("../Data/CCR_", area_nm, ".Rdata")

# Parameters of the linear model
nu      <- 1          # regularity coefficient of the Matern covariance
nm_var  <- "S_cap"    # Name of the variable to be predicted and observed (surface)
nm_fac  <- "NB_cap"   # Name of the factor (number of rooms)
nm_nature <- "C1NATURE" # Name of the nature (1 = House / 2 = Flat)

# Options for the optimizer (maximization of the likelihood)
optim_maxit = 100
optim_method = "Nelder-Mead"

#------------------------------------------------------
# Defining the cases and the table to store the results
#------------------------------------------------------
case_nature    <- c("house", "house", "flat", "flat" )
case_factor  <- c("{N}",     "{1,N}",   "{N}",    "{1,N}" )
case_model   <- c(
        "Global LM",
        "Local LM (500)",
        "Local LM (50)",
        "SPDE M1 (regular grid)",
        "SPDE M1 (adaptive meshing)",
        "SPDE M2 (adaptive meshing)",
        "SPDE M3 (adpative meshing)"
      )

#------------------------------------------------------
# Selection of the case to be treated
#------------------------------------------------------
sel_case  <- c(TRUE, TRUE, TRUE, TRUE)
sel_model <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

#------------------------------------------------------
# Table of the RMSE
#------------------------------------------------------
res_RMSE <- matrix(NaN, nrow = length(case_model), ncol = length(case_nature))
c_name   <- c()
for (c in (1:length(sel_case))){
  c_name <- c(c_name, paste0(case_nature[c], "-", case_factor[c]))
}
colnames(res_RMSE) <- c_name
rownames(res_RMSE) <- case_model

#------------------------------------------------------
# Table of the MLL
#------------------------------------------------------
res_MLL <- matrix(NaN, nrow = length(case_model), ncol = length(case_nature))
c_name   <- c()
for (c in (1:length(sel_case))){
  c_name <- c(c_name, paste0(case_nature[c], "-", case_factor[c]))
}
colnames(res_MLL) <- c_name
rownames(res_MLL) <- case_model


#------------------------------------------------------
# Loading the input data
#------------------------------------------------------
load(file = nm_input_RData)
dat   <- filter_dat(dat,filter_list)
folds <- sample(x = 1:n.fold, size = dat$nech, replace = TRUE)

#------------------------------------------------------
# Meshing the domain for SPDE (according to the data selection)
#------------------------------------------------------

adapted_mesh <- compute_mesh(
  poly = dpts[[1]]$polygon, 
  dat  = dat, 
  nodes = c(50,50), 
  margin = c(10,10),
  lambda = 0.04, N = NA
)

print(paste0("Number of vertices in the meshing = ", adapted_mesh$npoints))

#------------------------------------------------------
# Processing the different cases
#------------------------------------------------------
for (m in seq_along(case_model)[sel_model]) {
  # -------------------------------------------
  # Defining the modeling method
  # -------------------------------------------
  if (m == 1) { # Global linear model
    source (file = "./functions_LM.R")
    fn_estim <- fn_estim_LM1
  }
  if (m == 2) { # Local linear model (N = 500)
    N <- 500
    source (file = "./functions_LM.R")
    fn_estim <- fn_estim_LM2
  }
  if (m == 3) { # Local linear model (N = 500)
    N <- 50
    source (file = "./functions_LM.R")
    fn_estim <- fn_estim_LM2
  }
  if (m >= 4) { # SPDE meshing
    nodes   <- c(100,100) # Number of nodes in the main directions
    margin  <- 0.2        # Proportion to extend the bounding box of the domain
    param   <- list(margin = margin, nodes = nodes, nu = nu)
    if (m == 4) {
      mesh    <- NA
    } else {
      mesh   = adapted_mesh
    }
  }

  for (c in seq_along(case_nature)[sel_case]){
    print(paste0("---> Processing case: ", case_nature[c], "+", case_factor[c], "+", case_model[m]))

    nature  <- case_nature[c]
    factor  <- case_factor[c]
    
    if(factor == "{1}") {
      model = "OK"
    } else if(factor == "{N}") {
        model = "Mono"
        } else if (factor == "{1,N}"){
      if (m <= 5)      { model = "M1"}
      else if (m == 6) { model = "M2"}
      else if (m == 7) { 
        model  = "M3"
        nature = "both"}
    }
    nm_sel <- paste0("sel_", nature)
    type   <- nature # for compatibility
    s_est  <- 1.0

    # Cross validation (for the last model m == 7, houses and flats are evaluated jointly)
    if ((m < 7)|((m == 7)&(c <= 2))){
    
      # Inference of the model SPDE by MLE
      if (m > 3) {
        
        # Initialisation of the SPDE model
        if (factor == "{1}") {
          if (nature %in% c("house", "flat")){
            source(file = paste0("./param_", "Mono",".R"))
          }
          if (nature %in% c("both")){
            source(file = paste0("./param_", "Mono", "Bivar",".R"))
          }
          
        } else if (factor == "{N}"){
          if (nature %in% c("house", "flat")){
            source(file = paste0("./param_", model,".R"))
          }
          if (nature %in% c("both")){
            source(file = paste0("./param_", "Mono", "Bivar",".R"))
          }
          
        } else if (factor == "{1,N}"){
          source(file = paste0("./param_", model,".R"))
        }
    
        # Parameter inference of the SPDE model by MLL
        paste0("**** PARAMETER INFERENCE ****")
        cont$maxit = optim_maxit
        datMLL <- prepareDb(dat, nm_fac, nm_var, nm_sel, nm_nature, nature, model)
        fn_MLL <- SPDE_LMOK_init_MLL(datMLL, param2model = modelMulti, verbose = FALSE,mesh=mesh)
        res    <- optim(par = x_ini, fn = fn_MLL, verbose = F, method=optim_method, control = cont)
        res_MLL[m,c] <- res$value
        m_est  <- modelMulti(res$par, verbose = FALSE)$model
        s_est  <- modelMulti(res$par, verbose = FALSE)$sigma
        # estimation function
        fn_estim <- function(dbin, dbout, radix = "kfold_LMOK"){
          SPDE_LMOK_krigsim(dbout = dbout, dbin = dbin, 
                            model = m_est,
                            mesh  = mesh,
                            nsim  = 0, radix = radix, verbose = FALSE)
        }
      }
      datXval <- prepareDb(dat,nm_fac,nm_var,nm_sel,nm_nature, nature, model, sigma = s_est)
      datXval <- db.add(datXval, folds, locnames = "code")
      datXval <- kfold_compute(datXval, fn_estim = fn_estim, radix = "kfold_LMOK", 
                               verbose = FALSE)
      datXval <- db.locate(datXval, "C1NATURE", loctype = "code")
      tabXval <- kfold_results(datXval, name_real = nm_var, name_esti = "kfold_LMOK.S.estim", 
                               stats = c("number", "min", "max", "mean", "MAE", "RMSE"),
                               flag.fig = FALSE, plot_fn = NA)
      if (m < 7) {
        rmse   <- tabXval[nrow(tabXval), ncol(tabXval)]
        res_RMSE[m,c] <- rmse
      } else { # m == 7
        res_RMSE[m,c]   <- tabXval[1, ncol(tabXval)] # Houses
        res_RMSE[m,c+2] <- tabXval[2, ncol(tabXval)] # Flats
      }
      
    } # case (m < 7)|((m == 7)&(c <= 2))
  } # loop over the cases (c in 1:4)
} # loop over the models (m in 1:7)

res_RMSE <- res_RMSE[sel_model, sel_case]
res_MLL  <- res_MLL[sel_model, sel_case]
#------------------------------------------------------
# Table of the results - RMSE
#------------------------------------------------------
knitr::kable(round(res_RMSE, 2), 
             caption = paste0("RMSE of the cross validation (K = ", n.fold,")"))
write.csv(res_RMSE, file = "LMOK_Xvalidation_69_table_of_RMSE.csv")

#------------------------------------------------------
# Table of the results - MLL
#------------------------------------------------------
knitr::kable(round(res_MLL, 2), 
             caption = paste0("MLL for the parameter inference"))
write.csv(res_MLL, file = "LMOK_Xvalidation_69_table_of_MLL.csv")
