#' setting the input data file for the estimation using a linear model with locally varying parameters
#' 
#'@param dat The db-class structure containing the data file
#'@param nm_fac A string containing the name of the factors
#'@param nm_var A string containing the name of the variable to be estimated
#'@param nm_sel A string containing the name of the factors
#'@param nm_type A string containing the name of variable defining the type of sample
#'@param type A string defining the estimation problem. Possible values are:
#' - "house": factors for the surface of houses are 1 and nm_fac
#' - "flat" : factors for the surface of flats are 1 and nm_fac
#' - "both" : factors for the surface of houses and flats,nm are 1_H, 1_F, 1_H x nm_fac, 1_F x nm_fac
#'@param model A string defining the model to be used. Possible values are: 
#' - M1     : the coefficients are independent (type should be house or flat)
#' - M2     : the coefficients are correlated  (type should be house or flat)
#' - M3     : the coefficients are correlated  (type should be both)
#' - Mono   : a single coefficient is used.
#'@param sigma A numerical vectors with the standard deviation of the modeling error.
#'A single value if the mono-variate cases and two values in the bi-variate cases.
#'@return The input db where the role of each required variable has been setup.
prepareDb = function(dat,nm_fac,nm_var,nm_sel,nm_type,type,model,sigma=1)
{
  # re-initialisation de la base de données
  dat  <- db.locerase(dat, "z")
  dat  <- db.locerase(dat, "f")
  dat  <- db.locerase(dat, "code")
  dat  <- db.delete(dat, "sel")
  dat  <- db.delete(dat,"Un*")
  dat  <- db.delete(dat, "sigma2")
  dat  <- db.delete(dat,c("N.a*","N.b*","Un.a*","Un.b*"))
  
  # Selection et creation des variables
 
  if(type %in% c("house","flat")){
    if (model %in% c("M1", "M2")){
      dat <- db.add(dat, Un = 1.)
      dat <- db.locate(dat, name = c("Un", nm_fac), loctype = "f")
    } else if (model %in% c("Mono")){
      dat <- db.locate(dat, name = c(nm_fac), loctype = "f")
    } else if (model %in% c("OK")){
      dat <- db.add(dat, Un = 1.)
      dat <- db.locate(dat, name = c("Un"), loctype = "f")
    } else {
      stop(paste0("Not allowed model = ", model))
    }
    dat <- db.add(dat, sigma2 = sigma^2, loctype = "v")
  } else  if(type %in% c("both")){
    if (model %in% c("M3")){
      Un.a <- as.numeric(dat[,nm_type] == 1)
      Un.b <- as.numeric(dat[,nm_type] == 2)
      N.a  <- dat[,nm_fac]; N.a[dat[,nm_type] != 1] <- 0.0 
      N.b  <- dat[,nm_fac]; N.b[dat[,nm_type] != 2] <- 0.0 
      dat <- db.add(dat, Un.a, Un.b, N.a, N.b,loctype="f")
    } else if(model %in% "Mono") {
      N.a  <- dat[,nm_fac]; N.a[dat[,nm_type] != 1] <- 0.0 
      N.b  <- dat[,nm_fac]; N.b[dat[,nm_type] != 2] <- 0.0 
      dat <- db.add(dat, N.a, N.b,loctype="f")
    } else if(model %in% "OK") {
      Un.a <- as.numeric(dat[,nm_type] == 1)
      Un.b <- as.numeric(dat[,nm_type] == 2)
      dat <- db.add(dat, Un.a, Un.b,loctype="f")
    } else {
      stop(paste0("Not allowed model = ", model))
    }
    dat_sigma <- sigma[dat[, nm_type]]
    dat <- db.add(dat, sigma2 = dat_sigma^2, loctype = "v")
    dat <- db.locate(dat, names = nm_type, loctype = "code", flag.locnew = TRUE)
  }

  # selection
  if(type == "house"){
    tempsel <- (dat[,nm_type] == 1)
  } else if (type == "flat"){
    tempsel <- (dat[,nm_type] == 2)
  } else if (type == "both"){
    tempsel <- rep(1, dat$nech) 
  } else {
    stop(paste0("Not allowed type = ", type))
  }
  dat = db.add(dat, tempsel = tempsel)
  dat = db.rename(dat,"tempsel",nm_sel)
  dat = db.locate(dat,nm_sel,"sel")
  
  # target variable     
  dat <- db.locate(dat, name = c(nm_var), loctype = "z", flag.locnew = TRUE)
  
  return(dat)
}


filter = function(dat,selection)
{
  temp = db.reduce(db.sel(dat,selection))
  temp = db.delete(temp,"sel")
  temp$locators=dat$locators
  return(temp)
}

filter_dat = function(dat,filter_list)
{
  sel = rep(TRUE,dat$nech)
  for(selmaker in filter_list)
  {
    sel = sel & selmaker(dat) 
  }
  
  filter(dat,sel)
}

