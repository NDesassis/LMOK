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
#' - "both" : factors for the surface of houses and flats, are 1_H, 1_F, 1_H x nm_fac, 1_F x nm_fac
#' - mono_house: factor for the surface of houses is nm_fac
#' - mono_flat:  factor for the surface of flats is nm_fac
#' - mono_bivar: factor for the surface of houses and flats, 1_H x nm_fac, 1_F x nm_fac
#'@param sigma A numerical vectors with the standard deviation of the modeling error.
#'A single value if the mono-variate cases and two values in the bi-variate cases.
#'@return The input db where the role of each required variable has been setup.
prepareDb = function(dat,nm_fac,nm_var,nm_sel,nm_type,type,sigma=1)
{
  # re-initialisation de la base de données
  dat  <- db.locerase(dat, "z")
  dat  <- db.locerase(dat, "f")
  dat  <- db.locerase(dat, "code")
  dat  <- db.delete(dat, "sel")
  dat  <- db.delete(dat,"Un*")
  dat  <- db.delete(dat, "sigma2")
  # Selection des variables
 
  if(type%in%c("flat","house"))
  {
    tempsel = dat[,"C1NATURE"] == ifelse(type=="flat",2,1)
    dat = db.add(dat, Un = 1.)
    dat <- db.locate(dat, name = c("Un", nm_fac), loctype = "f")
    dat <- db.add(dat, sigma2 = sigma^2, loctype = "v")
  } else if (type %in% c("mono_flat", "mono_house")){
    
    if (type == "mono_house"){
      tempsel = dat[,"C1NATURE"] == 1
    } else {
      tempsel = dat[,"C1NATURE"] == 2
    } 
    dat <- db.locate(dat, name = c(nm_fac), loctype = "f")
    dat <- db.add(dat, sigma2 = sigma^2, loctype = "v")
  }else  if(type %in% c("both"))
  {
    tempsel = rep(1,dat$nech)
    dat = db.delete(dat,c("N.a*","N.b*","Un.a*","Un.b*"))
    Un.a <- as.numeric(dat[,nm_type] == 1)
    Un.b <- as.numeric(dat[,nm_type] == 2)
    N.a  <- dat[,nm_fac]; N.a[dat[,nm_type] != 1] <- 0.0 
    N.b  <- dat[,nm_fac]; N.b[dat[,nm_type] != 2] <- 0.0 
    dat <- db.add(dat, Un.a, Un.b, N.a, N.b,loctype="f")
    dat <- db.locate(dat, names = nm_type, loctype = "code", flag.locnew = TRUE)
    dat_sigma <- sigma[dat[, nm_type]]
    dat  <- db.add(dat, sigma2 = dat_sigma^2, loctype = "v")

  }else  if(type %in% c("mono_bivar"))
  {
    tempsel = rep(1,dat$nech)
    dat = db.delete(dat,c("N.a*","N.b*"))
    N.a  <- dat[,nm_fac]; N.a[dat[,nm_type] != 1] <- 0.0 
    N.b  <- dat[,nm_fac]; N.b[dat[,nm_type] != 2] <- 0.0 
    dat <- db.add(dat, N.a, N.b,loctype="f")
    dat <- db.locate(dat, names = nm_type, loctype = "code", flag.locnew = TRUE)
    dat_sigma <- sigma[dat[, nm_type]]
    dat  <- db.add(dat, sigma2 = dat_sigma^2, loctype = "v")
  }
  # selection
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

