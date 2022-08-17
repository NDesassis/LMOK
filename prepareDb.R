prepareDb = function(dat,nm_fac,nm_var,nm_sel,nm_type,type,sigma2=1)
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
    dat = db.add(dat,Un = 1.)
    dat <- db.locate(dat, name = c("Un",nm_fac), loctype = "f")
  }else
  {
    tempsel = rep(1,dat$nech)
    dat = db.delete(dat,c("N.a*","N.b*","Un.a*","Un.b*"))
    Un.a <- as.numeric(dat[,nm_type] == 1)
    Un.b <- as.numeric(dat[,nm_type] == 2)
    N.a  <- dat[,nm_fac]; N.a[dat[,nm_type] != 1] <- 0.0 
    N.b  <- dat[,nm_fac]; N.b[dat[,nm_type] != 2] <- 0.0 
    dat <- db.add(dat, Un.a, Un.b, N.a, N.b,loctype="f")
    dat <- db.locate(dat, names = nm_type, loctype = "code", flag.locnew = TRUE)
    if (length(sigma2) > 1){
      idx_in_c <- db.getcols(dat, loctype = "code", rank.match = 1)
      code <- db.extract(dbin, names = idx_in_c, flag.compress = TRUE)
      sigma2v <- sigma2[code]
    }
  }
  
  
  dat = db.add(dat, tempsel = tempsel)
  dat = db.rename(dat,"tempsel",nm_sel)
  dat = db.locate(dat,nm_sel,"sel")
                  
  dat <- db.locate(dat, name = c(nm_var), loctype = "z", flag.locnew = TRUE)

 dat = db.add(dat,sigma2v,loctype = "v")
  #dat <- db.locate(dat, name = c("sel_train"), loctype = "sel", flag.locnew = TRUE)
  return(dat)
}
