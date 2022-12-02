# Regression linéaire globale
fn_estim_LM1 = function(dbin,dbout,radix = "kfold_LM")
{
  Y = db.extract(dbin,db.getname(dbin,"z",1),flag.compress = T)
  X = as.data.frame(db.extract(dbin,db.getname(dbin,"f",NA),flag.compress = T))
  reslm = lm(Y~0+.,data=X)
  for(i in 1:ncol(X))
  {
    nm_estim <- paste(radix, paste0("Z",i), "estim", sep = ".")
    dbout = db.replace(dbout, name = nm_estim, tab = rep(reslm$coefficients[i],dbout$nactive),flag.sel = T)
  }
  
  return(dbout)
}

# Regression linéaire par unité administrative (nombre de biens > N)
fn_estim_LM2 = function(dbin,dbout,radix = "kfold_LM_admin")
{
  Y = db.extract(dbin,db.getname(dbin,"z",1),flag.compress = T)
  X = as.data.frame(db.extract(dbin,db.getname(dbin,"f",NA),flag.compress = T))
  cityTrain = db.extract(dbin,"COM_INSEE",flag.compress = T)
  cityTest = db.extract(dbout,"COM_INSEE",flag.compress = T)
  reslmGlobal = lm(Y~0+.,data=X)
  res = matrix(NA,nrow=dbout$nactive,ncol=ncol(X))
  for(i in 1:ncol(X))
  {
    res[,i] = reslmGlobal$coefficients[i]
  }
  for(i in unique(cityTrain))
  {
    ind = cityTrain==i
    if(sum(ind)>N)
    {
      x = as.data.frame(X[ind,])
      y = Y[ind]
      reslm = lm(y~0+.,data=x)
      indTest = cityTest == i
      for(j in 1:ncol(x))
      {
        res[indTest,j] = reslm$coefficients[j]
      }
    }
  }
  
  for(i in 1:ncol(X))
  {
    nm_estim <- paste(radix, paste0("Z",i), "estim", sep = ".")
    dbout = db.replace(dbout, name = nm_estim, tab = res[,i],flag.sel = T)
  }
  
  return(dbout)
}
