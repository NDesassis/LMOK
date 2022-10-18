good_cedantes = function(dat)
{
  !dat[,"CED_CODE_ANON"]%in%c(5,9,16,14,22,23,28,35)
}

lyon_suburb = function(dat)
{
  #Centroid communes
  centr = aggregate(dat[][c("X_TEMP","Y_TEMP")],by = list(dat[]$COM_INSEE),mean)
  
  #3ème arondissement de Lyon d'après :
  #https://www.galichon.com/codesgeo/insee.php?commune=lyon&comm=1
  lyon_insee  = 69383
  coordslyon = as.numeric(centr[centr$Group.1==lyon_insee,c("X_TEMP","Y_TEMP")])
  distlyon = sqrt((centr$X_TEMP-coordslyon[1])^2 + (centr$Y_TEMP-coordslyon[2])^2)
  suburbs = centr$Group.1[distlyon < 10000]
  dat[]$COM_INSEE%in%suburbs
}