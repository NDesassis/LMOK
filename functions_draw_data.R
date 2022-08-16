# Des fonctions de dessin pour les departements et les préfectures
drawDept <- function(departements, idx = NA, flag.add = FALSE, flag.name=TRUE, col = "red"){
  if(is.na(idx)) {idx <- (1:length(departements))}
  for (i in seq_along(idx)){
    plot(departements[[i]]$polygon, add = flag.add)
    if(flag.name){
      points(departements[[i]]$city_xy[1], departements[[i]]$city_xy[2], pch = 19, col = col)
      text(
        departements[[i]]$city_xy[1], 
        departements[[i]]$city_xy[2]+5000, 
        departements[[i]]$city, col = col)
    }
  }
  NULL
}

boundingboxDept <- function(departements, idx = NA){
  bb <- matrix(NaN, 2,2)
  if(is.na(idx)) {idx <- (1:length(departements))}
  for (i in seq_along(idx)){
    b <- departements[[i]]$polygon$limits
    if (i == 1) {bb <- b} 
    else {
      bb[1,1] <- min(bb[1,1], b[1,1])
      bb[1,2] <- min(bb[1,2], b[1,2])
      bb[2,1] <- max(bb[2,1], b[2,1])
      bb[2,2] <- max(bb[2,2], b[2,2])
      }
    }
  colnames(bb) <- c("X","Y")
  rownames(bb) <- c("Min","Max")
  bb
}
