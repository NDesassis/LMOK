# functions to draw polygon and the main cities of the french departments

#' draw the limits of french departments
#' 
#'@param department The list of departments. For each department the following attributes are defined
#' - polygon A polygon-class RGeostats structure containing the limits of the area
#' - city A string containing the name of the main city (i.e. the french *préfecture*)
#' - city_xy The real vector of the coordinates of the main city
#'@param idx The integer vector containing the indices in the list of the departments to be drawn
#' If idx == NA, all the elements of the list are taken into account
#'@param flag.add A Boolean variable. 
#'If flag.add == TRUE, the drawing is performed on the active figure. 
#'If flag.add == FALSE, a new figure is open.
#'@param flag.name A Boolean variable.
#'Only if flag.name == TRUE, the main city name 
#'@param col
#'@return NULL value
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

#' compute the bounding box of a list of french departments
#' 
#'@param department The list of departments. For each department the following attributes are defined
#' - polygon A polygon-class RGeostats structure containing the limits of the area
#' - city A string containing the name of the main city (i.e. the french *préfecture*)
#' - city_xy The real vector of the coordinates of the main city
#'@param idx The integer vector containing the indices in the list of the departments to be drawn
#' If idx == NA, all the elements of the list are taken into account
#'@return A 2x2 table containing the coordinates of the bounding box
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
