#' compute the adpated meshing for SDPE modeling (LMOK)
#' 
#' The SPDE approach is used to compute the estimations of the latent fields.
#' 
#'@param poly a polygon-class structure containing the limits of the domain.
#'@param dat The db-class structure containing the data file
#'@param nodes The number of nodes to discretize the domain
#'@param margin The number of the nodes defining the margins
#'@param lambda The Poisson density controlling the number of vertices
#'@param N The average number of vertices of the mesh (lambda or N should be given)
#'@return The computed mesh is return
#'
compute_mesh <- function(poly, dat, nodes, margin, lambda=NA, N=NA){
  # limits of the mesh from polygon and data set
  grd=db.grid.init(dat,nodes = c(50,50),margin=c(10,10))
  grd=db.polygon(grd, poly)
  grd = db.locate(grd,"Polygon","z")
  # limits are inflated using morphological operators
  grd = morpho(grd,0.5,1.5,oper = "close")
  grd = morpho(grd,0.5,1.5,oper = "dilate",niter = 10)
  grd = morpho(grd,0.5,1.5,oper = "erode")
  # coarse grid
  grdCoarse = db.create(flag.grid = TRUE,
                        nx = grd$nx/2,
                        dx = grd$dx*2,
                        x0=grd$x0+grd$dx/2)
  grdCoarse = migrate(dbin = grd, dbout = grdCoarse, radix="")
  grdCoarse = db.rename(grdCoarse,grdCoarse$natt,"morpho")
  grdCoarse = db.stat.grid(dat,grdCoarse,"x1",fun="num")
  grdCoarse = db.polygon(grdCoarse,polygon = poly)
  coords = NULL
  if(is.na(lambda)){
    lambda <- N / sum(grdCoarse[,"Stats.rank"], na.rm = TRUE)
    print(paste0("Lambda = ", lambda))
  }
  n = rpois(grdCoarse$nech,lambda * grdCoarse[,"Stats.rank"])
  n[grdCoarse[,"Polygon"]==1] = apply(cbind(n[grdCoarse[,"Polygon"]==1],2),1,max)
  n[grdCoarse[,"Polygon"]==0 & grdCoarse[,"morpho"]==1] = 1
  for(i in 1:grdCoarse$nech) {
    x = runif(n[i]) * grdCoarse$dx[1] + grdCoarse[i,"x1"]
    y = runif(n[i]) * grdCoarse$dx[2] + grdCoarse[i,"x2"]
    coords = rbind(coords,cbind(x,y))
  }
  pts  <- db.create(x1=coords[,1],x2=coords[,2])
  meshing(pts,triswitch="nqQ")
}
