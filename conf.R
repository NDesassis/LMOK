source("filter_examples.R")
source("functions_prepareDB.R")
# General parameters
set.seed(seed=125875)  # Initialisation of the random generator (to reproduce simulations)
nsim = 1000            # Number of simulations to compute the conditional expectation
n.fold = 10            # Number of folds to compute the cross-validation


# zone à traiter
area_nm   <- "Rhone" # Selection de la zone dans "BPL", "Creuse", "Oise", "Rhone"
nm_input_raw   <- paste0("../Data/CCR_data_", area_nm, ".Rdata")
nm_input_RData <- paste0("../Data/CCR_", area_nm, ".Rdata")
nm_input_mesh  <- paste0("../Data/CCR_", area_nm, "_meshCoords.ascii")
nm_output_mesh <- nm_input_mesh
# sauvegarde du modèle
nm_output_model<- paste0("../Data/CCR_", area_nm, "_model_SPDE.Rdata")

# Parameters of the linear model
nu      <- 1          # regularity coefficient of the Matern covariance
nm_var  <- "S_cap"    # Name of the variable to be predicted and observed (surface)
nm_fac  <- "NB_cap"   # Name of the factor (number of rooms)
nm_type <- "C1NATURE" # Name of the type (1 = House / 2 = Flat)

# Definition of the meshing of the domain used for SPDE
nodes   <- c(150,150) # Number of nodes in the main directions
margin  <- 0.2        # Proportion to extend the bouding box of the domain
param   <- list(margin = margin, nodes = nodes, nu = nu)

# The mesh is computed by the script *makeMeshing.Rmd*
coords = read.csv(file = nm_input_mesh)
pts = db.create(x1=coords[,2],x2=coords[,3])
#mesh = meshing(pts,triswitch="nqQ")
mesh = NA
#plot(mesh,lwd=0.1)
#plot(poly,add=T)

# Choose model type among the next combinations
type   = "both" # "house", "flat", "both"
model  = "M3"    # "Mono", "M1", "M2", "M3" 
nm_sel = paste0("sel_",type)
nm_model = paste0(type,"_",model)

# number of data under which global regression is used (for regression by city)
N_adm = 500

source(file = paste0("./param_", model,".R"))

# Options for the optimizer (maximization of the likelihood)
optim_maxit = 1
optim_method = "Nelder-Mead"


#filters: a list which contains all the filters to apply. The filters are some functions
#They can be written in the file "filter_examples.R" which already contains some examples.

filter_list = list(lyon_suburb)
# filter_list = list() # no filter
