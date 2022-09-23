



# Paramètres du modèle

nu      <- 1        # Coefficient de la covariance de Matern
nodes   <- c(100,100) # Valeur par default pour le meshing
margin  <- 0.2      # Proportion pour l'extension de la zone à mailler

set.seed(seed=125875)

nsim = 1000 # Number of simulations to compute CE
np = 100    # Number of samples in the training data base
n.fold = 10 # Number of folds for the cross-validation
nm_var = "S_cap"
nm_fac = "NB_cap"
nm_type = "C1NATURE"

coords = read.csv("../Data/meshCoords.ascii")
pts = db.create(x1=coords[,2],x2=coords[,3])
mesh = meshing(pts,triswitch="nqQ")
#mesh = NA
#plot(mesh,lwd=0.1)
#plot(poly,add=T)


# Choose model type among the next combinations
type = "house" #"flat" "house" "mono_house" "mono_flat" or "both" (only valid for M3)
model = "M1" #Mono M1 M2 MonoBivar

nm_sel = paste0("sel_",type)
nm_model = paste0(type,"_",model)

source(paste0(model,".R"))
