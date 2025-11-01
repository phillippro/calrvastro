Nw <- 1#fit models to multiple wavelength data sets simultaneously
prior.type <- 'mt'
rjd <- TRUE

ids <- ids[1:Nw]

par.global <- c('trend')
#par.global <- ''
peta <- 1
#offset <- 'small'
CHIB <- FALSE
Nsubsamp <- 10
offset <- 'large'
gp.single <- TRUE
Njitter <- 1
Nb <- 1
detrend <- FALSE
detrend.type <- 'poly'
detrend.plot <- FALSE
##########################################
##############part I: global parameters
##########################################
##############check if the parameters can be passed from previous runs
time.unit <- 365.24
if(any(grepl('U1',ids))){
    time.unit <- 1
}
frac.burn <- 0.5
frac.adapt <- 1
index.more <- FALSE
dRV.mode <- 'force'#center; seq
#independence <- TRUE#assuming that the Keplerian signal, noise components in different aperture data sets are independent
independence <- FALSE
if(!exists('quantify')){
    quantify <- FALSE
}
#if(Np==0) quantify <- TRUE
workflow <- TRUE
tree <- FALSE
tree.burn <- FALSE
#chain.type <- 'parallel'#'section', 'adapt', 'parallel','normal'
sim.data <- FALSE
data.mode <- 'single'#'combined'#combined different data sets for one target
scaled <- TRUE
Nupdate <-  1#adaptive frequency
Nburn = round(Niter/Ncores*(frac.burn))#burning fraction
n0 = max(round(Niter/Ncores*(1-frac.adapt)+2),10)#where the adaptive cov. used#at least start from the second interation
cat('n0=',n0,'\n')
#n0 = Niter
eps = 1e-8##optional covariance parameter to avoid negative determinate
Ntry = 1000#the times limit to try
if(exists('prior.type0')){
    prior.type <- prior.type0
    rm(prior.type0)
}
verbose <- FALSE
#verbose <- TRUE
Nverbose <- round(Niter/10)#frequency of showing
#the parameters are stored in the following list object
par.data <- list()
par.data <- lapply(1:Nw,function(i) par.data[[ids[i]]] <- list())
names(par.data) <- ids
