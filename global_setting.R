chain.type <- 'optimal'# or 'multiple'
sig.type <- 'posterior'# likelihood
noMCMC <- TRUE
Prot <- NA
tau <- NA
#GPtype <- 'MA+GP'
MPonly <- FALSE
alpha.dP <- 1e-6
alpha.dTc <- 1e-4
period.par <- 'logP'
mode <- 10
Npoly <- floor(mode/10)
Npoly.sub <- mode%%10
prior.type <- 'mt'
sr <- 1
fl <- 1
fp <- 1
eta <- 1
CHIB <- FALSE
if(target=='HD70642'){
    par.global <- c('trend')
}else{
    par.global <- c('')
}
time.unit <- 365.25#for trend m/s/year

###global parameters
Nact.max <- 2#maximum activity signals investigated
xi <- 10#remove outlier
####hyper parameters for prior
s0 <- 1#hyper par of s prior
K0 <- 1#hyper par of K prior
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
####other parameters
Nboot <- 1e1
noise.types <- c('white','MA','AR')
#noise.types <- c('white','MA')
#noise.types <- c('white','MA','AR')
gp.par <- c(NA,log(50),log(100))#amplitude,logProt,logtau
Nkeppar <- 5
eps <- 1e-6
if(prior.type=='e0'){
    Nkeppar <- 3
}
Ntype <- length(noise.types)
pq <- 1
Nsamp <- 2#sampling of BFP
#refine <- TRUE
refine <- FALSE
eTc <- Tc <- NULL
sampling <- 'combined'
ind.proxy <- 0#
Ptransit <- ePransit <- NA
#transit.target<- c('HIP47103','HD30219','HD86226','HD86226PFS','EPIC249223471')
trans <- read.table('toi_info.txt',header=TRUE)
transit.target<- trans[,1]
#ind <- which(target==transit.target)
ind <- c()
###if use transiting planet parameters
#Ntr0 <- Ntr <- length(ind)
###if not using transiting planet parameters
Ntr0 <- Ntr <- 0
if(Ntr==0){
   ind.transit <- 0
}else{
   ind.transit <- 1:Ntr
   Ptransit <- trans[ind,'P']
   ePtransit <- trans[ind,'eP']
   Tc <- trans[ind,'Tc']+2457000
   eTc <- trans[ind,'eTc']
   basis <- 'linear2'
   basis2 <- 'linear1'
}
noise.only <- FALSE

###model comparison parameters
Nar.max <- 0
Nma.max <- 5

####mcmc parameters
Nbin.per <- 1
nbin.per <- 1
