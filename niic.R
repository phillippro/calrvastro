###input parameters
nrank <- nbeta <- 4
Nchain <- Ncores
#betas <- c(0.001, 0.01, 0.1, 1)
betas <- 10^(seq(-7,0,length.out=Nchain))
tune.ladder <- 1
NstopTuneLadder <- 3e5
scale.tune.ladder <- 0.2
zero.stretch <- 0.1
Niter.stack <- 5e3#total iterations in a stack of batches
Niter.batch.base <- 20# number of iteration in a batch before judging should a swap be taken
Niter.batch.rand <- 5# random +/- range of n-iter-a-batch-base (must < n-iter-a-batch-base)
Nswap <- 1
Swapmode <- 0
Nstoptune <- 4e5#stop tuning proposals
Nbegintune <- 3e5
Niter.tune <- 1e3#tuning sigma_prop
ar.ok.lower <- 0.1
ar.ok.upper <- 0.45
ar.best <- 0.23
ar.accept.diff <- 0.1
sigma.scale.half.ratio <- 20

sigma.scale.min <- 5e-9
sigma.scale.max <- 0.2
sigma.jumpin.ratio <- 300
i.save.begin <- 1e5

init.rand.seed <- 20
init.gp.ratio <- 0.1# ratio of the init gaussian_proposal_sigma to the allowed prior range

a.min <- -20
a.max <- 20
b.min <- -50
b.max <- 50
d.min <- 0.01
d.max <- 30

###Alloc a N_beta*N_param transit array for mpi_init and mpi_swap
###S(Nrank x Nparam); step size matrix
root.rank <- 1#another name: root chain
sigma.RanksParm.root <- transit.BetaParm.root <- array(NA,dim=c(nbeta,Npar))
logpost <- rep(NA,nbeta)
###initiate
#mcmc <- foreach(kk=1:Ncores,.errorhandling = 'pass') %dopar% {
##Generate a initial random N_beta*N_parm parm array
for(j in 1:Nchain){
    source('initial_condition.R',local=TRUE)
    logpost[j] <- posterior(chain[1,],tem=betas[kk],bases=bases,RVonly=RVonly)
#    c(logpost=logpost,beta=betas[kk])
    sigma.gaussian.prop <- rep(NA,Npar)
}
