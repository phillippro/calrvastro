###This routine is used to identify and constrain signals the data using MCMC method for a given noise model

##step 1:set up
Nsig <- 0
data1 <- data0 <- tab[,1:3]#data0 is the raw data; data 1 is data0-signals
assign.par <- FALSE

###step 2: pure noise model; raw data; cold chain
data <- data1
source('setting_up.R')
Np <- 0
source('initial_condition.R')
##cold chain
tem <- 1
out <- AMH(round(Niter/2/Ncores),Niter/Ncores)$out
##save data
output.all[[noise]][[1]] <- out
ll0 <- max(out[,ncol(out)])

###step 3:noise+1planet model; raw data; hot chain
#data <- data1
#source('setting_up.R')
out.all <- c()
Np <- 1
source('initial_condition.R')
source('hot_chain.R')
out.all <- rbind(out.all,out)
startvalue <- out[which.max(out[,Npar+1]),1:Npar]

###step 4:noise+1planet model; raw data; cold chain
tem <- 1
out <- AMH(round(Niter/2/Ncores),Niter/Ncores)$out
par.opt <- out[which.max(out[,ncol(out)]),1:Npar]
out.all <- rbind(out.all,out)
output.all[[noise]][[2]] <- out.all
source('plot_BFP.R')
ll0 <-max(lls)

###step 5-?:noise+Nplanet model
if(max(logbf)>5 | mcmc.type=='test'){
    for(k2 in 2:Nmax){
        out.all <- list()
        cat('\nIdentify and constrain',j,'signal!\n')

        ##noise+1planet; raw-(j-1)planet; hot chain
        data1[,2] <- data0[,2]-RV.kepler(pars.kep=par.opt,kep.only=TRUE,Np.kep=1)$rv[[1]]
        data <- data1
        source('setting_up.R')
        Np <- 1
        source('initial_condition.R')
        source('hot_chain.R')
        out.all[[1]] <- out.tmp
        par.1sig <- as.numeric(out[which.max(out[,Npar+1]),1:Nkeppar])

        ##cold chain
        data <- data0
        source('setting_up.R')
        Np <- k2
        source('initial_condition.R')
        startvalue <- assign.names(c(par.1sig,as.numeric(par.opt)),Np=Np,p=Nar.opt,q=Nma.opt)
#        startvalue[-(1:Nkeppar)] <- par.opt
        tem <- 1
        out <- AMH(round(Niter/2/Ncores),Niter/Ncores)$out
        par.opt <- out[which.max(out[,Npar+1]),1:Npar]
        out.all[[2]] <- out
        output.all[[noise]][[j+1]] <- out.all
        source('plot_BFP.R')
        ll0 <-max(lls)
    }
}
##
acceptance  <- (1-mean(duplicated(mcmc.out)))*100
out <- amh$out
cat('acceptance percentage:',acceptance,'\n')

