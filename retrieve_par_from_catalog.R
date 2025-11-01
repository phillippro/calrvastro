args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    star1 <- star <- args[1]
    if(length(args)>1) star1 <- args[2]
}else{
    star <- '55 Cnc'
    star1 <- 'HD75732'
}
source('mcmc_func.R')
tab <- read.csv('exoplanet_eu220320.csv')
ii <- which(tab[,'star_name']==star)
msini <- tab[ii,'mass_sini']
mp <- tab[ii,'mass']#mjup
jj <- which(is.na(msini))
inc <- tab[ii,'inclination']/180*pi
if(length(jj)>0) msini[jj] <- mp[jj]*sin(inc[jj])
Rp <- tab[ii,'radius']
e <- tab[ii,'eccentricity']
p <- tab[ii,'orbital_period']
omega <- tab[ii,'omega']/180*pi
lambda <- tab[ii,'lambda_angle']
Mstar <- 1.015
eMstar <- 0.051
K <- msini2K(msini,p,e,Mstar)
ns <- c('per','K','e','omega','Mo')
out <- c()
for(j in 1:length(ii)){
    out <- c(out,paste0('per',j,' ',log(p[j])))
    out <- c(out,paste0('K',j,' ',K[j]))
    out <- c(out,paste0('e',j,' ',e[j]))
    if(!is.na(omega[j])){
        og <- omega[j]
    }else{
        og <- 0
    }
    out <- c(out,paste0('omega',j,' ',og))
    if(!is.na(inc[j])){
        I <- inc[j]
    }else{
        I <- pi/2
    }
    out <- c(out,paste0('Inc',j,' ',I))
}
fout <- paste0('pars/',star1,'.par')
cat(fout,'\n')
write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
