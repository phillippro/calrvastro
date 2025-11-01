args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    star1 <- star <- args[1]
    type <- args[2]
    if(length(args)>2) star1 <- args[3]
}else{
    star <- '55 Cnc'
    type <- 'RV'
    star1 <- 'HD75732'
}
source('mcmc_func.R')
tab <- read.csv('exoplanet_eu220320.csv')
ii <- which(tab[,'star_name']==star)
msini <- tab[ii,'mass_sini']
mp <- tab[ii,'mass']#mjup
jj <- which(is.na(msini))
inc <- tab[ii,'inclination']/365.25
if(length(jj)>0) msini[jj] <- mp[jj]*sin(inc[jj])
Rp <- tab[ii,'radius']
e <- tab[ii,'eccentricity']
p <- tab[ii,'orbital_period']
omega <- tab[ii,'omega']
lambda <- tab[ii,'lambda_angle']
Mstar <- 1.015
eMstar <- 0.051
K <- msini2K(msini,p,e,Mstar)
ns <- c('per','K','e','omega','Mo')
out <- c()
for(j in 1:length(ii)){
    out <- c(out,paste0('per',j,' ',p[j]))
    out <- c(out,paste0('K',j,' ',K[j]))
    out <- c(out,paste0('e',j,' ',e[j]))
    out <- c(out,paste0('omega',j,' ',omega[j]))
    out <- c(out,paste0('Mo',j,' ',0))
}
if(type=='RV'){
    fout <- paste0('~/Documents/projects/dwarfs/red3/pars/',star1,'.par')
}else if(type=='astro'){
    fout <- paste0('~/Documents/projects/dwarfs/calrvastro/pars/',star1,'.par')
}
cat(fout,'\n')
write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
