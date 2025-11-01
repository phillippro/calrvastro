source('mcmc_func.R')
args <- commandArgs(trailingOnly=TRUE)
dt <- as.numeric(args[1])
dt2 <- as.numeric(args[2])
###HIP 24186
obs <- c(ra=077.89672072,dec=-45.00448692,parallax=255.66,pmra=6505.08,pmdec=-5730.84,radial_velocity=245.2)
obs1 <- obs.lin.prop(obs,c(0,dt)*365.25)
#obs1 <- vec.lin.prop(obs,c(0,dt)*365.25)
#dra <- diff(obs1[,'ra'])*3.6e6*cos(obs1[1,'dec']/180*pi)-obs1[1,'pmra']*dt
#dra <- diff(obs1[,'ra'])*3.6e6*cos(mean(obs1[,'dec'])/180*pi)-obs1[1,'pmra']*dt
dra <- diff(obs1[,'ra'])*3.6e6*cos(obs1[2,'dec']/180*pi)-obs1[1,'pmra']*dt
ddec <- diff(obs1[,'dec'])*3.6e6-obs1[1,'pmdec']*dt
ang <- sqrt(dra^2+ddec^2)
print(obs)
cat('obs1[1,]-obs=',obs1[1,]-obs,'\n')
cat('theta=',ang,'mas for ',dt,'years!\n')
#obs1 <- obs.lin.prop(obs,c(0,dt2)*365.25)
obs1 <- obs.lin.prop(obs,c(-dt2/2,dt2/2)*365.25)
dpmra <- diff(obs1[,'pmra'])
dpmdec <- diff(obs1[,'pmdec'])
cat('pmra/pmdec=',obs['pmra']/obs['pmra'],'\n')
cat('dpmra/dpmdec=',dpmra/dpmdec,'\n')
mudot <- diff(sqrt(obs1[,'pmra']^2+obs1[,'pmdec']^2))/dt2
mu <- sqrt(obs['pmra']^2+obs['pmdec']^2)
#mudot <- sqrt(diff(obs1[,'pmra'])^2+diff(obs1[,'pmdec'])^2)/dt2
cat('mudot=',mudot,'mas/yr^2\n')
mur <- obs['radial_velocity']/4.74047*obs['parallax']
cat('mudot=-2mu*mur=',-2*mu*mur/3.6e6*(pi/180),'mas/yr^2\n')


