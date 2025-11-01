###In this new file, I have chagned the calcuation memthod of mu1 and correct a bug in the expression of Vr.ma;
####add celerite SHO kernel
library(fields)
library(MASS)
library(foreach)
library(doMC)
library(Matrix)
library(doParallel)
library(parallel)
library(e1071)
library(kernlab)
library(glasso)
library(JPEN)
library(matrixcalc)
library(mvtnorm)
library(MASS)
library(rootSolve)
Cauyr <- 173.144632674*365.25
Cauday <- 173.144632674
auyr2kms <- 4.74047
SRS <- 1.97412574336e-8#au
CMPS <- 299792458.0#speed of light in m/s
s2j <- 1047.348644#Solar mass to Jupiter mass
pc2au <- 206264.81
Me2s <- 3.003e-6#Earth mass in solar unit
Mj2s <- 1/1048
Me2j <- Me2s/Mj2s
####https://link.springer.com/article/10.1007/s10509-022-04066-1#data-availability
###http://adsabs.harvard.edu/abs/2013ApJS..208....9P
###http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
#m22 <- read.table('../data/combined/mamajek22.txt',header=TRUE,check.names=TRUE)
m22 <- read.table('~/Documents/projects/dwarfs/data/combined/mamajek22.txt',header=TRUE,check.names=TRUE)
m22[m22=='...'|m22=='....'|m22=='.....'|m22=='......'] <- NA
ms <- as.numeric(m22[,'Msun'])
mg <- as.numeric(gsub(':','',m22[,'M_G']))
ind <- which(!is.na(ms) & !is.na(mg))
mrl.m22 <- approxfun(ms[ind],mg[ind])
mlow.m22 <- min(ms[ind])
mup.m22 <- max(ms[ind])

CalibrateJitter <- function(p,a=1,b=2.16){
###this is to calibrate jitter
    p0 <- p
    p['ra_error'] <- sqrt((a*p0['ra_error'])^2+b^2)
    p['dec_error'] <- sqrt((a*p0['dec_error'])^2+b^2)
    p['parallax_error'] <- sqrt((a*p0['parallax_error'])^2+b^2)
    p['pmra_error'] <- sqrt((a*p0['pmra_error'])^2+b^2)
    p['pmdec_error'] <- sqrt((a*p0['pmdec_error'])^2+b^2)
    p['ra_dec_cov'] <- p['ra_dec_cov']/p0['ra_error']/p0['dec_error']*p['ra_error']*p['dec_error']
    p['ra_parallax_cov'] <- p['ra_parallax_cov']/p0['ra_error']/p0['parallax_error']*p['ra_error']*p['parallax_error']
    p['ra_pmra_cov'] <- p['ra_pmra_cov']/p0['ra_error']/p0['pmra_error']*p['ra_error']*p['pmra_error']
    p['ra_pmdec_cov'] <- p['ra_pmdec_cov']/p0['ra_error']/p0['pmdec_error']*p['ra_error']*p['pmdec_error']
    p['dec_parallax_cov'] <- p['dec_parallax_cov']/p0['dec_error']/p0['parallax_error']*p['dec_error']*p['parallax_error']
    p['dec_pmra_cov'] <- p['dec_pmra_cov']/p0['dec_error']/p0['pmra_error']*p['dec_error']*p['pmra_error']
    p['dec_pmdec_cov'] <- p['dec_pmdec_cov']/p0['dec_error']/p0['pmdec_error']*p['dec_error']*p['pmdec_error']
    p['parallax_pmra_cov'] <- p['parallax_pmra_cov']/p0['parallax_error']/p0['pmra_error']*p['parallax_error']*p['pmra_error']
    p['parallax_pmdec_cov'] <- p['parallax_pmdec_cov']/p0['parallax_error']/p0['pmdec_error']*p['parallax_error']*p['pmdec_error']
    p['pmra_pmdec_cov'] <- p['pmra_pmdec_cov']/p0['pmra_error']/p0['pmdec_error']*p['pmra_error']*p['pmdec_error']
    p
}
###ref:https://www.aanda.org/articles/aa/pdf/2023/04/aa44144-22.pdf
###Lindegren et al. 2018
GaiaCalibrate <- function(astrometry,par,dt){
    ex <- par[1]
    ey <- par[2]
    ez <- par[3]
    wx <- par[4]
    wy <- par[5]
    wz <- par[6]
    dplx <- par[7]
    eX <- ex+wx*dt
    eY <- ey+wy*dt
    eZ <- ez+wz*dt
    alpha <- astrometry['ra']/180*pi
    ddec <- as.numeric(- sin(alpha)*eX + cos(alpha)*eY)/3.6e6#deg
    astrometry['dec'] <- astrometry['dec'] + ddec

    delta <- astrometry['dec']/180*pi
    dra <- as.numeric((cos(alpha)*sin(delta)*eX+sin(alpha)*sin(delta)*eY-cos(delta)*eZ)/3.6e6)#deg
    astrometry['ra'] <- astrometry['ra'] + dra/cos(delta)

    dpmra <- as.numeric(cos(alpha)*sin(delta)*wx+sin(alpha)*sin(delta)*wy-cos(delta)*wz)
    astrometry['pmra'] <- astrometry['pmra'] + dpmra

    dpmdec <- as.numeric(- sin(alpha)*wx + cos(alpha)*wy)
    astrometry['pmdec'] <- astrometry['pmdec'] + dpmdec

    astrometry['parallax'] <- astrometry['parallax'] + dplx

    cat('dra=',dra*3.6e6,'mas;ddec=',ddec*3.6e6,'mas;dpmra=',dpmra,'mas/yr;dpmdec=',dpmdec,'\n')
    astrometry
}
###ref:Feng+2024
Calibrate2GDR3 <- function(astrometry,par,dt,gamma=-1){
###gamma=-1 (subtract bias from the given catalog); gamma=1 (add bias to the given catalog)
    ex <- par[1]
    ey <- par[2]
    ez <- par[3]
    wx <- par[4]
    wy <- par[5]
    wz <- par[6]
    dplx <- par[7]
    eX <- ex+wx*dt
    eY <- ey+wy*dt
    eZ <- ez+wz*dt

    alpha <- as.numeric(astrometry['ra']/180*pi)#rad
    ddec <- as.numeric(- sin(alpha)*eX + cos(alpha)*eY)/3.6e6#deg
    astrometry['dec'] <- astrometry['dec'] + gamma*ddec

    delta <- as.numeric(astrometry['dec']/180*pi)#rad
    dra.star <- as.numeric((cos(alpha)*sin(delta)*eX+sin(alpha)*sin(delta)*eY-cos(delta)*eZ)/3.6e6)#deg
    dra <- dra.star/cos(delta)
    astrometry['ra'] <- astrometry['ra'] +  gamma*dra

    dpmra <- as.numeric(cos(alpha)*sin(delta)*wx+sin(alpha)*sin(delta)*wy-cos(delta)*wz)
    astrometry['pmra'] <- astrometry['pmra'] + gamma*dpmra

    dpmdec <- as.numeric(- sin(alpha)*wx + cos(alpha)*wy)
    astrometry['pmdec'] <- astrometry['pmdec'] + gamma*dpmdec

    astrometry['parallax'] <- astrometry['parallax'] + gamma*dplx

    cat('dra.star=',dra.star*3.6e6,'mas;dra=',dra*3.6e6,'mas;ddec=',ddec*3.6e6,'mas;dplx=',dplx,'mas;dpmra=',dpmra,'mas/yr;dpmdec=',dpmdec,'mas/yr\n')
    astrometry
}

convert.angle <- function(angle,mid=pi,max=2*pi){
    angle1 <- angle <- angle%%max
    angle1[angle1>mid] <- angle1[angle1>mid]-max
    if(sd(angle1)<sd(angle)) angle <- angle1
    angle
}

###ref:https://www.aanda.org/articles/aa/pdf/2023/04/aa44144-22.pdf
###Lindegren et al. 2018
frame.rotation <- function(astrometry,rot,dt){
    wx <- rot[1]
    wy <- rot[2]
    wz <- rot[3]
    alpha <- astrometry['ra']/180*pi
    ddec <- as.numeric(- sin(alpha)*wx + cos(alpha)*wy)*dt/3.6e6#deg
    astrometry['dec'] <- astrometry['dec'] + ddec

    delta <- astrometry['dec']/180*pi
    dra <- as.numeric((cos(alpha)*sin(delta)*wx+sin(alpha)*sin(delta)*wy-cos(delta)*wz)*dt/3.6e6)#deg
    astrometry['ra'] <- astrometry['ra'] + dra

    dpmra <- as.numeric(cos(alpha)*sin(delta)*wx+sin(alpha)*sin(delta)*wy-cos(delta)*wz)
    astrometry['pmra'] <- astrometry['pmra'] + dpmra

    dpmdec <- as.numeric(- sin(alpha)*wx + cos(alpha)*wy)
    astrometry['pmdec'] <- astrometry['pmdec'] + dpmdec

    cat('dra=',dra*3.6e6,'mas;ddec=',ddec*3.6e6,'mas;dpmra=',dpmra,'mas/yr;dpmdec=',dpmdec,'\n')
    astrometry
}

cor2cov.full <- function(p){
  n1 <- c('ra_dec_corr','ra_parallax_corr','ra_pmra_corr','ra_pmdec_corr')
  n2 <- c('dec_parallax_corr','dec_pmra_corr','dec_pmdec_corr')
  n3 <- c('parallax_pmra_corr','parallax_pmdec_corr')
  n4 <- 'pmra_pmdec_corr'
  n5 <- c('ra_error','dec_error','parallax_error','pmra_error','pmdec_error')
  n <- c(n1,n2,n3,n4)
  cor1 <- cov1 <- rep(0,10)
  names(cor1) <- names(cov1) <- n
  cor1[n] <- p[n]
  cor0 <- array(0,dim=c(5,5))
  cor0[1,2:5] <- as.numeric(p[n1])
  cor0[2,3:5] <- as.numeric(p[n2])
  cor0[3,4:5] <- as.numeric(p[n3])
  cor0[4,5] <- as.numeric(p[n4])
  err0 <- as.numeric(p[n5])
###make it symmetric
  cor0 <- t(cor0) + cor0
  diag(cor0) <- 1
  cov0 <- cor2cov(cor0,err0)
  cov1[n1] <- cov0[1,2:5]
  cov1[n2] <- cov0[2,3:5]
  cov1[n3] <- cov0[3,4:5]
  cov1[n4] <- cov0[4,5]
  return(list(cor=cor0,cov=cov0,cov1=cov1,cor1=cor1))
}
cov2vec <- function(cov){
    n1 <- c('ra_dec_cov','ra_parallax_cov','ra_pmra_cov','ra_pmdec_cov')
    n2 <- c('dec_parallax_cov','dec_pmra_cov','dec_pmdec_cov')
    n3 <- c('parallax_pmra_cov','parallax_pmdec_cov')
    n4 <- 'pmra_pmdec_cov'
    n5 <- c('ra_error','dec_error','parallax_error','pmra_error','pmdec_error')
    vec <- c(cov[1,2:5],cov[2,3:5],cov[3,4:5],cov[4,5],sqrt(diag(cov)))
    names(vec) <- c(n1,n2,n3,n4,n5)
    return(vec)
}
vec2cov <- function(vec){
    n1 <- c('ra_dec_cov','ra_parallax_cov','ra_pmra_cov','ra_pmdec_cov')
    n2 <- c('dec_parallax_cov','dec_pmra_cov','dec_pmdec_cov')
    n3 <- c('parallax_pmra_cov','parallax_pmdec_cov')
    n4 <- 'pmra_pmdec_cov'
    n5 <- c('era','edec','eparallax','epmra','epmdec')
    cov <- diag(vec[n5]^2)
    cov[1,2:5] <- cov[2:5,1] <- as.numeric(vec[n1])
    cov[2,3:5] <- cov[3:5,2] <- as.numeric(vec[n2])
    cov[3,4:5] <- cov[4:5,3] <- as.numeric(vec[n3])
    cov[4,5] <- cov[5,4] <- as.numeric(vec[n4])
#    colnames(cov) <- rownames(cov) <- c('ra','dec','parallax','pmra','pmdec')
    return(cov)
}
cor2vec <- function(cor){
    n1 <- c('ra_dec_cor','ra_parallax_cor','ra_pmra_cor','ra_pmdec_cor')
    n2 <- c('dec_parallax_cor','dec_pmra_cor','dec_pmdec_cor')
    n3 <- c('parallax_pmra_cor','parallax_pmdec_cor')
    n4 <- 'pmra_pmdec_cor'
#    n5 <- c('ra_error','dec_error','parallax_error','pmra_error','pmdec_error')
    vec <- c(cor[1,2:5],cor[2,3:5],cor[3,4:5],cor[4,5])
    names(vec) <- c(n1,n2,n3,n4)
    return(vec)
}
cov2cor.full <- function(p){
  n1 <- c('ra_dec_cov','ra_parallax_cov','ra_pmra_cov','ra_pmdec_cov')
  n2 <- c('dec_parallax_cov','dec_pmra_cov','dec_pmdec_cov')
  n3 <- c('parallax_pmra_cov','parallax_pmdec_cov')
  n4 <- 'pmra_pmdec_cov'
  n5 <- c('ra_error','dec_error','parallax_error','pmra_error','pmdec_error')
  n <- c(n1,n2,n3,n4)
  cor1 <- cov1 <- rep(0,10)
  names(cor1) <- names(cov1) <- n
  cov1[n] <- p[n]
  cov0 <- array(0,dim=c(5,5))
  cov0[1,2:5] <- as.numeric(p[n1])
  cov0[2,3:5] <- as.numeric(p[n2])
  cov0[3,4:5] <- as.numeric(p[n3])
  cov0[4,5] <- as.numeric(p[n4])
  err0 <- as.numeric(p[n5])
###make it symmetric
  cov0 <- t(cov0) + cov0
  diag(cov0) <- err0^2
  cor0 <- cov2cor(cov0,err0)
  names(cor1) <- n
  cor1[n1] <- cor0[1,2:5]
  cor1[n2] <- cor0[2,3:5]
  cor1[n3] <- cor0[3,4:5]
  cor1[n4] <- cor0[4,5]

  return(list(cor=cor0,cov=cov0,cov1=cov1,cor1=cor1))
}

u2cov <- function(ut,err,u){
    cors <- covs <- array(NA,dim=c(nrow(ut),10))
    for(j in 1:nrow(ut)){
        c1 <- as.numeric(c(ut[j,1],rep(0,4)))
        c2 <- as.numeric(c(ut[j,2:3],rep(0,3)))
        c3 <- as.numeric(c(ut[j,4:6],rep(0,2)))
        c4 <- as.numeric(c(ut[j,7:10],rep(0,1)))
        c5 <- as.numeric(ut[j,11:15])
        U <- cbind(c1,c2,c3,c4,c5)
        cov.mat <- solve(t(U)%*%U)
        if(u[j]>1) cov.mat <- cov.mat*u[j]^2
        cor.mat <- cov.mat/outer(as.numeric(err[j,]),as.numeric(err[j,]),"*")
        covs[j,] <- as.numeric(c(cov.mat[1,2:5],cov.mat[2,3:5],cov.mat[3,4:5],cov.mat[4,5]))
        cors[j,] <- as.numeric(c(cor.mat[1,2:5],cor.mat[2,3:5],cor.mat[3,4:5],cor.mat[4,5]))
    }
    return(list(cov=covs,cor=cors))
}

cor2cov <- function(cor,var){
    sd1 <- replicate(length(var),var)
    sd2 <- t(sd1)
    sd1*cor*sd2
}
cov2cor <- function(cov,var){
    w1 <- replicate(length(var),1/var)
    w2 <- t(w1)
    w1*cov*w2
}

time_Doy2jd <- function(yr,doy){
####################################
## Convert integer Year and nubmer of days since the begining of the year to 2-part Julian Date
## ref: https://www.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox
##
## Input:
##   yr - Integer Year
##   doy - Number of days since the begining of the year
##
## Output:
##   jd - 2-part Julian Date
####################################
    jd <- time_Cal2JD(cbind(yr,1,0))
    cbind(jd[,1]+doy,jd[,2])
}

time_Yr2jd <- function(yr){
####################################
## Convert Year to Julian Day
## ref: https://www.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox
##
## Input:
##   yr -  Year
##
## Output:
##   jd - Julian Date
####################################
    iyr <-  floor(yr)
    jd0  <-  time_Cal2JD(cbind(iyr,1,1))
    days  <-  rowSums(time_Cal2JD(cbind(iyr+1,1,1)) - jd0)
    doy  <-  (yr-iyr)*days + 1
    time_Doy2jd(iyr,doy)
}

time_Cal2JD <- function(cal){
####################################
## Gregorian Calendar date to Julian Date.
## ref: https://uk.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox?focused=6513448&tab=function
## CAL2JD  Converts calendar date to Julian date using algorithm
## from "Practical Ephemeris Calculations" by Oliver Montenbruck
##   (Springer-Verlag, 1989). Uses astronomical year for B.C. dates
##   (2 BC = -1 yr). Non-vectorized version. See also DOY2JD, GPS2JD,
##   JD2CAL, JD2DOW, JD2DOY, JD2GPS, JD2YR, YR2JD.
##
## Input:
##   cal - Calendar date with day fraction
##
## Output:
##   JD - 2-part or 2D Julian Date
####################################
    y <- yr <- cal[,1]
    m <- mn <- cal[,2]
    dy <- cal[,3]
    ind <- which(mn<=2)
    y[ind] <- y[ind]-1
    m[ind] <- m[ind]+12
    date1 <- 4.5+31*(10+12*1582)# Last day of Julian calendar (1582.10.04 Noon)
    date2 <- 15.5+31*(10+12*1582)# First day of Gregorian calendar (1582.10.15 Noon)
    date <- dy+31*(mn+12*yr)
    ind1 <- which(date<=date1)
    ind2 <- which(date>=date2)
    b <- y
    b[ind1] <- -2
    b[ind2] <- trunc(y/400) - trunc(y/100)
    if(length(ind1)==0 & length(ind2)==0){
        cat('Dates between October 5 & 15, 1582 do not exist!\n')
    }
    ind1 <- which(y>0)
    ind2 <- which(y<0)
    jd <- y
    jd[ind1] <- trunc(365.25*y[ind1]) + trunc(30.6001*(m[ind1]+1)) + b[ind1] + 1720996.5 + dy[ind1]
    jd[ind2] <- trunc(365.25*y[ind2]-0.75) + trunc(30.6001*(m[ind2]+1)) + b[ind2] + 1720996.5 + dy[ind2]
#    return(cbind(DJM0,jd-DJM0))
    return(cbind(jd%/%1,jd%%1))
}

rd2ps <- function(dra,edra,ddec,eddec){
###convert dra, ddec to pa and separation
    sep <- sqrt(dra^2+ddec^2)
    esep <- sqrt((dra^2*edra^2+ddec^2*eddec^2)/(dra^2+ddec^2))
    pa <- (atan2(dra,ddec)%%(2*pi))*180/pi#or rho in deg
    epa <- sqrt((dra^2*eddec^2+ddec^2*edra^2)/(dra^2+ddec^2))/sep*180/pi#deg
    cbind(sep=sep,esep=esep,pa=pa,epa=epa)
}

ps2rd <- function(rho,erho,pa,epa){
##########
##purppose: rho, PA to Dra and Ddec
##input:
##rho: separation in mas
##erho: error of rho
##pa: projection angle in degree
##epa: error of pa
##########
    pa1 <- pa/180*pi
    epa1 <- epa/180*pi
    dec.mas <- rho*cos(pa1)
    ra.mas <- rho*sin(pa1)
    edec.mas <- sqrt((rho*sin(pa1)*epa1)^2+(erho*cos(pa1))^2)
    era.mas <- sqrt((rho*cos(pa1)*epa1)^2+(erho*sin(pa1))^2)
    cbind(dra=ra.mas,edra=era.mas,ddec=dec.mas,eddec=edec.mas)
}
tol1 <- 1e-16
yr2d <- 365.25
###modeL Of Keplerian Motions Of planets
###According to http://w.astro.berkeley.edu/~kclubb/pdf/RV_Derivation.pdf and Ford 2006
###The mass is from Berger, D. H.; et al. (2006). "First Results from the CHARA Array. IV. The Interferometric Radii of Low-Mass Stars".
###Ms=0.404
###assign names to the parameter vectors of a RV model
#Naming system: Keplerian parameters: {{per},{K},{e},{omega},{Mo}}; trend pars: {a, b}; noise pars (white noise: {s}, Gaussian process red noise: {sigma.white,l}, ARMA red noise:{{phi},alpha,{w},beta}>)
###Note that the x,y,z here should be in the heliocentric frame, and x axis point to the GC while y axis point to the rotation direction
xyz2bl.vec <- function(x,y,z){
    b <- atan2(z,sqrt(x^2+y^2))
    ind <- which(b>pi/2)
    if(length(ind)>0) b[ind] <- b[ind]-pi
    l <- atan2(y,x)%%(2*pi)
    return(cbind(b,l))
}
##kepler solver
kep.mv <- function(m,e){
#    tol = 1e-6
    tol = tol1
    m <- m%%(2*pi)
    E0 <- m+e*sin(m)+e^2*sin(m)*cos(m)+0.5*e^3*sin(m)*(3*cos(m)^2-1)#initial value
    Ntt <- 1e4
    for(k in 1:Ntt){
        t1 <- cos(E0)
        t2 <- -1+e*t1
        t3 <- sin(E0)
        t4 <- e*t3
        t5 <- -E0+t4+m
        t6 <- t5/(0.5*t5*t4/t2+t2)
        E1 = E0-t5/((1/2*t3 - 1/6*t1*t6)*e*t6+t2)
        if(all(abs(E1-E0)<tol)) break()
        if(k==Ntt) cat('Keplerian solver does not converge!\n')
        E0 <- E1
    }
    return(E1%%(2*pi))
}
deg2hdms <- function(RAdeg,DEdeg){
    val <- RAdeg/360*24#hour
    RAh <- floor(val)
    RAm <- floor((val%%1)*60)
    RAs <- (((val%%1)*60)%%1)*60
    sig <- sign(DEdeg)
    DEdeg <- abs(DEdeg)
    DEd <- sig*floor(DEdeg)
    DEm <- floor((DEdeg%%1)*60)
    DEs <- (((DEdeg%%1)*60)%%1)*60
    return(cbind(RAh,RAm,RAs,DEd,DEm,DEs))
}

hdms2deg <- function(rah,ram,ras,ded,dem,des){
    RA.deg <- (rah+ram/60+ras/3600)/24*360
    DEC.deg <- sign(ded)*(abs(ded)+dem/60+des/3600)
    RA.rad <- RA.deg/180*pi
    DEC.rad <- DEC.deg/180*pi
    cbind(RA.deg,DEC.deg,RA.rad,DEC.rad)
}
kep.mt <- function(m,e){
#    tol = 1e-8
    tol <- tol1
    E <- rep(NA,length(m))
    Ntt <- 1000
    for(j in 1:length(m)){
        E0 <- m[j]
        for(k in 1:Ntt){
            E1 = E0-(E0-e*sin(E0)-m[j])/(sqrt((1-e*cos(E0))^2-(E0-e*sin(E0)-m[j])*(e*sin(E0))))
            if(abs(E1-E0)<tol) break()
            if(j==Ntt) cat('Keplerian solver does not converge!\n')
            E0 <- E1
        }
        E[j] <- E1%%(2*pi)
    }
    return(E)
}
getphase <- function(e,omega,type='primary'){
##e: eccentricity
##omega: argument of ascending node
##T: true anomaly
    if(type=='primary') theta <- pi/2-omega
    if(type=='secondary') theta <- 3*pi/2-omega
    if(type=='periastron') theta <- 0
    if(type=='ascendingnode') theta <- -omega
    if(type=='descendingnode') theta <- pi-omega
    if(type=='l4') theta <- 5*pi/6-omega
    if(type=='l5') theta <- pi/6-omega
     E <- 2*atan(sqrt((1-e)/(1+e))*tan(theta/2))
     M <- E-e*sin(E)
     M/(2*pi)
}
getM0 <- function(e,omega,P,T,T0,type='primary'){
    Tp <- T-getphase(e,omega)*P
    ((T0-Tp)%%P)*2*pi/P
}
kep.mt2 <- function(m,e){
    tol  <-  1e-8
#    tol <- tol1
    E0 <- m
    Ntt <- 1e3
    for(k in 1:Ntt){
        E1 = E0-(E0-e*sin(E0)-m)/(sqrt((1-e*cos(E0))^2-(E0-e*sin(E0)-m)*(e*sin(E0))))
        if(all(abs(E1-E0)<tol)) break()
#        if(k==Ntt) cat('Keplerian solver does not converge:',e,m,E0,E1,'!\n')
        E0 <- E1
    }
    if(k==Ntt){
        cat('Keplerian solver does not converge!\n')
        cat('length(which(abs(E1-E0)>tol))=',length(which(abs(E1-E0)>tol)),'\n')
    }
    return(E1)
}
##refer to Murison A Practical Method for Solving the Kepler Equation
eps3 <- function(m,e,x){
    t1 <- cos(x)
    t2 <- -1+e*t1
    t3 <- sin(x)
    t4 <- e*t3
    t5 <- -x+t4+m
    t6 <- t5/(1/2*t5*t4/t2+t2)
    return(t5/((1/2*t3-1/6*t1*t6)*e*t6+t2))
}
KeplerStart3 <- function(m,e){
    t34 <- e^2
    t35 <- e*t34
    t33 <- cos(m)
    return(m+(-1/2*t35+e+(t34+3/2*t33*t35)*t33)*sin(m))
}
kep.murison2 <- function(m,e,tol=tol1){
    Mnorm <- m%%(2*pi)
    E0 <- KeplerStart3(Mnorm,e)
    Ntry <- 1000
    for(k in 1:Ntry){
        E1 <- E0-eps3(Mnorm,e,E0)
        if(k==Ntry) 'Kepler solver failed to converge!\n'
        if(all(abs(E1-E0)<tol)) break()
        E0 <- E1
    }
    return(E1)
}
##R package: uniroot
kep.R <- function(m,e){
    E <- rep(NA,length(m))
    for(j in 1:length(m)){
        y <- function(x) x-e*sin(x)-m[j]%%(2*pi)
        E[j] <- uniroot(y,interval=c(0,2*pi))$root
    }
    return(E)
}
divide.pars <- function(pars){
    #select Keplerian pars
    if(Np>0){
        ind.kep <- 1:(Nkeppar*Np)
    }else{
        ind.kep <- c()
    }
    #select noise parameters for different aperture
#    ind.noise <- c()
    ind.trend <- sort((1:length(pars))[grepl('\\da',names(pars)) | grepl('\\db',names(pars)) | grepl('\\ds',names(pars))])
    ind.noise <- (1:length(pars))[-c(ind.kep,ind.trend)]
    return(list(kep=ind.kep,trend=ind.trend,noise=ind.noise))
}

calc.eta <- function(m1,m2,band='G'){
###https://iopscience.iop.org/article/10.1088/0004-6256/145/3/81/pdf
###Arenou et al. 2000
###m1: primary; m2: secondary
###eta=(a1-a0)/a0 where a0 and a1 are the angular semimajor axes of the photocentric and the primary orbits.
    eta <- 0
    if(band=='Hp' & m2>0.1 & m1>0.1){
        dhp <- -13.5*(log10(m2)-log10(m1))
        eta <- (1+10^(0.074*dhp))/(10^(0.4*dhp)-10^(0.074*dhp))
    }else if(band=='G' & m1>mlow.m22 & m1<mup.m22 & m2>mlow.m22 & m2<mup.m22){
###http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
        dg <- mrl.m22(m2)-mrl.m22(m1)
        eta <- (1+m1/m2)/(10^(0.4*dg)-m1/m2)
    }
    eta
}

assign.names <- function(par.vector,Np,p=0,q=0,n=0,bases=rep('natural',10)){
    nams <- c()
    #names for keplerian part
    if(Np>0 | out$Nm>0){
        for(j in 1:(Np+out$Nm)){
            if(bases[j]=='natural'){
                if(prior.type=='e0'){
                    nams.kep <- c('per','K','omega')
                }else{
                    nams.kep <- c('per','K','e','omega','Mo')
                }
            }else if(bases[j]=='linear1'){
                if(prior.type=='e0'){
                    nams.kep <- c('per','lnK','Tc')
                }else{
                    nams.kep <- c('per','lnK','sqresinw','sqrecosw','Tc')
                }
            }else if(bases[j]=='linear2'){
                if(prior.type=='e0'){
                    nams.kep <- c('dP','lnK','dTc')
                }else{
                    nams.kep <- c('dP','lnK','sqresinw','sqrecosw','dTc')
                }
            }
            if(length(out$astrometry)>0){
                nams.kep <- c(nams.kep,'Inc','Omega')
            }
            if(j>Np){
                nams <- c(nams,paste0('m',nams.kep,j))
            }else{
                nams <- c(nams,paste0(nams.kep,j))
            }
        }
    }
    if(length(p)>0 & length(out$Nrv)>0){
        for(i1 in 1:length(p)){
            if(i1==1 & par.global=='trend'){
                nams <- c(nams,paste0('a',i1,1:Npoly))
            }
            if(offset) nams <- c(nams,paste0('b',i1))
            nams <- c(nams,paste0('J_',i1))
            if(p[i1]>0){
                nams.ar <- 'phi'
                for(j in 1:p[i1]){
                    nams <- c(nams,paste0(nams.ar,i1,j))
                }
                nams <- c(nams,paste0('alpha',i1))
            }
            if(q[i1]>0){
                nams.ma <- 'w'
                for(j in 1:q[i1]){
                    nams <- c(nams,paste0(nams.ma,i1,j))
                }
                nams <- c(nams,paste0('beta',i1))
            }
            if(ns[i1]>0){
                nams <- c(nams,paste0('c',i1,1:ns[i1]))
            }
        }
    }
#    if(any(names(out)=='astrometry') & !grepl('-',astro.type)){
    if(length(out$astrometry)>0){
        nams <- c(nams,'J_gaia')
        if(length(out$ins.epoch)>0){
            nams <- c(nams,paste0('J_',out$ins.epoch))
        }else if(!is.null(out$ihip) & !any(grepl('hip',out$ins.epoch))){
            nams <- c(nams,'J_hip')
        }
        if(length(out$gost)>0 | length(out$data.epoch)>0){
            nams <- c(nams,c('dra','ddec','dplx','dpmra','dpmdec'))
        }else{
            nams <- c(nams,c('dra','ddec','dpmra','dpmdec'))
        }
    }
#    if(any(names(out)=='rel'))  nams <- c(nams,c('J_rel'))
    if(any(names(out)=='rel')){
        nams <- c(nams,paste0('J_',rel.name),'Mstar')
    }
    if(out$Nrvc>0){
        if(rvc.type=='relative'){
            nams <- c(nams,paste0('bc',out$Irvc))
        }else{
            nams <- c(nams,'rvb')
        }
    }
    if(out$relativity){
        nams <- c(nams,'gamma')
        nn <- c('dra','ddec','dplx','dpmra','dpmdec','drv')
        ii <- which(is.na(match(nn,nams)))
        if(length(ii)>0) nams <- c(nams,nn[ii])
    }
####reflex motion
   if(length(out$data.epoch)>0){
       nams <- c(nams,paste0('J_',out$ins.epoch))
       ind <- which(!grepl('hip',out$ins.epoch))
       if(length(ind)>0){
           nams <- c(nams,paste0('dra_',out$ins.epoch[ind]))
           nams <- c(nams,paste0('ddec_',out$ins.epoch[ind]))
       }
   }

#   if(length(out$data.ref)>0){
#       nams <- c(nams,paste0('J_',out$ins.ref))
#   }
    if(photovary){
       nams <- c(nams,'eta')
    }
    if(length(out$data.binary)>0){
        if(offset) nams <- c(nams,paste0('b_',out$ins.binary))
        nams <- c(nams,paste0('s_',out$ins.binary))
    }
    if(FALSE){
#    if(TRUE){
        cat('length(nams)=',length(nams),'\n')
        cat('nams=',nams,'\n')
        cat('length(par.vector)=',length(par.vector),'\n')
    }
    names(par.vector) <- nams
    return(par.vector)
}
plot.labels.simple <- function(par.name,pn='b'){
    pat.kep <- c('per','K','e','omega','Mo','Omega','Inc','Tc','dP','dTc','lnK','sqrecosw','sqresinw','Mstar')
    unit.kep <- c('day','m/s','','rad','rad','rad','rad','day','day','day','','','','M[sun]')
    par.kep <- c(pat.kep,paste0('m',pat.kep))
    pat.noise <- c('a','b','beta','alpha','w','m','tau','s','c','J','logJ','jitter','dra','ddec','dpmra','dpmdec','dplx','drv')
    unit.noise <- c('m/s/yr','m/s','','','m/s','m/s','day','m/s','m/s','m/s','','m/s','mas','mas','mas/yr','mas/yr','mas','km/s')
    pats <- c(pat.kep,pat.noise)
    units <- c(unit.kep,unit.noise)
    new.name <- par.name
    for(i in 1:length(pats)){
        pat  <- pats[i]
        ind <- grep(paste0('^',pat),par.name)
        if(pat=='b') ind <- grep(paste0('^',pat,'[0-9]'),par.name)
        if(pat=='c') ind <- grep(paste0('^',pat,'[0-9]'),par.name)
        if(length(ind)>0){
            ii <- gsub(paste0('^',pat,'|^',pat,'\\.|0$'),'',par.name[ind])
            ips <- iups <- ilows <- c()
            for(k in 1:length(ii)){
                ss <- unlist(strsplit(ii[k],''))
                if(grepl('[1-9]$',ii[k])){
                    if(length(ss)>1){
                        iup <- ss[length(ss)-1]
                        ilow <- ss[length(ss)]
                    }else if(any(pat==pat.noise)){
                        ilow <- ''
                        iup <- ss
                        if(grepl('[1-9]',iup)) iup <- out$ins.rv[as.numeric(iup)]
#                        if(iup=='HARPSpre') iup <- 'pre'
#                        if(iup=='HARPSpost') iup <- 'post'
                    }else{
                        ilow <- ss
                        iup <- ''
                    }
                }else{
                    iup <- ii[k]
                    ilow <- ''
                }
                iups <- c(iups,iup)
                if(!any(pat.kep==pat)){
                    ilows <- c(ilows,ilow)
                    ips <- c(ips,'')
                }else{
                    ilows <- c(ilows,'')
                    ips <- c(ips,c('b','c','d','e','f','g','h','i','j')[as.numeric(ilow)])
#                    if(pat=='omega') ilow <- paste0("'*'",ilow)
                }
            }
#            if(grepl('jitter',pat)) pat <- 'J'
            if(grepl('J',pat)) pat <- 'J'
#            if(grepl('dra|ddec|dpmra|dpmdec',pat)) pat <- gsub('^d','Delta*',pat)
            if(pat=='per') pat <- 'P'
            if(pat=='Inc') pat <- 'I'
            if(pat=='b') pat <- 'gamma'
            if(pat=='Mo') pat <- 'M[0]'
            if(pat=='s') pat <- 'J'
            if(pat=='Mstar') pat <- "M['*']"
            if(pat=='dra') pat <- 'Delta*alpha'
            if(pat=='ddec') pat <- 'Delta*delta'
            if(pat=='dplx') pat <- 'Delta*varpi'
            if(pat=='dpmra') pat <- 'Delta*mu[alpha]'
            if(pat=='dpmdec') pat <- 'Delta*mu[delta]'
            if(pat=='drv') pat <- 'Delta*RV'
            for(j in 1:length(ind)){
                    if(ilows[j]!='' & iups[j]!='') new.name[ind[j]] <- paste0("",pat,"[",ilows[j],"]^{",iups[j],"}")
                    if(ilows[j]!='' & iups[j]=='') new.name[ind[j]] <- paste0(pat,' [',ilows[j],']')
                    if(ilows[j]=='' & iups[j]!='') new.name[ind[j]] <- paste0(pat,'^{',iups[j],'}')
                    if(ilows[j]=='' & iups[j]=='') new.name[ind[j]] <- pat
                    if(ips[j]!='') new.name[ind[j]] <- paste0(pat,' [',ips[j],']')
                    if(units[i]!='') new.name[ind[j]] <- paste0(new.name[ind[j]],"*' [",units[i],"]'")
            }
        }
    }
#    sapply(new.name, function(x) eval(parse(text=x)))
    new.name
#    nams <- vector('expression',length(new.name))
#    for(k in 1:length(nams)){
#        labs[k] <- substitute(lab.name,list(lab.name=nams[k]))
#    }
#    labs
}

cal.trend <- function(a,b,t){
    trend <- rep(b,length(t))
#cat('range(a)=',range(a),'\n')
#cat('range(b)=',range(b),'\n')
#cat('range(t)=',range(t),'\n')
    for(j3 in 1:length(a)){
        trend <- trend+a[j3]*t^j3#t[,j3]
    }
    return(trend)
}
extract.par <- function(pars.kep,bases=NULL){
    Np <- length(grep('dP|per|logP|^P',names(pars.kep)))
    if(is.null(bases)) bases <- rep('natural',Np)
    Ps <- Ks <- es <- Mos <- omegas <- Omegas <- Incs <- arcs <- apms <- c()
    bc <- list()
    if(Np>0){
        for(k in 1:(Np+out$Nm)){
            basis <- bases[k]
            if(basis=='natural'){
                if(k>Np){
                    P  <-  exp(pars.kep[grep(paste0('^mper',k),names(pars.kep))])
                    K <- pars.kep[grep(paste0('^mK',k),names(pars.kep))]
                    e <- pars.kep[grep(paste0('^me',k),names(pars.kep))]
                    Mo <- pars.kep[grep(paste0('^mMo',k),names(pars.kep))]
                    omega  <-  pars.kep[grep(paste0('^momega',k),names(pars.kep))]
                }else{
                    if(any(dP==k)){
                        P  <-  Prefs[k]+pars.kep[paste0('dP',k)]/60/24#day
                    }else{
                        P  <-  exp(pars.kep[grep(paste0('^per',k),names(pars.kep))])
                    }
                    K <- pars.kep[grep(paste0('^K',k),names(pars.kep))]
                    e <- pars.kep[grep(paste0('^e',k),names(pars.kep))]
                    Mo <- pars.kep[grep(paste0('^Mo',k),names(pars.kep))]
                    omega  <-  pars.kep[grep(paste0('^omega',k),names(pars.kep))]
                }
            }else{
                if(k>Np){
                    if(basis=='linear2'){
                        dP = pars.kep[grep(paste0('^mdP',k),names(pars.kep))]*alpha.dP
                        P <- dP+Ptransit[k]
                        dTc = pars.kep[grep('^mdTc([[:digit:]]{1})',names(pars.kep))]*alpha.dTc
                        Tc1 <- dTc+Tc[k]
                    }else if(basis=='linear1'){
                        P  <-  exp(pars.kep[grep(paste0('^mper',k),names(pars.kep))])
                        Tc1  <-  pars.kep[grep(paste0('^mTc',k),names(pars.kep))]
                    }
                    K <- exp(pars.kep[grep(paste0('^mlnK',k),names(pars.kep))])
                    sw  <-  pars.kep[grep(paste0('^msqresinw',k),names(pars.kep))]
                    cw  <-  pars.kep[grep(paste0('^msqrecosw',k),names(pars.kep))]
                }else{
                    if(basis=='linear2'){
                        dP = pars.kep[grep(paste0('^dP',k),names(pars.kep))]*alpha.dP
                        P <- dP+Ptransit[k]
                        dTc = pars.kep[grep('^dTc([[:digit:]]{1})',names(pars.kep))]*alpha.dTc
                        Tc1 <- dTc+Tc[k]
                    }else if(basis=='linear1'){
                        P  <-  exp(pars.kep[grep(paste0('^per',k),names(pars.kep))])
                        Tc1  <-  pars.kep[grep(paste0('^Tc',k),names(pars.kep))]
                    }
                    K <- exp(pars.kep[grep(paste0('^lnK',k),names(pars.kep))])
                    sw  <-  pars.kep[grep(paste0('^sqresinw',k),names(pars.kep))]
                    cw  <-  pars.kep[grep(paste0('^sqrecosw',k),names(pars.kep))]
                }
                e <- sw^2+cw^2
                if(cw!=0){
                    omega <- atan(sw/cw)
                }else{
                    omega <- atan(sw/1e-3)
                }
                Mo <- getM0(e=e,omega=omega,P=P,T=Tc1,T0=tref,type='primary')
            }
            if(k>Np){
                Omega <- pars.kep[grep(paste0('^mOmega',k),names(pars.kep))]
                Inc <- pars.kep[grep(paste0('^mInc',k),names(pars.kep))]
            }else{
                Omega <- pars.kep[grep(paste0('^Omega',k),names(pars.kep))]
                Inc <- pars.kep[grep(paste0('^Inc',k),names(pars.kep))]
            }
            arc <- apm <- c()
            if(any(grepl('^arc',names(pars.kep)))) arc <- pars.kep[grep(paste0('^arc',k),names(pars.kep))]
            if(any(grepl('^apm',names(pars.kep)))) apm <- pars.kep[grep(paste0('^apm',k),names(pars.kep))]
            Ps <- c(Ps,P)
            Ks <- c(Ks,K)
            es <- c(es,e)
            omegas <- c(omegas,omega)
            Omegas <- c(Omegas,Omega)
            Incs <- c(Incs,Inc)
            Mos <- c(Mos,Mo)
            arcs <- c(arcs,arc)
            apms <- c(apms,apm)
        }
    }
    oo <- list(P=Ps,K=Ks,e=es,Mo=Mos,omega=omegas,Omega=Omegas,Inc=Incs,arc=arcs,apm=apms)
    ind <- grep('bc',names(pars.kep))
    if(length(ind)>0){
        ns <- names(pars.kep)[ind]
        for(n in ns){
            oo[[n]] <- pars.kep[n]
        }
    }
    return(oo)
}

#calc.astro <- function(K,P,e,inc,omega,Omega,E,plx,Mstar,pars,eta=NA,state=FALSE,band='G',comp=1){
calc.astro <- function(par,index,E,plx,Mstar,pars,eta=NA,state=FALSE,band='G',comp=1){
####calculate orbital motion
    P <- par$P[index]
    e <- par$e[index]
    inc <- par$Inc[index]
    omega <- par$omega[index]
    Omega <- par$Omega[index]
    if(length(par$K)>0){
        K <- par$K[index]
        alpha0 <- K/sin(inc)/1e3/4.74047#au/yr
        beta0 <- P/365.25*(K*1e-3/4.74047)*sqrt(1-e^2)/(2*pi)/sin(inc)#au
        beta <- beta0*plx#mas
        T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
        alpha <- alpha0*plx#proper motion in mas/yr
        mp <- k2m(K,P,e,Mstar,Inc=inc)$ms
        if(all(is.na(eta))){
            if(any(names(pars)=='eta')){
                eta <- pars['eta']
            }else{
                eta <- calc.eta(Mstar,mp)
            }
        }
        xi <- 1/(eta+1)
        if(comp==2) xi <- -Mstar/mp
    }else{
        beta <- par$apm[index]#mas
        xi <- 1
        alpha <- alpha0 <- 0
    }
    ##semi-major axis is the astrometric signature in micro-arcsec
    A <- cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(inc)
    B <- sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(inc)
    F <- -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(inc)
    G <- -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(inc)
    C <- sin(omega)*sin(inc)
    H <- cos(omega)*sin(inc)
    ##    C <- sin(omega)*sin(inc)
    ##    H <- cos(omega)*sin(inc)
    Vx <- -sin(T)
    Vy <- cos(T)+e

###calculate POS
    X <- cos(E)-e
    Y <- sqrt(1-e^2)*sin(E)
    ##    T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    ##    alpha <- (P/365.25)^{2/3}*mp/(out$Mstar+mp)^{2/3}*plx#semi-major axis in mas
    ##    beta <- (P*3600*24)*K*sqrt(1-e^2)/(2*pi)*plx/sin(inc)*6.68459e-12#mas
    raP <- beta*(B*X+G*Y)
    decP <- beta*(A*X+F*Y)
    plxP <- -beta*(C*X+H*Y)*plx/206265e3#parallax change
    pmraP <- alpha*(B*Vx+G*Vy)
    pmdecP <- alpha*(A*Vx+F*Vy)
    ##    rvP.epoch <- alpha*(C*Vx+H*Vy)
    raP <- raP*xi
    decP <- decP*xi
    plxP <- plxP*xi
    pmraP <- pmraP*xi
    pmdecP <- pmdecP*xi
    if(FALSE){
#    if(TRUE){
        cat('E=',head(E),'\n')
        cat('C=',C,'\n')
        cat('H=',H,'\n')
        cat('X=',head(X),'\n')
        cat('Y=',head(Y),'\n')
        cat('eta=',eta,'\n')
        cat('beta0=',beta0,'\n')
        cat('beta=',beta,'\n')
        cat('eta=',eta,';xi=',xi,';mp=',mp,';Mstar=',Mstar,'\n')
    }
    rv <- alpha0*(C*Vx+H*Vy)#km/s
    if(state){
        BT <- cbind(raP/plx,decP/plx,beta0*(C*X+H*Y),pmraP/plx,pmdecP/plx,rv)#au,au/yr
        return(BT)
    }else{
        return(cbind(ra=raP,dec=decP,plx=plxP,pmra=pmraP,pmdec=pmdecP,rv=rv*4.74047))
    }
}

fm2m <- function(fm,inc,m1){
    m20 <- fm
    fm0 <- fm
    for(j in 1:100){
#        m2 <- sqrt((m20*sin(inc))^3/fm)-m1
        m2 <- ((m1+m20)^2*fm)^(1/3)/sin(inc)
        if(abs(m2-m20)<1e-6){
            break()
        }else{
            m20 <- m2
        }
    }
    m2
}

calcVp <- function(K,e,inc,w,Omega,E){
    if(inc==0) inc <- 1e-8
    T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    plx <- out$plx
    alpha <- K/sin(inc)/1e3/4.74047*plx#proper motion in mas/yr
    ##semi-major axis is the astrometric signature in micro-arcsec
    A <- alpha*(cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(inc))
    B <- alpha*(sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(inc))
    F <- alpha*(-cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(inc))
    G <- alpha*(-sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(inc))
    Vx <- -sin(T)
    Vy <- cos(T)+e
    pmraP <- B*Vx+G*Vy
    pmdecP <- A*Vx+F*Vy
    return(cbind(pmraP,pmdecP))
}

pqu <- function(alpha,delta){
###alpha: RA in rad
###delta: DEC in rad
###ref. to https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu3ast/sec_cu3ast_cali/ssec_cu3ast_cali_source.html
    p <- cbind(-sin(alpha),cos(alpha),0)
    q <- cbind(-sin(delta)*cos(alpha),-sin(delta)*sin(alpha),cos(delta))
    u <- cbind(cos(delta)*cos(alpha),cos(delta)*sin(alpha),sin(delta))
    list(p=p,q=q,u=u)
}

calcRp <- function(K,P,e,inc,w,Omega,E){
    if(inc==0) inc <- 1e-8
    X <- cos(E)-e
    Y <- sqrt(1-e^2)*sin(E)
    plx <- out$plx
    mp <- k2m(K,P,e,Mstar,Inc=inc)$ms#in unit of solar mass
#    T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    alpha <- (P/365.25)^{2/3}*mp/(Mstar+mp)^{2/3}*plx#semi-major axis in mas
    ##semi-major axis is the astrometric signature in micro-arcsec
    A <- alpha*(cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(inc))
    B <- alpha*(sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(inc))
    F <- alpha*(-cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(inc))
    G <- alpha*(-sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(inc))
    raP <- B*X+G*Y
    decP <- A*X+F*Y
    return(cbind(raP,decP))
}

calcMo <- function(Ps,es,Mos,tt2,dt){
    Np <- length(Ps)
    ms <- c()
    for(j in 1:Np){
        e <- es[j]
        Mo <- (Mos[j]+2*pi*dt/Ps[j])%%(2*pi)#mean anomaly of the astrometry reference point
        m <- (Mo+2*pi*tt2/Ps[j])%%(2*pi)#mean anomaly
        ms <- cbind(ms,m)
    }
    return(ms)
}
pa <- function(obs,t){
##t: in unit of days
##PM PA ref: The Hundred Thousand Proper Motions Project: Mignard 2009; GAIA-C3-TN-OCA-FM-040-01;
##https://arxiv.org/pdf/1208.3048.pdf
##parallax PA ref: https://www.aanda.org/articles/aa/pdf/2021/08/aa41344-21.pdf
##RV PA ref:https://www.aanda.org/articles/aa/pdf/2009/38/aa12479-09.pdf
##when comparing pa(...) with vectorized propagations, please use the following:
##db <- bary0-bary1;db[,'ra'] <- db[,'ra']*3.6e6*cos(bary0[,'dec']/180*pi);db[,'dec'] <- db[,'dec']*3.6e6;db[,'pmra'] <- db[,'pmra']*cos(bary0[,'dec']/180*pi)
##where bary0 considered PA and bary0 does not.
###output unit: mas or mas/yr
    if(!is.null(obs)) obs <- obs[1,]
    ra <- as.numeric(obs['ra'])/180*pi#rad
    dec <- as.numeric(obs['dec'])/180*pi
    plx <- as.numeric(obs['parallax'])
    pmra <- as.numeric(obs['pmra'])
    pmdec <- as.numeric(obs['pmdec'])
    rv <- as.numeric(obs['radial_velocity'])
    pmrv <- rv/4.74047*plx#mas
#    qra <- -2*(pmrv*pmra-tan(dec)*pmra*pmdec)
    qra <- -2*pmrv*pmra
#    qdec <- -2*(pmrv*pmdec+0.5*tan(dec)*pmra^2)
    qdec <- -2*pmrv*pmdec
    r2m <- 206265e3#1rad/1mas
    dT <- t/365.25#yr
    dra <- 0.5*qra*dT^2/r2m#dalpha*
    ddec <- 0.5*qdec*dT^2/r2m#rad
    dplx <- -plx*pmrv*dT/r2m
    dpmra <- qra*dT/r2m
    dpmdec <- qdec*dT/r2m
    pm <- sqrt(pmra^2+pmdec^2)*1e-3#arcsec/yr
    drv <- 22.98*pm^2/plx*t/365.25#m/s/yr*yr=m/s
    tmp <- cbind(dra,ddec,dplx,dpmra,dpmdec,drv)
    colnames(tmp) <- c('dra','ddec','dplx','dpmra','dpmdec','drv')
    tmp
}

obs.lin.prop <- function(obs,t,PA=TRUE){
##whether to consider perspective acceleration
##t is in unit of day
##obs: ra(deg), dec(deg),plx(mas),pmra(mas/yr),pmdec(mas/yr)
    kpcmyr2auyr <- 1e3*206265/1e6
    pc2au <- 206265
    kpcmyr2kms <- kpcmyr2auyr*auyr2kms
    ra <- as.numeric(obs['ra'])/180*pi#rad
    dec <- as.numeric(obs['dec'])/180*pi
    plx <- as.numeric(obs['parallax'])
    pmra <- as.numeric(obs['pmra'])
    pmdec <- as.numeric(obs['pmdec'])
    rv <- as.numeric(obs['radial_velocity'])
    if(PA){
#### obs-> initial state: propagation observables to states
        d <- 1/plx#kpc
        r <- as.numeric(bl2xyz(dec,ra))*d*1e3
        x <- r[1]
        y <- r[2]
        z <- r[3]
        vde <- pmdec*d
        vra <- pmra*d
        vp <- sqrt(vra^2+vde^2)
        vr <- rv/auyr2kms#au/yr
        vx.equ <- vr*cos(dec)*cos(ra)-vde*sin(dec)*cos(ra)-vra*sin(ra)##note: vr is positive if the star is moving away from the Sun
        vy.equ <- vr*cos(dec)*sin(ra)-vde*sin(dec)*sin(ra)+vra*cos(ra)
        vz.equ <- vr*sin(dec)+vde*cos(dec)
        x1 <- x+vx.equ*t/365.25/pc2au
        y1 <- y+vy.equ*t/365.25/pc2au
        z1 <- z+vz.equ*t/365.25/pc2au

### propagation: convert time-varying states back to observables
        dec.ra <- xyz2bl.vec(x1,y1,z1)
        d1 <- sqrt(x1^2+y1^2+z1^2)*1e-3#kpc
        ra1.rad <- dec.ra[,2]#rad
        dec1.rad <- dec.ra[,1]
        ra1 <- ra1.rad*180/pi#deg
        dec1 <- dec1.rad*180/pi

###state -> obs: velocity to pm
        vequ <- array(NA,dim=c(length(t),3))
        for(j in 1:length(t)){
            rotz <- matrix(data=c(cos(ra1.rad[j]),sin(ra1.rad[j]),0.0,-sin(ra1.rad[j]),cos(ra1.rad[j]),0.0,0.0,0.0,1.0),nrow=3,ncol=3,byrow=TRUE)#o-xyz -> o-x'y'z'
            roty <- matrix(data=c(cos(dec1.rad[j]),0.0,sin(dec1.rad[j]),0.0,1.0,0.0,-sin(dec1.rad[j]),0.0,cos(dec1.rad[j])),nrow=3,ncol=3,byrow=TRUE)
            vequ[j,] <- roty%*%rotz%*%as.numeric(c(vx.equ,vy.equ,vz.equ))
        }
        pmra1 <- vequ[,2]/d1#mas/yr
        pmdec1 <- vequ[,3]/d1#mas/yr
        rv1 <- vequ[,1]*auyr2kms
        out <- cbind(ra1,dec1,1/d1,pmra1,pmdec1,rv1)
    }else{
        decs <- dec+pmdec*t/365.25/206265e3#rad
        ras <- ra+pmra*t/365.25/cos(decs)/206265e3#rad
        out <- cbind(ras*180/pi,decs*180/pi,rep(plx,length(t)),rep(pmra,length(t)),rep(pmdec,length(t)),rep(rv,length(t)))
    }
    colnames(out) <- c('ra','dec','parallax','pmra','pmdec','radial_velocity')
    return(out)
}

###calculate UWE according to the formulae given by El-Badry+2024
calc.uwe <- function(dabs,sfov,Nfov,Nbin,Npar=5){
##dabs: abscissa residual;
##sfov: FOV transit uncertainty, including photon noise and attitude noise
    chi2.bin <- sum(dabs^2/sfov^2)*Nfov/length(dabs)#scaling to consider bad FOV transits
    chi2.unbin <- chi2.bin+Nfov*(Nbin-1)
    sqrt(chi2.unbin/(Nfov*Nbin-Npar))
}

vec.lin.prop <- function(obs,t){
##t is in unit of day
##obs: ra(deg), dec(deg),plx(mas),pmra(mas/yr),pmdec(mas/yr)
    t <- t/365.25#yr
    kpcmyr2auyr <- 1e3*206265/1e6
    pc2au <- 206265
    kpcmyr2kms <- kpcmyr2auyr*auyr2kms
    if(is.null(dim(obs)))  obs <- t(obs)
    if(nrow(obs)==1 & length(t)>1){
        obs <- t(replicate(length(t),obs[1,]))
    }
    ra <- obs[,'ra']/180*pi#rad
    dec <- obs[,'dec']/180*pi#rad
    plx <- obs[,'parallax']#mas
    pmra <- obs[,'pmra']#mas/yr
    pmdec <- obs[,'pmdec']#mas/yr
    rv <- obs[,'radial_velocity']#km/s
#cat('rv=',rv,'\n')
#### obs-> initial state: propagation observables to states
    d <- 1/plx#kpc
    r <- bl2xyz(dec,ra)*d*1e3#pc
    vde <- pmdec*d
    vra <- pmra*d
    vp <- sqrt(vra^2+vde^2)
    vr <- rv/auyr2kms#au/yr
    C <- pqu(ra,dec)
    v <- vra*C$p+vde*C$q+vr*C$u#velocity in cartecian coordinates
    r1 <- r+v*t/pc2au

###propagation: convert time-varying states back to observables
    dec.ra <- xyz2bl.vec(r1[,1],r1[,2],r1[,3])
    d1 <- sqrt(rowSums(r1^2))*1e-3#kpc
    ra1.rad <- dec.ra[,2]#rad
    dec1.rad <- dec.ra[,1]
    ra1 <- ra1.rad*180/pi#deg
    dec1 <- dec1.rad*180/pi

###state -> obs: velocity to pm
    C1 <- pqu(ra1.rad,dec1.rad)
    vequ <- cbind(rowSums(v*C1$p),rowSums(v*C1$q),rowSums(v*C1$u))
    pmra1 <- vequ[,1]/d1#mas/yr
    pmdec1 <- vequ[,2]/d1#mas/yr
    rv1 <- vequ[,3]*auyr2kms
    out <- cbind(ra1,dec1,1/d1,pmra1,pmdec1,rv1)
    colnames(out) <- c('ra','dec','parallax','pmra','pmdec','radial_velocity')
    return(out)
}

hip2gaia.mag <- function(V,VI){
    V-0.0257-0.0924*VI-0.1623*VI^2+0.009*VI^3
}

show.par <- function(par){
    cat(deparse(substitute(par)),'=',head(unlist(par)),'\n')
}

####Note that bl2xyz return coordinates in the heliocentric frame and the x axis point to the GC, yaxis point to the rotation of the Galaxy
bl2xyz <- function(b.rad,l.rad){
    x <- cos(b.rad)*cos(l.rad)
    y <- cos(b.rad)*sin(l.rad)
    z <- sin(b.rad)
    return(cbind(x,y,z))
}
binary.orbit <- function(pars.kep,tt=NULL,bases=rep('natural',10),Ein=NULL){
    if(!any(names(pars.kep)=='Mstar')){
        Mstar <- out$Mstar
    }else{
        Mstar <- pars.kep['Mstar']
    }
    planet <- astrometry.kepler(pars.kep,tt=tt,bases=bases)$planet
    Nsig <- length(grep('^per',names(pars.kep)))
    astro <- list()
    for(j in 1:Nsig){
        Popt <- exp(par.opt[paste0('per',j)])#day
        Mp <- k2m(K,P,e,Mstar,Inc=par.opt[paste0('Inc',j)])
        eta <- Mp/(Mstar+Mp)
#        ra <- -planet$ra/eta
#        dec <- -planet$dec/eta
#        pmra <- -planet$pmra/eta
#        pmdec <- -planet$pmdec/eta
#        astro[['ra']][[paste0('p',j)]] <- ra
#        astro[['dec']][[paste0('p',j)]] <- dec
#        astro[['pmra']][[paste0('p',j)]] <- pmra
#        astro[['pmdec']][[paste0('p',j)]] <- pmdec
        astro[[paste0('p',j)]] <- -planet/eta
    }
    astro
}
astrometry.rel <- function(pars.kep,tt=NULL,bases=rep('natural',10),pp=NULL,eta=0,band='G'){
    if(is.null(pp))  pp <- extract.par(pars.kep,bases=bases)
    pmraB <- pmdecB <- raB <- decB <- list()
    pmraR <- pmdecR <- raR <- decR <- c()
    Np.kep <- length(grep('omega',names(pars.kep)))
    if(!any(names(pars.kep)=='Mstar')){
        Mstar <- out$Mstar
    }else{
        Mstar <- pars.kep['Mstar']
    }
    if(Np.kep>0){
        ii <- sort(pp$P,index.return=TRUE)$ix
        for(kk in 1:length(ii)){
            j  <-  ii[kk]
####add inner planet mass to mass center
            mm <- Mstar
            if(kk>1){
                jj <- ii[1:(kk-1)]
                mp.tot <- sum(k2m(pp$K[jj],pp$P[jj],pp$e[jj],Ms=Mstar,Inc=pp$Inc[jj])$ms)
                mm <- mm + mp.tot
            }

####relative astrometry reference: https://arxiv.org/pdf/2109.10671.pdf
            if(!is.null(tt)){
                ms <- (pp$Mo[j]+2*pi*(tt-out$tref)/pp$P[j])%%(2*pi)
            }else{
                ms <- (pp$Mo[j]+2*pi*(out$trel-out$tref)/pp$P[j])%%(2*pi)
            }
            E <- kep.mt2(ms,pp$e[j])
            plx <- out$plx
            if(any(names(pars.kep)=='dplx')) plx <- plx-pars.kep['dplx']
#            parin <- c(pp$K[j],pp$P[j],pp$e[j],pp$Inc[j],pp$omega[j],pp$Omega[j],E,plx=plx,Mstar=mm,pars=pars.kep,eta=eta)
#            cat(parin,'\n')
            tmp <- calc.astro(pp,index=j,E,plx=plx,Mstar=mm,pars=pars.kep,eta=eta)
            Mp <- k2m(pp$K[j],pp$P[j],pp$e[j],Ms=mm,Inc=pp$Inc[j])$ms
            xi <- Mp/(mm+Mp)
            raR <- cbind(raR,tmp[,1])#reflex motion
            decR <- cbind(decR,tmp[,2])
            pmraR <- cbind(pmraR,tmp[,4])#reflex motion
            pmdecR <- cbind(pmdecR,tmp[,5])
            n <- paste0('p',j)
            inss <- names(out$indrel[[n]])
            if(!is.null(tt)){
                raB[[n]] <- -tmp[,1]/xi
                decB[[n]] <- -tmp[,2]/xi
                pmraB[[n]] <- -tmp[,4]/xi
                pmdecB[[n]] <- -tmp[,5]/xi
                if(kk>1){
                    raB[[n]] <- raB[[n]] - rowSums(raR[,1:(kk-1),drop=FALSE])
                    decB[[n]] <- decB[[n]] - rowSums(decR[,1:(kk-1),drop=FALSE])
                    pmraB[[n]] <- pmraB[[n]] - rowSums(pmraR[,1:(kk-1),drop=FALSE])
                    pmdecB[[n]] <- pmdecB[[n]] - rowSums(pmdecR[,1:(kk-1),drop=FALSE])
                }
            }else if(length(out$indrel[[n]])>0){
                for(i in inss){
                    raB[[n]][[i]] <- -tmp[out$indrel[[n]][[i]],1]/xi
                    decB[[n]][[i]] <- -tmp[out$indrel[[n]][[i]],2]/xi
                    pmraB[[n]][[i]] <- -tmp[out$indrel[[n]][[i]],4]/xi
                    pmdecB[[n]][[i]] <- -tmp[out$indrel[[n]][[i]],5]/xi
                    if(kk>1){
                        raB[[n]][[i]] <- raB[[n]][[i]] - rowSums(raR[out$indrel[[n]][[i]],1:(kk-1),drop=FALSE])
                        decB[[n]][[i]] <- decB[[n]][[i]] - rowSums(decR[out$indrel[[n]][[i]],1:(kk-1),drop=FALSE])
                        pmraB[[n]][[i]] <- pmraB[[n]][[i]]- rowSums(pmraR[out$indrel[[n]][[i]],1:(kk-1),drop=FALSE])
                        pmdecB[[n]][[i]] <- pmdecB[[n]][[i]] - rowSums(pmdecR[out$indrel[[n]][[i]],1:(kk-1),drop=FALSE])
                    }
                }
            }
        }
    }
    list(ref=list(dra=raR,ddec=decR,dpmra=pmraR,dpmdec=pmdecR),rel=list(dra=raB,ddec=decB,dpmra=pmraB,dpmdec=pmdecB))
}

astrometry.ref <- function(pars.kep,tt=NULL,bases=rep('natural',10),pp=NULL,eta=0,Pmin=0,band='G'){
    if(is.null(pp))  pp <- extract.par(pars.kep,bases=bases)
    ra.planet <- dec.planet <- pmra.planet <- pmdec.planet <- 0
    Np.kep <- length(grep('omega',names(pars.kep)))
    if(!any(names(pars.kep)=='Mstar')){
        Mstar <- out$Mstar
    }else{
        Mstar <- pars.kep['Mstar']
    }
    plx <- out$plx
    if(any(names(pars.kep)=='dplx')) plx <- plx-pars.kep['dplx']
    dpmdec.reflex <- dpmra.reflex <- ddec.reflex <- dra.reflex <- dplx.reflex <- 0
    if(Np.kep>0){
        for(j in 1:Np.kep){
            if(pp$P[j]>Pmin){
                if(!is.null(tt)){
                    ms <- (pp$Mo[j]+2*pi*(tt-out$tref)/pp$P[j])%%(2*pi)
                }else{
                    ms <- (pp$Mo[j]+2*pi*(out$data.ref[,1]-out$tref)/pp$P[j])%%(2*pi)
                }
                E <- kep.mt2(ms,pp$e[j])
                tmp <- calc.astro(pp,j,E,plx=plx,Mstar=Mstar,pars=pars.kep,band=band,eta=eta0[j])
                dra.reflex <- dra.reflex+tmp[,1]
                ddec.reflex <- ddec.reflex+tmp[,2]
                dplx.reflex <- dplx.reflex+tmp[,3]
                dpmra.reflex <- dpmra.reflex+tmp[,4]
                dpmdec.reflex <- dpmdec.reflex+tmp[,5]
            }
        }
    }
    list(ref=data.frame(dra=dra.reflex,ddec=ddec.reflex,dplx=dplx.reflex,dpmra=dpmra.reflex,dpmdec=dpmdec.reflex))
}

astrometry.epoch <- function(pars.kep,tt=NULL,bases=rep('natural',10),pp=NULL,eta=0,barycenter=NULL,band='G',comp=1){
    if(is.null(pp))  pp <- extract.par(pars.kep,bases=bases)
    ra.planet <- dec.planet <- pmra.planet <- pmdec.planet <- 0
    Np.kep <- length(grep('omega',names(pars.kep)))

    if(!any(names(pars.kep)=='Mstar')){
        Mstar <- out$Mstar
    }else{
        Mstar <- pars.kep['Mstar']
    }
    plx <- out$plx
    if(any(names(pars.kep)=='dplx')) plx <- plx-pars.kep['dplx']
    drv.reflex <- dpmdec.reflex <- dpmra.reflex <- dplx.reflex <- ddec.reflex <- dra.reflex <- mps <- 0
    if(Np.kep>0){
        for(j in 1:Np.kep){
                ms <- (pp$Mo[j]+2*pi*(tt-out$tref)/pp$P[j])%%(2*pi)
                E <- kep.mt2(ms,pp$e[j])
                if(band=='G') eta <- etaG[j]
                if(band=='H') eta <- etaH[j]
                tmp <- calc.astro(pp,j,E,plx=plx,Mstar=Mstar,pars=pars.kep,band=band,eta=eta,comp=comp)
                dra.reflex <- dra.reflex+tmp[,1]
                ddec.reflex <- ddec.reflex+tmp[,2]
                dplx.reflex <- dplx.reflex+tmp[,3]
                dpmra.reflex <- dpmra.reflex+tmp[,4]
                dpmdec.reflex <- dpmdec.reflex+tmp[,5]
                drv.reflex <- drv.reflex+tmp[,6]#km/s
        }
    }
    list(epoch=data.frame(dra=dra.reflex,ddec=ddec.reflex,dplx=dplx.reflex,dpmra=dpmra.reflex,dpmdec=dpmdec.reflex,drv=drv.reflex))
}

delay <- function(pars.kep,t,E,Mstar=NULL,bases=rep('natural',10)){
    pp <- extract.par(pars.kep,bases=bases)
    plx <- out$plx
    if(any(names(pars.kep)=='dplx')) plx <- plx-pars.kep['dplx']
    Np.kep <- length(grep('omega',names(pars.kep)))
#####iterations to calulate BT
    roemerT <- 0
    for(j in 1:Np.kep){
###calculate BT
        if(length(pp$K)>0){
            bt <- calc.astro(pp,j,E[[j]],plx=plx,Mstar=Mstar,state=TRUE,eta=0,pars=pars.kep)
            roemerT <- roemerT+bt[,3]/Cauday
        }else{
###https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.4136B/abstract
            ltte <- pp$arc[j]*sin(pp$Inc[j])*(sqrt(1-pp$e[j]^2)*sin(E[[j]])*cos(pp$omega[j])+(cos(E[[j]])-pp$e[j])*sin(pp$omega[j]))/3600/24#s to day
            roemerT <- roemerT+ltte
        }
    }
    tauT <- t-roemerT#barycentric time to emission time; without relativitistic correction
    return(data.frame(tau=tauT,roemer=roemerT))
}

###This is to account for the relativitiy-induced RV variation
RV.relativity <- function(pars.kep,t=NULL,indP=NULL,indC=NULL,Mstar=NULL,Mps=NULL,as=NULL,bases=rep('natural',10),pp=NULL,gamma=1,sim.kep=FALSE,tol=1e-5){
####initial settings
    if(is.null(pp))  pp <- extract.par(pars.kep,bases=bases)
    binarity <- FALSE
    if(length(out$data.binary)>0) binarity <- TRUE
    plx <- out$plx
    if(any(names(pars.kep)=='dplx')) plx <- plx-pars.kep['dplx']
    rv <- 0
    Np.kep <- length(grep('omega',names(pars.kep)))
    BT <- 0
    EC <- NULL
    CT <- list()
    EP <- BT.all <- list()
    if(any(names(pars.kep)=='gamma')){
        gamma <- pars.kep['gamma']
#        cat('gamma=',gamma,'\n')
    }
    imax <- which.max(pp$P)
    for(j in 1:Np.kep){
        m <- (pp$Mo[j]+2*pi*(t[indP]-out$tref)/pp$P[j])%%(2*pi)
        EP[[j]] <- kep.mt2(m,pp$e[j])%%(2*pi)
    }
    if(binarity){
        m <- (pp$Mo[imax]+2*pi*(t[indC]-out$tref)/pp$P[imax])%%(2*pi)
        EC <- kep.mt2(m,pp$e[imax])%%(2*pi)
    }
#####iterations to calulate BT
    Ntry <- 10
    roemerC <- roemerT <- roemerT0 <- roemerC0 <- 0
    TC <- BC <- BT.all <- NULL
    for(k in 1:Ntry){
        roemerT <- 0
        for(j in 1:Np.kep){
###calculate BT
            if(pp$K[j]>0){
                bt <- calc.astro(pp,j,EP[[j]],plx=plx,Mstar=Mstar,state=TRUE,eta=0,pars=pars.kep)
                BT.all[[j]] <- bt
                roemerT <- roemerT+bt[,3]/Cauday
            }else{
                BT.all[[j]] <- 0
            }
###decide whether to repeat
        }
        if(binarity){
            bt <- calc.astro(pp,imax,EC,plx=plx,Mstar=Mstar,state=TRUE,eta=0,pars=pars.kep)
            BC <- -bt*Mstar/Mps[imax]
            roemerC <- BC[,3]/Cauday
        }
        if(max(roemerT-roemerT0)<tol & max(roemerC-roemerC0)<tol) break()
        roemerT0 <- roemerT
        roemerC0 <- roemerC
    }

####Roemer delay
    tau <- t#barycentric time
    tauC <- NULL
#    cat('roemerT=',sd(roemerT),'\n')
    tauT <- t[indP]-roemerT#barycentric time to emission time; without relativitistic correction
    if(binarity){
        if(indP[1]==indC[1]){
            tauC <- t-roemerC#barycentric time to emission time; without relativitistic correction
        }else{
            tauC <- tau[indC] <- t[indC]-roemerC
        }
    }
####eccentrici anomalty
    BT <- 0
    for(j in 1:Np.kep){
            m <- (pp$Mo[j]+2*pi*(tauT-out$tref)/pp$P[j])%%(2*pi)
            EP[[j]] <- kep.mt2(m,pp$e[j])%%(2*pi)
            if(!is.null(BT.all)) BT <- BT+BT.all[[j]]
    }
    if(binarity){
        m <- (pp$Mo[imax]+2*pi*(tauC-out$tref)/pp$P[imax])%%(2*pi)#assume that the first companion is identified as stellar binary component
        EC <- kep.mt2(m,pp$e[imax])%%(2*pi)
    }

####relativistic effect
    zgP <- c()
    for(j in 1:Np.kep){
        if(pp$K[j]>0){
            eta <- Mps[j]/(Mstar+Mps[j])
            CT[[j]] <- BT.all[[j]]/eta
            RCT <- sqrt(rowSums(CT[[j]]^2))
            zgP <- cbind(zgP,SRS*Mps[j]*(1+gamma)/4*(1/RCT-1/as[j]))#relative RV
        }else{
            CT[j] <- list(NULL)
#            Eastro[j] <- list(NULL)
        }
    }

    RBT <- rowSums(BT[,1:3]^2)#au
#    VBT <- rowSums(BT[,4:6]^2)#au/yr
    if(length(out$astrometry)==0){
        vBx <- out$pmra/plx#au/yr
        vBy <- out$pmdec/plx#au/yr
        vBz <- out$srv/auyr2kms#au/yr
    }else{
        vBx <- out$astrometry[out$iref,'pmra']/plx#au/yr; in p direction
        vBy <- out$astrometry[out$iref,'pmdec']/plx#au/yr; in q direction
        vBz <- out$astrometry[out$iref,'radial_velocity']/auyr2kms#au/yr; in u direction
    }
    VST2 <- (vBx+BT[,4])^2+(vBy+BT[,5])^2+(vBz+BT[,6])^2#relative velocity
    VSB2 <- vBx^2+vBy^2+vBz^2
###relative relativistc RV for primary
    rvsT <- gamma*(VST2-VSB2)/Cauyr^2/2*CMPS#m/s
    rvgT <- rowSums(zgP)*CMPS#m/s
    rvgsT <- rvgT+rvsT
###relative relativistc RV for binary companion
    if(binarity){
        TC <- -bt*(Mstar+Mps[imax])/Mps[imax]
        RTC <- sqrt(rowSums(TC^2))
        zgC <- SRS*Mstar*(1+gamma)/4*(1/RTC-1/as[imax])#GR relativistic
        rvgC <- zgC*CMPS#m/s
        VSC2 <- (vBx+BC[,4])^2+(vBy+BC[,5])^2+(vBz+BC[,6])^2#relative velocity
        rvsC <- gamma*(VSC2-VSB2)/Cauyr^2/2*CMPS#m/s
        rvgsC <- rvgC+rvsC
    }else{
        rvsC <- rvgC <- rvgsC <- 0
    }
###calibrate the radial velocity by accounting for astrometric bias
    drvT <-  drvC <- 0
    if(any(grepl('dra$',names(pars.kep))) & length(out$astrometry)>0){
        obs0 <- unlist(out$astrometry[out$iref,c('ra','dec','parallax','pmra','pmdec','radial_velocity')])#reference astrometry
        obs <- obs0-c(pars.kep[c('dra','ddec')]/3.6e6,pars.kep[c('dplx','dpmra','dpmdec','drv')])#corrected astrometry
        ast1 <- obs.lin.prop(obs0,tau-out$astrometry[out$iref,'ref_epoch'])
        ast2 <- obs.lin.prop(obs,tau-out$astrometry[out$iref,'ref_epoch'])
        dast <- ast2-ast1
        du <- t(outer(out$pqu[,1],dast[,1]/180*pi*cos(obs['dec']/180*pi),'*'))+t(outer(out$pqu[,2],dast[,2]/180*pi,'*'))
###add reflex motion
        duP <- du[indP,]+BT[,1:3]/pc2au/(1000/out$astrometry[out$iref,'parallax'])#change of direction in rad
###add reflex motion to binary companion
        if(binarity){
            duC <- du[indC,]+BC[,1:3]/pc2au/(1000/out$astrometry[out$iref,'parallax'])#change of direction in rad
        }else{
            duC <- 0
        }
        if(!sim.kep){
            drvT <- rowSums(out$vSO[indP,]*duP*1e3)#m/s
            if(binarity) drvC <- rowSums(out$vSO[indC,]*duC*1e3)#m/s
        }else{
            drvT <- rowSums(out$vSO.sim*duP*1e3)
            if(binarity) drvC <- rowSums(out$vSO.sim*duC*1e3)
        }
#        cat('range(duP)=',range(duP),'\n')
#        cat('range(out$vSO.sim)=',range(out$vSO.sim),'km/s\n')
    }
    return(list(rvgsT=rvgsT,rvgT=rvgT,rvsT=rvsT,tau=tau,tauT=tauT,tauC=tauC,roemerT=roemerT,drvT=drvT,rvgsC=rvgsC,rvgC=rvgC,rvsC=rvsC,EP=EP,EC=EC,roemerC=roemerC,drvC=drvC))
}

astrometry.bary <- function(pars.kep,tt=NULL,bases=rep('natural',10),Ein=NULL,pa=FALSE){
###This is to give both absolute and relative astrometry if the absolute astrometry data is given
    data.astrometry <- out$astrometry
    if(is.null(tt)){
        if(!grepl('hgca',astro.type)){
            tt <- out$astrometry[,1]
        }else{
            tt <- c(as.numeric(out$astrometry[1,c('epoch_ra','epoch_dec')]),as.numeric(out$astrometry[2,c('epoch_ra','epoch_dec')]))
        }
        sim <- FALSE
    }else{
        sim <- TRUE
    }
    dt <- data.astrometry[out$iref,1]-out$tref#epoch offset between astrometry reference point and the RV reference point
#    if(is.null(tt)) tt <- data.astrometry[,'ref_epoch']
    DT <- tt - data.astrometry[out$iref,1]#relative to astrometry reference point
    if(grepl('hgca',astro.type)){
        dt <- as.numeric(data.astrometry[2,c('epoch_ra','epoch_dec')])-out$tref
        DT <- as.numeric(tt[1:2]-tt[3:4])
    }
    plx1 <- pmdec1 <- pmra1 <- rep(0,length(tt))
    Np.kep <- length(grep('per|logP|lnP|^P',names(pars.kep)))
    pp <- extract.par(pars.kep,bases=bases)
    es <- pp$e
    Omegas <- pp$Omega
    ws <- pp$omega
    Incs <- pp$Inc
    Ps <- pp$P
    Ks <- pp$K
    Mos <- pp$Mo
    dra <- ddec <- dpmra <- dpmdec <- 0
    dplx <- drv <- 0
    if(any(names(pars.kep)=='dra')) dra <- dra+pars.kep['dra']#as; note that this dra0 is ra
    if(any(names(pars.kep)=='ddec')) ddec <- ddec+pars.kep['ddec']
    if(any(names(pars.kep)=='dplx')) dplx <- dplx+pars.kep['dplx']
    if(any(names(pars.kep)=='dpmra')) dpmra <- dpmra+pars.kep['dpmra']#mas/yr
    if(any(names(pars.kep)=='dpmdec')) dpmdec <- dpmdec+pars.kep['dpmdec']
    if(any(names(pars.kep)=='drv')) drv <- drv+pars.kep['drv']
    if(any(names(out)=='rel')){
        Mstar <- pars.kep['Mstar']
    }else{
        Mstar <- out$Mstar
    }
###model parallax and mean proper motion; only assumption: constant heliocentric velocity; reference Gaia epoch
    obs <- unlist(data.astrometry[out$iref,c('ra','dec','parallax','pmra','pmdec','radial_velocity')])
###subtract the offset position and PM to get the initial condition for the barycentric motion
    obs['ra'] <- obs['ra']-dra/3.6e6/cos(obs['dec']/180*pi)#mas to deg
    obs['dec'] <- obs['dec']-ddec/3.6e6
    obs['parallax'] <- obs['parallax']-dplx#mas
    obs['pmra'] <- obs['pmra']-dpmra
    obs['pmdec'] <- obs['pmdec']-dpmdec
    obs['radial_velocity'] <- obs['radial_velocity']-drv#km/s

##propagation barycentric observables
    obs.lin.prop(obs,DT)
}

###calculate the astrometric difference between two epochs
AstroDiff <- function(obs1,obs2){
###obs1, obs2: ra[deg], dec[deg], parallax [mas], pmra [mas/yr], pmdec [mas/yr], rv [km/s]
    astro.name <- c('ra','dec','parallax','pmra','pmdec','radial_velocity')
    o1 <- as.numeric(obs1[astro.name])
    o2 <- as.numeric(obs2[astro.name])
    dobs <- o2-o1
    dobs[1:2] <- dobs[1:2]*3.6e6#deg to mas
    dobs[1] <- dobs[1]*cos(mean(c(o1[2],o2[2]))*pi/180)
    names(dobs) <- astro.name
    dobs
}

AstroDiff.array <- function(obs1,obs2){
###obs1, obs2: ra[deg], dec[deg], parallax [mas], pmra [mas/yr], pmdec [mas/yr], rv [km/s]
    astro.name <- c('ra','dec','parallax','pmra','pmdec','radial_velocity')
    o1 <- as.numeric(obs1[,astro.name])
    o2 <- as.numeric(obs2[,astro.name])
    dobs <- o2-o1
    dobs[,1:2] <- dobs[,1:2]*3.6e6#deg to mas
    dobs[,1] <- dobs[,1]*cos(o2[,2]*pi/180)
    colnames(dobs) <- astro.name
    dobs
}
####model of astrometry
astrometry.kepler <- function(pars.kep,tt=NULL,bases=rep('natural',10),Pmin=0){
    if(is.null(tt)){
        sim <- FALSE
    }else{
        sim <- TRUE
    }
    tmp <- list()
    if(length(out$astrometry)>0){
        tmp$barycenter <- astrometry.bary(pars.kep=pars.kep,tt=tt,bases=bases)
    }
    if(any(names(out)=='rel')){
        tmp$rel <- astrometry.rel(pars.kep,tt=tt,bases=bases)$rel
    }
    if(length(out$data.epoch)>0){
        for(jj in 1:length(out$ins.epoch)){
            i <-  out$ins.epoch[jj]
            if(!sim) tt <- out$data.epoch[[i]][,'BJD']
            if(grepl('hip',i)){
                band <- 'Hp'
            }else{
                band <- 'G'
            }
            tmp$epoch[[i]] <- astrometry.epoch(pars.kep,tt=tt,bases=bases,band=band,comp=out$comp.epoch[jj])$epoch
####calculate the dra and ddec relative to the Gaia DR3 reference epoch
            if(!grepl('hip|TYC',i)){
                if(!sim){
                    plx.ra <- out$pf[out$ind.all$epoch[[i]],'ra']
                    plx.dec <- out$pf[out$ind.all$epoch[[i]],'dec']
                }else{
                    plx.ra <- out$pf.sim[,'ra']
                    plx.dec <- out$pf.sim[,'dec']
                }
                dt <- (tt-out$data.epoch[[i]][1,1])/365.25#year
                dra <- tmp$barycenter[out$iref,'pmra']*dt+plx.ra*tmp$barycenter[out$iref,'parallax']+pars.kep[paste0('dra_',i)]
                ddec <- tmp$barycenter[out$iref,'pmdec']*dt+plx.dec*tmp$barycenter[out$iref,'parallax']+pars.kep[paste0('ddec_',i)]
                tmp$epoch[[i]] <- data.frame(dra=tmp$epoch[[i]][,'dra']+dra,ddec=tmp$epoch[[i]][,'ddec']+ddec)
            }
        }
    }

    if(length(out$gost)>0){
        reflex <- astrometry.epoch(pars.kep,tt=out$gost[,'BJD'],bases=bases,band='G')$epoch
        obs0 <- unlist(tmp$barycenter[out$iref,])
#        bary <- obs.lin.prop(obs0,t=out$gost[,'BJD']-out$astrometry[out$iref,'ref_epoch'])
        bary <- obs.lin.prop(obs0,t=out$gost[,'BJD']-out$astrometry[out$iref,'ref_epoch'],PA=FALSE)
        dec <- bary[,'dec']+reflex[,'ddec']/3.6e6#deg
        dra <- (bary[,'ra']-out$astrometry[out$iref,'ra'])*cos(dec/180*pi)*3.6e6+reflex[,'dra']#mas
        ddec <- (dec-out$astrometry[out$iref,'dec'])*3.6e6#mas
        gabs <- dra*sin(out$gost[,'psi'])+ddec*cos(out$gost[,'psi'])+(bary[,'parallax']+reflex[,'dplx'])*out$gost[,'parf']#mas
#        gabs <- dra*sin(out$gost[,'psi'])+ddec*cos(out$gost[,'psi'])+bary[,'parallax']*out$gost[,'parf']
        dabs <- list()
        binarycats <- cats <- c()
        for(k in 1:length(out$cats)){
            data <- data.frame(gabs=gabs[out$cat.ind[[k]]],a1=out$a1[[k]],a2=out$a2[[k]],a3=out$a3[[k]],a4=out$a4[[k]],a5=out$a5[[k]])
            if(out$cats[k]=='GDR1'){
                ast <- lm(gabs~a1+a2+0,data=data)$coefficients
                res <- out$astro.gost[k,1:2]-ast
                res <- c(res,rep(0,3))
            }else if(out$cats[k]=='TYC'){
                t <- out$astrometry[out$ityc,'ref_epoch']
                reflex1 <- astrometry.epoch(pars.kep,tt=t,bases=bases,band='H')$epoch
                bary1 <- obs.lin.prop(obs0,t=t-out$astrometry[out$iref,'ref_epoch'])
                dastro <- AstroDiff(bary1[1,],out$astrometry[out$ityc,])
                res <- dastro[c('ra','dec')]
                if(is.na(out$astrometry[out$ityc,'parallax'])){
                    res <- res-reflex1[,c('ra','dec')]
                }
                res <- c(res,rep(0,3))
            }else{
                fit <- lm(gabs~a1+a2+a3+a4+a5+0,data=data)
                ast <- fit$coefficients
                dabs[[out$cats[k]]] <- fit$residuals
                res <- out$astro.gost[k,]-ast
            }
            cats <- rbind(cats,res)
        }
        tmp$cats <- cats
        tmp$dabs <- dabs
        if(length(out[[binary]]$gost)>0){
            Np <- length(grep('^per',names(pars.kep)))
            pars.kep2 <- pars.kep
            if(Np>1){
                ind.rm <- 8:(Np*7)
                pars.kep2 <- pars.kep[-ind.rm]
            }
            mb <- k2m(pars.kep2['K1'],exp(pars.kep2['per1']),pars.kep2['e1'],pars.kep2['Mstar'],Inc=pars.kep2['Inc1'])$ms
            mratio <- mb/pars.kep2['Mstar']
            reflex2 <- astrometry.epoch(pars.kep2,tt=out$gost[,'BJD'],bases=bases,band='G')$epoch
            reflex <- -reflex2/mratio
            dec <- bary[,'dec']+reflex[,'ddec']/3.6e6#deg
            dra <- (bary[,'ra']-out[[binary]]$astrometry[out[[binary]]$iref,'ra'])*cos(dec/180*pi)*3.6e6+reflex[,'dra']#mas
            ddec <- (dec-out[[binary]]$astrometry[out[[binary]]$iref,'dec'])*3.6e6#mas
            gabs <- dra*sin(out[[binary]]$gost[,'psi'])+ddec*cos(out[[binary]]$gost[,'psi'])+(bary[,'parallax']+reflex[,'dplx'])*out[[binary]]$gost[,'parf']#mas
            for(k in 1:length(out[[binary]]$cats)){
                data <- data.frame(gabs=gabs[out[[binary]]$cat.ind[[k]]],a1=out[[binary]]$a1[[k]],a2=out[[binary]]$a2[[k]],a3=out[[binary]]$a3[[k]],a4=out[[binary]]$a4[[k]],a5=out[[binary]]$a5[[k]])
                if(out[[binary]]$cats[k]=='GDR1'){
                    ast <- lm(gabs~a1+a2+0,data=data)$coefficients
                    res <- out[[binary]]$astro.gost[k,1:2]-ast
                    res <- c(res,rep(0,3))
                }else if(out[[binary]]$cats[k]=='TYC'){
                    t <- out[[binary]]$astrometry[out[[binary]]$ityc,'ref_epoch']
                    reflex1 <- astrometry.epoch(pars.kep2,tt=t,bases=bases,band='H')$epoch
                    bary1 <- obs.lin.prop(obs0,t=t-out[[binary]]$astrometry[out[[binary]]$iref,'ref_epoch'])
                    dastro <- AstroDiff(bary1[1,],out[[binary]]$astrometry[out[[binary]]$ityc,])
                    res <- dastro[c('ra','dec')]
                    if(is.na(out[[binary]]$astrometry[out[[binary]]$ityc,'parallax'])){
                        res <- res-reflex1[,c('ra','dec')]
                    }
                    res <- c(res,rep(0,3))
                }else{
                    ast <- lm(gabs~a1+a2+a3+a4+a5+0,data=data)$coefficients
                    res <- out[[binary]]$astro.gost[k,]-ast
                }
                binarycats <- rbind(binarycats,res)
            }
        }
        tmp$binarycats <- binarycats
    }
    tmp
}

####functions for mcmc algorithms for detection of keplerian signal in RV data
####Keplerian model with 1 planet
RV.kepler <- function(pars.kep,tt=NULL,prior.kep=prior.type,injection=FALSE,kep.only=FALSE,noise.only=FALSE,bases=rep('natural',10),rv.pm=0,nqp=NULL){
    if(is.null(tt)){
        sim.kep <- FALSE
    }else{
        sim.kep <- TRUE
    }
    binarity <- FALSE
    if(length(out$data.binary)>0) binarity <- TRUE
    if(!sim.kep){
#        tt <- trv.all
        tt <- out$tiall[,1]
    }
#    cat('length(tt)=',length(tt),'\n')
#    if(!is.null(trv.all)){
#        if(max(tt)<24e5) tt <- tt+24e5
#        if(tmin<24e5) tmin <- tmin+24e5
#    }
#
    nqp.all <- nqp
    Np.kep <- length(grep('^dP|^per',names(pars.kep)))
#    tt <- tt-tmin
    Eastro <- list()
    Eb <- rvg <- rvs <- c()
    rvm <- rvc <- list()
    drvT <- rvgP <- rvsP <- rep(0,length(tt))
    if(any(names(pars.kep)=='Mstar')){
        Mstar <- pars.kep['Mstar']
    }else{
        Mstar <- out$Mstar
    }
    EC <- rvgC <- rvsC <- tauC <- drvC <- rvC <- 0
    tauT <- tau <- tt
    if(Np.kep>0){
        Ein <- list()
        pp <- extract.par(pars.kep,bases=bases)
###add relativistic effects
        if(binarity){
            if(!sim.kep){
                tauT <- tt
                indC <- length(tt)+1:nrow(out$data.binary)
                tauC <- out$data.binary[,1]
                tau <- c(tt,tauC)
            }else{
                indC <- 1:length(tt)
                tau <- tauC <- tt
            }
        }else{
            indC <- NULL
            tauT <- tauC <- tau <- tt
        }

        as <- Mps <- c()
        if(length(pp$K)>0){
            for(h in 1:Np.kep){
                if(pp$K[h]>0){
                    mp <- k2m(pp$K[h],pp$P[h],pp$e[h],Ms=Mstar,Inc=pp$Inc[h],more=TRUE)
                    Mps <- c(Mps,mp$ms)
                    as <- c(as,mp$a)
                }else{
                    Mps <- c(Mps,0)
                    as <- c(as,0)
                }
            }
        }
        imax <- which.max(pp$P)
        rvP <- rep(0,length(tt))
        tt0 <- c()
###initial eccentric anomaly without Roemer correction
#        if(out$Nrv>0){
        if(length(tt)>0){
            indP <- 1:length(tt)
#            cat('length(tau)=',length(tau),'\n')
            if(out$relativity & out$Nrv>0){
                tmp <- RV.relativity(pars.kep,t=tau,indP=indP,indC=indC,Mstar=Mstar,Mps=Mps,as=as,pp=pp,sim.kep=sim.kep)
                if(nepoch>0){
                    indP <- 1:length(tt)
                    for(kk in 1:Np.kep){
                        Eastro[[kk]] <- tmp$EP[[kk]][1:nepoch]
                        Ein[[kk]] <- tmp$EP[[kk]][indP]
                    }
                }else{
                    indP <- 1:length(tt)
                    Ein <- tmp$EP
                }
                rvgP <- tmp$rvgT
                rvsP <- tmp$rvsT
                drvT <- tmp$drvT
                rvP <- rvP+rvgP+rvsP+tmp$drvT#account for relativistic and astroemtric effect
#cat('GR+SR:',head(tmp$rvgs),'\n')
                tau <- tmp$tau
                tauT <- tmp$tauT
                if(binarity){
                    EC <- tmp$EC
                    tauC <- tmp$tauC
                    rvgC <- tmp$rvgC
                    rvsC <- tmp$rvsC
                    drvC <- tmp$drvC
                    rvC <- rvC+tmp$rvgsC+tmp$drvC
                }
            }else{
                for(h in 1:Np.kep){
                    m <- (pp$Mo[h]+2*pi*(tt-out$tref)/pp$P[h])%%(2*pi)#mean anomaly
                    if(prior.kep!='e0'){
                        Ein[[h]] <- kep.mt2(m,pp$e[h])%%(2*pi)
                    }else{
                        Ein[[h]] <- m
                    }
                }
                if(!is.null(out$timing)){
                    tauT <- delay(pars.kep,t=tt,E=Ein,Mstar=Mstar)$tau
                }
                if(binarity){
                    EC <- (pp$Mo[imax]+2*pi*(tauC-out$tref)/pp$P[h])%%(2*pi)
                }
            }
###
            if(out$Nrv>0){
                for(h in 1:Np.kep){
                    if(prior.kep!='e0'){
                        T <- 2*atan(sqrt((1+pp$e[h])/(1-pp$e[h]))*tan(Ein[[h]]/2))
                        rv <- pp$K[h]*(cos(pp$omega[h]+T)+pp$e[h]*cos(pp$omega[h]))
                    }else{
                        rv <- pp$K[h]*cos(Ein[[h]])
                    }
                    if(any(is.na(rv))) stop('rv is na!\n')
                    rvP <- rvP+rv
                }
                if(binarity){
                    T <- 2*atan(sqrt((1+pp$e[imax])/(1-pp$e[imax]))*tan(EC/2))
                    rv <- pp$K[imax]*(cos(pp$omega[imax]+T)+pp$e[imax]*cos(pp$omega[imax]))
                    rvC <- rvC-rv*Mstar/Mps[imax]
                }
            }
        }
        if(out$Nrvc>0){
            for(i in 1:Np.kep){
                if(any(i==out$Irvc)){
                    if(pp$K[i]>0){
###planet veloctiy
                        eta <- (Mstar+Mps[i])/Mps[i]
                        if(!sim.kep){
                            m <- (pp$Mo[i]+2*pi*(out$rvc[[paste0('p',i)]][,1]-out$tref)/pp$P[i])%%(2*pi)#mean anomaly
                        }else{
                            m <- (pp$Mo[i]+2*pi*(tt-out$tref)/pp$P[i])%%(2*pi)
                        }
                        if(prior.kep!='e0'){
                            E <- kep.mt2(m,pp$e[i])%%(2*pi)
                            T <- 2*atan(sqrt((1+pp$e[i])/(1-pp$e[i]))*tan(E/2))
                            omegac <- pp$omega[i]+pi
                            rv <- pp$K[i]*eta*(cos(omegac+T)+pp$e[i]*cos(omegac))
                        }else{
                            rv <- pp$K[i]*eta*cos(m+pi)
                        }

                        rvc[[paste0('p',i)]] <- rv/1e3
#                        cat('head(rvc)=',head(rvc[[paste0('p',i)]]),'km/s\n')
                        if(any(grepl('bc',names(pp)))) rvc[[paste0('p',i)]] <- rvc[[paste0('p',i)]]+pp[[paste0('bc',i)]]#km/s
                        if(any(names(pars.kep)=='rvb')) rvc[[paste0('p',i)]] <- rvc[[paste0('p',i)]]+pars.kep['rvb']#km/s
#
#                        cat('2head(rvc)=',head(rvc[[paste0('p',i)]]),'km/s\n')
#                        cat('rvc=',unlist(rvc),'km/s\n')
###add moon perturbation
                        if(out$Nm>0 & out$moon[[i]][1]>0){
                            inds <- (out$moon[[i]][2]:out$moon[[i]][3])+Np.kep
###moon mass
                                        #                                eta <- (msini$ms+Mp)/Mp
                            for(k3 in inds){
                                if(!sim.kep){
                                    m <- (pp$Mo[k3]+2*pi*(out$rvc[[paste0('p',i)]][,1]-out$tref)/pp$P[k3])%%(2*pi)#mean anomaly
                                }else{
                                    m <- (pp$Mo[k3]+2*pi*(tt-out$tref)/pp$P[k3])%%(2*pi)
                                }
                                if(prior.kep!='e0'){
                                    E <- kep.mt2(m,pp$e[k3])%%(2*pi)
                                    T <- 2*atan(sqrt((1+pp$e[k3])/(1-pp$e[k3]))*tan(E/2))
                                    rv <- pp$K[k3]*(cos(pp$omega[k3]+T)+pp$e[k3]*cos(pp$omega[k3]))
                                }else{
                                    rv <- pp$K[k3]*cos(m)
                                }
                                rvm[[i]] <- rv#moon-induced RV wobble of planet
                                rvc[[paste0('p',i)]] <- rvc[[paste0('p',i)]]+rv/1e3#km/s
                            }
                        }
                    }else{
                        rvc[[paste0('p',i)]] <- rep(0,length(tt))
                    }
                }
            }
        }
    }else{
        rvP <- rep(0,length(tt))
        if(binarity) rvC <- rep(0,nrow(out$data.binary))
    }
    ind.na <- which(is.na(tt))
    if(length(ind.na)>0){
        rvP[ind.na] <- NA
    }
###calculate noise component
    dVr.kep <- list()
    if(sim.kep){
        Nins <- 1
    }else{
        Nins <- length(out$ins.rv)
    }

    if(out$Nrv>0){
        for(j2 in 1:Nins){
###assign parameters and data
            ii <- out$ind.all$rv[[out$ins.rv[j2]]]
            if(is.null(nqp.all)){
                nqp <- out[[out$ins.rv[j2]]]$noise$nqp
            }else{
                nqp <- nqp.all[[out$ins.rv[j2]]]
            }
            if(!sim.kep){
#                dVr.kep[[out$ins.rv[j2]]] <- rvP[out[[out$ins.rv[j2]]]$index]
                dVr.kep[[out$ins.rv[j2]]] <- rvP[ii]
            }else{
                dVr.kep <- rvP
            }
            if(!kep.only){
                ind <- grep('^a',names(pars.kep))
                if(length(ind)>0){
                    at <- pars.kep[ind]
                }else{
                    at <- 0
                }
                b <- 0

                if(offset) b <- pars.kep[paste0('b_',out$ins.rv[j2])]
                if(!sim.kep){
                    trend <- cal.trend(a=at,b=b,t=(tau[ii]-out$tref)/time.unit)
                    dVr.kep[[out$ins.rv[j2]]] <- dVr.kep[[out$ins.rv[j2]]]+trend
                }else{
                    trend <- cal.trend(a=at,b=b,t=(tsim-out$tref)/time.unit)
                    dVr.kep <- dVr.kep+trend
                }

###check whether PM-induced RV trend is subtracted; perspective acceleration
                if(grepl('VLC|LC|CFHT',out$ins.rv[j2]) & length(out$astrometry)>0){
                    t3 <- out[[out$ins.rv[j2]]]$RV[,1]-out$astrometry[out$iref,1]
                    obs0 <- unlist(out$astrometry[out$iref,c('ra','dec','parallax','pmra','pmdec','radial_velocity')])#reference astrometry
                    tmp <- obs.lin.prop(obs0,t3)
                    rv.pm <- (tmp[,'radial_velocity']-tmp[1,'radial_velocity'])*1e3#m/s
                    dVr.kep[[out$ins.rv[j2]]] <- dVr.kep[[out$ins.rv[j2]]]+rv.pm
                }

                if(nqp[1]>0){
                    ind <- grep(paste0('c'),1:nqp[1])
                    if(length(ind)>0){
                        cs <- pars.kep[ind]
                        proxies <-out[[j2]]$noise$proxy.opt
                        if(length(cs)>1){
                            tmp <- rowSums(t(cs*t(proxies)))
                        }else{
                            tmp <- cs*proxies
                        }
                        dVr.kep[[out$ins.rv[j2]]] <- dVr.kep[[out$ins.rv[j2]]] + tmp
                    }
                }
            }
        }
    }

#####binary companion RV
    if(binarity & offset & !sim.kep){
        Nins <- length(out$ins.binary)
        for(j in 1:Nins){
            ins <- out$ins.binary[j]
            rvC[out$ind.binary[[ins]]] <- rvC[out$ind.binary[[ins]]]+pars.kep[paste0('b_',ins)]
        }
    }
    return(list(tt=tt,E=Ein,rv=dVr.kep,rvg=rvgP,drvT=drvT,rvs=rvsP,rvsC=rvsC,rvgC=rvgC,drvC=drvC,tauC=tauC,rvc=rvc,rvm=rvm,tau=tau,tauT=tauT,rvC=rvC))
}

####red noise kernel
ker.red <- function(x,y,sigma.red,l,Pgp=1,lp=1,type='abs'){
    if(type=='abs'){
        return(sigma.red^2*exp(-abs(x-y)/l))
    }
    if(type=='sq'){
        return(sigma.red^2*exp(-(x-y)^2/(2*l^2)))
    }
    if(type=='qp'){
#        return(sigma.red^2*exp(-(x-y)^2/(2*l^2)-(sin(pi*(x-y)/Pgp))^2/(2*lp^2)))
        val <- sigma.red^2*exp(-(x-y)^2/(2*l^2)-2*(sin(pi*(x-y)/Pgp))^2/(lp^2))#Haywood
        return(val)
    }
}
ker.sho <- function(x,y,S0,logQ,logProt){
    tau <- abs(x-y)
    Q <- exp(logQ)
    w0 <- 2*pi/exp(logProt)
    eta <- abs(1-(4*Q^2)^(-1))^0.5
    A <- S0*w0*Q*exp(-w0*tau/(2*Q))
    if(Q>0 & Q<1/2){
        A*(cosh(eta*w0*tau)+sinh(eta*w0*tau)/(2*eta*Q))
    }else if(Q==1/2){
        A*(2*(1+w0*tau))
    }else if(Q>1/2){
        A*(cos(eta*w0*tau)+sin(eta*w0*tau)/(2*eta*Q))
    }
}
ker.celerite <- function(x,y,term){
    a <- term[,1]
    b <- term[,2]
    c <- term[,3]
    d <- term[,4]
    R <- length(a)
    J <- R/2
    N <- length(x)
    tau <- abs(x-y)
    a*exp(-c*tau)*cos(d*tau)+b*exp(-c*tau)*sin(d*tau)
}
sho.term <- function(S0,logQ,logProt){
#    w0 <- 2*pi/Prot
    Q <- exp(logQ)
    w0 <- 2*pi/exp(logProt)
    a <- b <- c <- d<- 0
    if(Q<0.5){
        f <- sqrt(1-4*Q^2)
        a <- 0.5*S0*w0*Q*c(1+1/f,1-1/f)
        c <- 0.5*w0/Q*c(1-f,1+f)
    }else{
        f <- sqrt(4.0*Q^2-1)
        a <- S0*w0*Q
        b <- S0*w0*Q/f
        c <- 0.5*w0/Q
        d <- 0.5*w0*f/Q
    }
    return(cbind(a,b,c,d))
}
##refer to cholesky.h in celerite python module
celerite <- function(t,dy,y,s,term){
    a <- term[,1]
    b <- term[,2]
    c <- term[,3]
    d <- term[,4]
###a, b, c, d are for complex terms
###ar, cr are for real terms
#    Jc <- length(a)
#    Jr <- length(ar)
#    J <- Jc+Jr
#    R <- 2*J-Jr
    R <- 2*length(a)
    J <- length(a)
    N <- length(t)
    phi <- Wp <- Up <- Vp <- array(NA,dim=c(N,R))
    A <- D <- array(0,dim=c(N,N))
    S <- array(NA,dim=c(N,R,R))
#    diag(D) <- s^2+sum(a)+sum(ar)
    js <- 1:J
    ns <- 1:N
    dt <- outer(d,t,'*')
    ct <- outer(c,t,'*')
    cdt <- cos(dt)
    sdt <- sin(dt)
    nct <- exp(-ct)
    pct <- exp(ct)
    diag(A) <- dy^2+s^2+sum(a)
####pre-conditioned variables
    Up[,2*js-1] <- a%*%cdt+b%*%sdt
    Up[,2*js] <- a%*%sdt-b%*%cdt
    Wp[,2*js-1] <- cdt
    Wp[,2*js] <- sdt
    ect <- outer(c,t[2:N]-t[1:(N-1)],'*')
    phi[2:N,2*js-1] <- phi[2:N,2*js] <- exp(-ect)
    phi[1,] <- 0
####calculate S, D and W
    S[1,,] <- 0
    D[1,1] <- A[1,1]
    Wp[1,] <- 1/D[1,1]*Wp[1,]
    for(n in 2:N){
        S[n,,] <- outer(phi[n-1,],phi[n-1,],'*')*(S[n-1,,]+D[n-1,n-1]*outer(Wp[n-1,],Wp[n-1,],'*'))
        D[n,n] <- A[n,n]-Up[n,]%*%(S[n,,]%*%Up[n,])
        Wp[n,] <- 1/D[n,n]*(Wp[n,]-Up[n,]%*%S[n,,])
    }
    lndetK <- sum(log(diag(D)))
    f <- array(NA,dim=c(N,R))
    z <- rep(NA,N)
    z[1] <- y[1]
    f[1,] <- 0
    yky <- y[1]^2/D[1,1]
    for(n in 2:N){
        f[n,] <- phi[n,]*(f[n-1,]+Wp[n-1,]*z[n-1])
        z[n] <- y[n]-sum(Up[n,]*f[n,])
        yky <- yky+z[n]^2/D[n,n]
    }
    return(list(lndetK=lndetK,yky=yky))
}
inv.cel <- function(t,dy,s,a,b,c,d,ar,cr,y){
    R <- length(a)
    N <- length(t)
    tmp <- mat.cel(t,dy,s,a,b,c,d,ar,cr)
    W <- tmp$W
    U <- tmp$U
    D <- tmp$D
    phi <- tmp$phi
    f <- array(NA,dim=c(N,R))
    z <- rep(NA,N)
    z[1] <- y[1]
    f[1,] <- 0
    out <- y[1]^2/D[1,1]
    for(n in 2:N){
        f[n,2:R] <- phi[n,2:R]*(f[n-1,]+W[n-1,]*z[n-1])
        z[n] <- y[n]-sum(U[n,]*f[n,])
        out <- out+z[n]^2/D[n,n]
    }
    return(out)
}
ker.rot <- function(B, C, L, tau, Prot){
    B/(2+C)*exp(-tau/L)*(cos(2*pi*tau/Prot)+(1+C))
}
sho.chol <- function(tau,S0,Q,w0){
    K <- ker.sho(tau,S0,Q,w0)
    Sn
}
###red noise model
rednoise.cov <- function(sigma.red,l,tt=trv,tol=1e-12,Pgp=1,lp=1,type='abs'){
#method 1: easy but not so many choices of kernels
#    ker <- laplacedot(1/l)
#    cov.red <- kernelMatrix(ker,tt)#+diag(eRV^2)#to guarantee that it is positive definite
#method 2: more time-consuming
#    cov.red <- outer(1:length(RV),1:length(RV),Vectorize(function(u,v) ker.red(sigma.red,l,tt[u],tt[v],Pgp=Pgp,lp=lp,type=gp.type)))
#method 3:
    SE <- function(x,y) ker.red(x,y,sigma.red=sigma.red,l=l,Pgp=Pgp,lp=lp,type=gp.type)
    cov.red <- outer(tt, tt, SE)
    return(cov.red)
}
##arma noise model
arma <- function(t,ymodel,ydata,pars,p,q,set='',ARMAtype='abs'){
    dv.arma <- 0
    ##AR(p.ar) model
    if(p>0){
        dVr.ar <- 0
        phi <- pars[paste0('phi',1:p,'_',set)]
        alpha <- pars[paste0('alpha_',set)]
        for(j in 1:p){
            ar <- phi[j]*exp((t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(alpha))*ydata[-(length(t)+1-(1:j))]
            ar1 <- phi[j]*exp(-abs(t[j+(1:j)]-t[(1:j)])/exp(alpha))*ydata[j+(1:j)]#backward propagation for initial j epochs
#            dVr.ar <- dVr.ar + c(rep(ar[1],j),ar)
            dVr.ar <- dVr.ar + c(ar1,ar)
        }
        dv.arma <- dv.arma + dVr.ar
    }
    ##MA(q) model
    if(q>0){
        dVr.ma <- 0
        w <- pars[paste0('w',1:q,'_',set)]
        beta <- pars[paste0('beta_',set)]
        res <- ydata-ymodel
        for(j in 1:q){
            ma <- w[j]*exp(-abs(t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(beta))*res[-(length(t)+1-(1:j))]
            ma1 <- w[j]*exp(-abs(t[j+(1:j)]-t[(1:j)])/exp(beta))*res[j+(1:j)]#backward propagation for initial j epochs
#            dVr.ma <- dVr.ma + c(rep(ma[1],j),ma)
            dVr.ma <- dVr.ma + c(ma1,ma)
        }
        dv.arma <- dv.arma + dVr.ma
    }
    return(dv.arma)
}
# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
 par(mar = c(3,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
     axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
     axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
          cex.axis=0.7)

 # Color Scale
 par(mar = c(3,2.5,2.5,2))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

 layout(1)
}
matrix.power <- function(mat,power){
    for(j1 in 1:ncol(mat)){
        ind <- which(mat[,j1]>1e-3)
        mat[ind,j1] <- mat[ind,j1]^power
    }
    return(mat)
}
#likelihood
loglikelihood <- function(pars,bases=rep('natural',10),RVonly=FALSE,nqp=NULL,prediction=FALSE,indrel=NULL,verbose=FALSE){
    if(verbose){
        cat('\n')
        for(i in 1:length(pars)) cat(names(pars)[i],' ',pars[i],'\n')
    }
    if(is.null(indrel)){
        indrel <- out$indrel
        typerel <- out$typerel
    }
    pp <- extract.par(pars,bases=bases)
    nqp.all <- nqp
    res.all <- list()
    if(out$Nrv>0 | out$Nrvc>0 | length(out$tiall)>0){
        rvp  <-  RV.kepler(pars.kep=pars,bases=bases,nqp=nqp)
        RV.kep <- rvp$rv
    }
    Np <- length(grep('omega',names(pars)))
    epoch <- barycenter <- ref <- rel <- hg <- NULL
#    if((out$Nastro>0 | any(names(out)=='rel') | length(out$data.ref)>0 | length(out$data.epoch)>0 | length(out$gost)>0) & !RVonly){
    if((out$Nastro>0 | any(names(out)=='rel') | length(out$data.epoch)>0 | length(out$gost)>0) & !RVonly){
        astro <- astrometry.kepler(pars.kep=pars,bases=bases,Pmin=Tmin)
        barycenter <- astro$barycenter
        rel <- astro$rel
        ref <- astro$ref
        epoch <- astro$epoch
    }
##    if(class(RV.kep)=='try-error') RV.kep <- RV.kepler(pars.kep=startvalue)
    logLike <- 0
#dev.new()
#par(mfrow=c(4,4))
    ins <- out$ins.rv
    if(out$Nrv>0){
        for(k3 in 1:length(ins)){
            trv <- out[[out$ins.rv[k3]]]$RV[,1]
            rv.data <- out[[out$ins.rv[k3]]]$RV[,2]
            erv <- out[[out$ins.rv[k3]]]$RV[,3]
            rv.model <- rv.kep <- RV.kep[[out$ins.rv[k3]]]
            if(is.null(nqp.all)){
                nqp <- out[[out$ins.rv[k3]]]$noise$nqp
            }else{
                nqp <- nqp.all[[out$ins.rv[k3]]]
            }

            s <- pars[paste0('J_',out$ins.rv[k3])]
            if(prediction) res.all[[out$ins.rv[k3]]][['kep']] <- rv.kep
            if(prediction) res.all[[out$ins.rv[k3]]][['res1']] <- rv.data-rv.kep
            if(nqp[2]>0 | nqp[3]>0){
                rv.arma <- arma(t=trv,ymodel=rv.kep,ydata=rv.data,pars=pars,set=out$ins.rv[k3],p=nqp[3],q=nqp[2])
                rv.model <- rv.model +rv.arma
                if(prediction) res.all[[out$ins.rv[k3]]][['red']] <- rv.arma
            }
            if(prediction) res.all[[out$ins.rv[k3]]][['res2']] <- rv.data-rv.model
            ll <- sum(dnorm(rv.data,mean=rv.model,sd=sqrt(erv^2+s^2),log=T))
            logLike <-  logLike +ll
            if(verbose) cat('ll for RV=',ll,'\n')
        }
    }
#    save(list=ls(all=TRUE),file='test1.Robj')
    ll.astro <- 0
    if(any(names(out)=='rel')){
        for(j in 1:Np){
            if(length(indrel[[j]])>0){
                n <- paste0('p',j)
                inss <- names(out$rel[[n]])
                inss <- inss[inss!='tot']
                for(i in inss){
                    ll <- 0
                    if(typerel[[n]][[i]]=='rd'){
                        res <- cbind(out$rel[[n]][[i]][,'dra']-rel$dra[[n]][[i]],out$rel[[n]][[i]][,'ddec']-rel$ddec[[n]][[i]])
                        erel <- cbind(out$rel[[n]][[i]][,'era'],out$rel[[n]][[i]][,'edec'])
                        if(any(colnames(out$rel[[n]][[i]])=='dpmra')){
                            res <- cbind(res,out$rel[[n]][[i]][,'dpmra']-rel$dpmra[[n]][[i]],out$rel[[n]][[i]][,'dpmdec']-rel$dpmdec[[n]][[i]])
                            erel <- cbind(erel,out$rel[[n]][[i]][,'epmra'],out$rel[[n]][[i]][,'epmdec'])
                        }
                        if(any(names(out)=='cov')){
                            cov.rel <- out$cov[[n]][[i]]
                        }
                        ##                        cov <- matrix(c(erel[k,1]^2,out$rel[[n]][[i]][k,'cov'],out$rel[[n]][[i]][k,'cov'],erel[k,2]^2),nrow=2)
                    }else{
                        ps <- rd2ps(rel$dra[[n]][[i]],0,rel$ddec[[n]][[i]],0)#model
                        res <- cbind(out$rel[[n]][[i]][,'sep']-ps[,'sep'],out$rel[[n]][[i]][,'pa']-ps[,'pa'])#residual
                        erel <- cbind(out$rel[[n]][[i]][,'esep'],out$rel[[n]][[i]][,'epa'])#error
#                        cov <- list(diag(erel^2))
                    }
                    s <- 1
                    nn <- paste0('J_',n,'.',i)
                    if(any(names(pars)==nn)) s <- pars[nn]
#                    if(any(names(pars)==nn)) s <- exp()
                    erel <- erel*s
                    ##cat('head(erel)=',head(erel),'\n')
#                    if(length(ind0)>0)
#                    if(length(ind1)>0){
#                        for(k in ind1){
                    if(!exists('cov.rel')){
#                        ll <- ll+sum(dnorm(res[ind0,], mean=0, sd=erel[ind0,], log=TRUE))
                        ll <- ll+sum(dnorm(res, mean=0, sd=erel, log=TRUE))
                    }else{
                        for(k in 1:nrow(out$rel[[n]][[i]])){
                            ll <- ll+mvtnorm::dmvnorm(res[k,], mean=rep(0,length(res[k,])), sigma=cov.rel[[k]]*s, log=TRUE)
                        }
                    }
#
                    if(prediction) res.all[[paste0('rel',j)]] <- res
                    logLike <- logLike+ll
                    if(verbose) cat('ll for relative astrometry=',ll,'\n')
                }
            }
        }
    }

    if(!is.null(epoch)){
        ll <- 0
###reflex motion induced position change
        for(jj in 1:length(out$ins.epoch)){
            i  <- out$ins.epoch[jj]
            c  <- out$comp.epoch[jj]
            dpmdec <- dpmra <- dplx <- 0
            n1 <- paste0('J_',i)
            n2 <- paste0('logJ_',i)
            s <- 1
            if(any(names(pars)==n1)){
                s <- pars[paste0('J_',i)]
            }
            if(any(names(pars)==n2)){
                s <- 1+exp(pars[paste0('logJ_',i)])
            }
            data.epoch <- out$data.epoch[[i]]
##contribution of reflex motion to target astrometry
            dra <- epoch[[i]][,'dra']
            ddec <- epoch[[i]][,'ddec']
###since plx, pmra and pmdec are parameters *fixed* at the reference epoch, we should not consider the "time-varying" contribution from reflex motion
#            dplx <- epoch[[i]][,'dplx']
#            dpmra <- epoch[[i]][,'dpmra']
#            dpmdec <- epoch[[i]][,'dpmdec']
            if(FALSE){
                cat('dra:',dra,'mas\n')
                cat('ddec:',ddec,'mas\n')
                cat('dplx:',dplx,'mas\n')
                cat('dpmra:',dpmra,'mas/yr\n')
                cat('dpmdec:',dpmdec,'mas/yr\n')
            }
            if(grepl('hip',i)){
                dastro <- AstroDiff(out$astrometry[out$ihip,],barycenter[out$ihip,])
##contribution of barycenter-offset to target astrometry
                dra <- dra+dastro['ra']
                ddec <- ddec+dastro['dec']
                dplx <- dplx+dastro['parallax']
                dpmra <- dpmra+dastro['pmra']
                dpmdec <- dpmdec+dastro['pmdec']
                if(i=='hip1'){
                    dabs.new <- data.epoch[,'pv_pra']*dra+data.epoch[,'pv_pdec']*ddec+data.epoch[,'pv_pplx']*dplx+data.epoch[,'pv_ppmra']*dpmra*data.epoch[,'dt']+data.epoch[,'pv_ppmdec']*dpmdec*data.epoch[,'dt']#model of absissa residual
                    res <- data.epoch[,'dv']-dabs.new
                    if(is.null(out$cov.epoch[[i]])){
                        ll <- ll+sum(dnorm(res, mean=0, sd=data.epoch[,'ev'], log=TRUE))
                    }else{
                        ll <- ll+mvtnorm::dmvnorm(res,mean=rep(0,length(res)), sigma=out$cov.epoch[[i]]*s, log=TRUE)
                    }
                }
                if(i=='hip2'){
                    dabs.new <- data.epoch[,'CPSI']*(dra+dpmra*data.epoch[,'EPOCH'])+data.epoch[,'SPSI']*(ddec+dpmdec*data.epoch[,'EPOCH'])+data.epoch[,'PARF']*dplx
                    res <- data.epoch[,'RES']-dabs.new
                    ll <- ll+sum(dnorm(res, mean=0, sd=sqrt(data.epoch[,'SRES']^2+s^2), log=TRUE))
                }
            }else{
                dra.res <- out$data.epoch[[i]][,'dra']-dra
                ddec.res <- out$data.epoch[[i]][,'ddec']-ddec
                res <- c(dra.res,ddec.res)
                ll.epoch <- sum(dnorm(res,mean=0,sd=sqrt(c(out$data.epoch[[i]][,'era'],out$data.epoch[[i]][,'edec'])^2+s^2),log=TRUE))
#                cat('ll for ',i,'epoch:',ll.epoch,'\n')
                ll <- ll+ll.epoch
            }
####https://aas.aanda.org/articles/aas/pdf/1998/10/ds1401.pdf
####https://aas.aanda.org/articles/aas/pdf/2000/13/ds1810.pdf
####https://www.cosmos.esa.int/documents/532822/552851/vol3_all.pdf/dca04df4-dc6f-4755-95f2-b1217e539926
            if(prediction){
                res.all[[paste0('epoch_',i)]] <- res
            }
            logLike <- logLike+ll
            if(verbose) cat('ll for',i,'epoch astrometry=',ll,'\n')
        }
    }
####Gaia GOST fit
    if(!is.null(out$gost)){
        ll.gost <- 0
        s <- 0
        nast <- length(out$astro.index)
        for(k in 1:nast){
            j <- out$astro.index[k]
            s <- 1
            if(grepl('GDR',toupper(out$cats[k]))){
                if(any(names(pars)=='logJ_gaia')){
                    s <- 1+exp(pars['logJ_gaia'])
                }
                if(any(names(pars)=='J_gaia')){
                    s <- pars['J_gaia']
                }
            }
#            cat('s=',s,'\n')
            if(grepl('TYC',toupper(out$cats[k]))) s <- pars['J_TYC']
            if(grepl('HIP',toupper(out$cats[k]))) s <- pars['J_hip']
            if(prediction){
                res.all[[out$cats[k]]] <- astro$cats[k,]
            }
            if(out$cats[k]=='TYC'){
                ll <- mvtnorm::dmvnorm(astro$cats[k,1:2], mean=rep(0,2), sigma=out$cov.astro[1:2,1:2,j]*s, log=TRUE)
            }else{
                ll <- mvtnorm::dmvnorm(astro$cats[k,], mean=rep(0,5), sigma=out$cov.astro[1:5,1:5,j]*s, log=TRUE)
            }
            ll.gost <- ll.gost+ll
        }
        logLike <- logLike+ll.gost
        if(verbose) cat('ll for gost astrometry=',ll.gost,'\n')
        ll.ruwe <- 0
        if(ruweDR==3){
            uwe3 <- calc.uwe(astro$dabs$GDR3,sfov3,Nfov3,Nbin,Npar=5)
            ll.ruwe <- ll.ruwe+dnorm(log(uwe3),mean=log(ruwe3),sd=0.14,log=TRUE)
        }
        if(ruweDR==2){
            uwe2 <- calc.uwe(astro$dabs$GDR2,sfov2,Nfov2,Nbin,Npar=5)
            ll.ruwe <- ll.ruwe+dnorm(log(uwe2),mean=log(ruwe3),sd=0.25,log=TRUE)
        }
        logLike <- logLike+ll.ruwe
        if(verbose) cat('ll for RUWE astrometry=',ll.ruwe,'\n')
    }

    if(!is.null(out[[binary]]$gost)){
        ll.gost <- 0
        s <- 0
        nast <- length(out[[binary]]$astro.index)
        for(k in 1:nast){
            j <- out[[binary]]$astro.index[k]
            s <- 1
            if(grepl('GDR',toupper(out[[binary]]$cats[k]))) s <- pars['J_gaia']
            if(grepl('TYC',toupper(out[[binary]]$cats[k]))) s <- pars['J_TYC']
            if(grepl('HIP',toupper(out[[binary]]$cats[k]))) s <- pars['J_hip']
            if(prediction){
                res.all[[out[[binary]]$cats[k]]] <- astro$binarycats[k,]
            }
            if(out[[binary]]$cats[k]=='TYC'){
                ll <- mvtnorm::dmvnorm(astro$binarycats[k,1:2], mean=rep(0,2), sigma=out[[binary]]$cov.astro[1:2,1:2,j]*s, log=TRUE)
            }else{
                ll <- mvtnorm::dmvnorm(astro$binarycats[k,], mean=rep(0,5), sigma=out[[binary]]$cov.astro[1:5,1:5,j]*s, log=TRUE)
            }
            ll.gost <- ll.gost+ll
        }
        logLike <- logLike+ll.gost
        if(verbose) cat('ll for gost astrometry for binary=',ll.gost,'\n')
    }

####timing data and model
    if(!is.null(out$timing)){
        if(verbose) ll0.time <- 0
        ll.time <- 0
        for(n1 in names(out$timing)){
            for(n2 in names(out$timing[[n1]])){
                index <- out$ind.all$timing[[paste0(n1,'_',n2)]]
                if(!grepl('pulsation|eb',n2)){
                    ts <- rvp$tt[index]
                    lte  <-  (rvp$tt[index]-rvp$tauT[index])*24*60#min
                    ip <- as.integer(gsub('p','',n1))
                    p <- pp$P[ip]#day
                    dt <- out$timing[[n1]][[n2]][,1]-out$ttv0
                    if(grepl('occult',n2)) dt <- dt+p/2
                    ttv <- (dt-round(dt/p)*p)*24*60
                    ettv <- out$timing[[n1]][[n2]][,2]*24*60
                    ll.time <- ll.time+sum(dnorm(ttv,lte+pars['t0'],ettv,log=TRUE))
                    if(verbose) ll0.time <- ll0.time+sum(dnorm(lm(ttv~ts,weights=1/ettv^2)$res,0,ettv,log=TRUE))
                }else{
                    n3 <- gsub('pulsation|eb','',n2)
#                    if(n3=='') n3 <- '0'
                    if(n3=='') n3 <- ''
                    dt <- rvp$tt[index]-rvp$tauT[index]#day
#                    cat('dt=',head(dt),'\n')
##baseline delay model is a parabola: dti=ti+gi*dt.yr+hi*dt.yr
#                    dt.yr <- (out$timing[[n1]][[n2]][,1]-out$timing[[n1]][[n2]][1,1])/365.25
                    dt.yr <- (out$timing[[n1]][[n2]][,1]-out$tref)/365.25
#                    save(list=ls(all=TRUE),file='test1.Robj')
                    if(grepl('pulsation',n2)){
                        dt.model <- dt*3600*24+pars[paste0('t',n3)]+pars[paste0('g',n3)]*dt.yr+pars[paste0('h',n3)]*dt.yr^2
                    }else if(grepl('eb',n2)){
                        dt.model <- dt+pars['t0']+pars['g0']*dt.yr+pars['h0']*dt.yr^2
                    }
                    ll1 <- sum(dnorm(out$timing[[n1]][[n2]][,'dt'],mean=dt.model,sd=out$timing[[n1]][[n2]][,'et'],log=TRUE))
                    ll.time <- ll.time+ll1
                }
            }
        }
        if(verbose){
            cat('ll.time=',ll.time,'\n')
            cat('ll.time-ll0.time=',ll.time-ll0.time,'\n')
        }
        logLike <- logLike+ll.time
    }

#####companion's RV
    if(out$Nrvc>0){
        for(j in 1:Np){
            if(any(j==out$Irvc)){
                n <- paste0('p',j)
                res <- out$rvc[[n]][,2]-rvp$rvc[[n]]
#                cat('out$rvc[[n]]=',out$rvc[[n]][,2],'\n')
#                cat('rvp$rvc[[n]]=',rvp$rvc[[n]],'\n')
#                cat('res=',res,'\n')
#                cat('ll=',ll,'\n')
                ll <- sum(dnorm(res, mean=0, sd=out$rvc[[n]][,3], log=TRUE))
                logLike <- logLike+ll
                if(verbose) cat('ll for relative RV=',ll,'\n')
            }
        }
    }

####binary
   if(length(out$data.binary)>0){
       Nins <- length(out$ins.binary)
       for(j in 1:Nins){
           ins <- out$ins.binary[j]
           res <- out$data.binary[out$ind.binary[[ins]],2]-rvp$rvC[out$ind.binary[[ins]]]
           jitter <- pars[paste0('s_',ins)]
           ll <- sum(dnorm(res, mean=0, sd=sqrt(out$data.binary[out$ind.binary[[ins]],3]^2+jitter^2), log=TRUE))
           logLike <- logLike+ll
           if(verbose) cat('ll for binary=',ll,'\n')
       }
    }
                                        #    cat('RV+abs+rel+rvc logLike=',logLike,'\n')
    if(!prediction){
        return(logLike)                                        #        return(list(ll=logLike,llastro=ll.astro))
    }else{
        return(list(loglike=logLike,res=res.all))
    }
}

cal.residual <- function(pars,bases=rep('natural',10),nqp=NULL){
    nqp.all <- nqp
    RV.all  <-  RV.kepler(pars.kep=pars,bases=bases,nqp=nqp)$rv
    RV.kep  <-  RV.kepler(pars.kep=pars,kep.only=TRUE,bases=bases,nqp=nqp)$rv
    trend <- arma1 <- all <- noise <- signal <- residual.all <- residual.sigtrend <- residual.noise <- residual.sig <- list()
    if(out$Nrv>0){
        for(k3 in 1:length(out$ins.rv)){
            trv <-out[[out$ins.rv[k3]]]$RV[,1]
            rv.data <- out[[out$ins.rv[k3]]]$RV[,2]
            erv <- out[[out$ins.rv[k3]]]$RV[,3]
            rv.all <- RV.all[[out$ins.rv[k3]]]
            rv.sig <- RV.kep[[out$ins.rv[k3]]]
            rv.trend <- rv.all-rv.sig
            if(is.null(nqp.all)){
                nqp <- out[[out$ins.rv[k3]]]$noise$nqp
            }else{
                nqp <- nqp.all[[out$ins.rv[k3]]]
            }
            if((nqp[2]>0 | nqp[3]>0)){
                rv.arma <- arma(t=trv,ymodel=rv.all,ydata=rv.data,pars=pars,set=out$ins.rv[k3],p=nqp[3],q=nqp[2])
                if(FALSE){
                    cat('head(rv.arma)=',head(rv.arma),'\n')
                    cat('head(rv.all)=',head(rv.all),'\n')
                    cat('head(rv.data)=',head(rv.data),'\n')
                    cat('head(rv.sig)=',head(rv.sig),'\n')
                    cat('head(rv.trend)=',head(rv.trend),'\n\n')
                }
                rv.all <- rv.all +rv.arma
            }else{
                rv.arma <-rep(0,length(trv))
            }
            res.all <- rv.data-rv.all
            res.sig <- rv.data-rv.sig
            res.sigtrend <- rv.data-rv.sig-rv.trend
            rv.noise <- rv.all-rv.sig
            res.noise <- rv.data-rv.noise
            residual.all[[out$ins.rv[k3]]] <- res.all
            residual.sig[[out$ins.rv[k3]]] <- res.sig
            residual.sigtrend[[out$ins.rv[k3]]] <- res.sigtrend
            residual.noise[[out$ins.rv[k3]]] <- res.noise
            signal[[out$ins.rv[k3]]] <- rv.sig
            noise[[out$ins.rv[k3]]] <- rv.noise
            all[[out$ins.rv[k3]]] <- rv.all
            arma1[[out$ins.rv[k3]]] <- rv.arma
            trend[[out$ins.rv[k3]]] <- rv.trend
        }
    }
    list(res.sig=residual.sig,res.sigtrend=residual.sigtrend,res.all=residual.all,res.noise=residual.noise,rv.signal=signal,rv.noise=noise,rv.all=all,rv.arma=arma1,rv.trend=trend)
}

tjar <- function(t=trv,x=Sindex[,1],phi.sab,alpha.sab,psi.sab,par.tj=1,symmetry=TJAR.sym){
    rv.tjar <- 0
    for(i0 in 1:length(phi.sab)){
        if(symmetry=='sym'){
            rv.tjar <- rv.tjar + par.tj*c(rep(0,i0),phi.sab[i0]*exp(alpha.sab/time.unit*(t[-(length(t)+1-(1:i0))]-t[-(1:i0)]))*x[-(length(t)+1-(1:i0))]) + par.tj*c(psi.sab[i0]*exp(-alpha.sab/time.unit*(t[-(1:i0)]-t[-(length(t)+1-(1:i0))]))*x[-(1:i0)],rep(0,i0))
        }else if(symmetry=='past'){
            rv.tjar <- rv.tjar + par.tj*c(rep(0,i0),phi.sab[i0]*exp(alpha.sab/time.unit*(t[-(length(t)+1-(1:i0))]-t[-(1:i0)]))*x[-(length(t)+1-(1:i0))])
        }else if(symmetry=='future'){
            rv.tjar <- rv.tjar + par.tj*c(phi.sab[i0]*exp(-alpha.sab/time.unit*(t[-(1:i0)]-t[-(length(t)+1-(1:i0))]))*x[-(1:i0)],rep(0,i0))
        }
    }
    return(rv.tjar)
}
prior.func <- function(pars,bases='natural',nqp=NULL){
    nqp.all  <- nqp
    logprior <- 0
    pp <- extract.par(pars,bases=bases)
#    Np <- length(grep('per',names(pars)))
    Np <- length(pp$P)
    Ntot <- Np+out$Nm
    logkep.prior <- 0
    if(Ntot>0){
        for(k in 1:Ntot){
            basis <- bases[k]
            ss <- ''
            if(k>Np) ss <- 'm'
            indI <- grep(paste0('^',ss,'Inc',k),names(pars))
            indap <- grep(paste0('^',ss,'apm',k),names(pars))
            indM <- grep('Mstar',names(pars))
            mprior <- 0
            if(length(indM)>0){
                Mstar <- pars['Mstar']
                if(priorf==1) mprior <- dnorm(Mstar,mean=out$Mstar,sd=out$eMstar,log=TRUE)
#                cat('Mstar=',Mstar,';mprior=',mprior,';out$Mstar',out$Mstar,';out$eMstar',out$eMstar,'\n')
            }else{
                Mstar <- out$Mstar
            }
            logkep.prior <- logkep.prior+mprior
#            cat('mprior=',mprior,'\n')
            if(basis=='natural'){
                indP <- grep(paste0('^',ss,'per',k),names(pars))
                indK <- grep(paste0('^',ss,'K',k),names(pars))
                inda <- grep(paste0('^',ss,'arc',k),names(pars))
                inde <- grep(paste0('^',ss,'e',k),names(pars))
                indMo <- grep(paste0('^',ss,'Mo',k),names(pars))
                indomega <- grep(paste0('^',ss,'omega',k),names(pars))
                e <- pars[inde]
                if(!out$NormMass){
                    logP.prior <- 0
                    if(length(indP)>0){
                        if(!resonance | k==1){
                            logP.prior <-  log(1/(par.max[indP]-par.min[indP]))
                        }else{
                                        #                        logP.prior <- dnorm(pars[indP],pars['per1']+log(laplace[k]),0.1,log=TRUE)
                            logP.prior <- dnorm(pars[indP],log(pp$P[1])+log(laplace[k]),0.2,log=TRUE)
                        }
                    }
                    if(target=='GaiaBH2'){
#                        logP.prior <- dnorm(exp(pars[indP]),1352,45,log=TRUE)
                    }
#                    cat('logP.prior=',logP.prior,'\n')
                    logkep.prior <- logkep.prior+logP.prior
                    if(length(indK)>0) logkep.prior <- logkep.prior + log(1/(par.max[indK]-par.min[indK]))
                    if(length(inda)>0) logkep.prior <- logkep.prior + log(1/(par.max[inda]-par.min[inda]))
                }else{
                    mp <- k2m(pars[indK],exp(pars[indP]),e,Mstar,Inc=pars[indI])$mj
                    lp <- dnorm(mp,out$mp[1,k],out$mp[2,k],log=TRUE)
                    logkep.prior <- logkep.prior + lp
                }
                if(Esd<1){
                    logkep.prior <- logkep.prior + log(2*dnorm(e,mean=0,sd=Esd))#normalized semi-Gaussian distribution
                }
                Mo <- pars[indMo]
                logkep.prior <- logkep.prior + log(1/(par.max[indMo]-par.min[indMo]))
                logkep.prior <- logkep.prior + log(1/(par.max[indomega]-par.min[indomega]))
                if(length(indap)>0){
                    logkep.prior <- logkep.prior + sum(log(1/(par.max[indap]-par.min[indap])))
                }
                if(length(indI)>0){
                    if(coplanar){
                        if(k==1){
                            logkep.prior <- logkep.prior + log(sin(pars[indI])/2)
                        }else{
#                            logp <- dnorm(pars[indI],pars['Inc1'],0.01,log=TRUE)
                            logp <- dnorm(pars[indI],pars['Inc1'],0.1,log=TRUE)
                            logkep.prior <- logkep.prior + logp
                        }
                    }else{
#                    if(target=='GaiaBH2'){
#                        logkep.prior <- dnorm(pars[indI],0.6085963,0.005934119,log=TRUE)
#                    }else{
                        if(target=='CPD-632495'){
#                            logkep.prior <- logkep.prior + dnorm(pars[indI],2.687807,0.05235988,log=TRUE)
                        }else if(target=='TYC3588-11669-1'){
#                            logkep.prior <- logkep.prior + dnorm(pars[indI],0.96,0.1,log=TRUE)
                        }else{
                            logkep.prior <- logkep.prior + log(sin(pars[indI])/2)
                        }
#                    }
                    }
                }
            }else{
                if(basis=='linear2'){
                    dP  <- pars[grep(paste0('^',ss,'dP',k),names(pars))]
                    logkep.prior <- logkep.prior+dnorm(dP,0,ePtransit[k]/alpha.dP,log=TRUE)
                    dTc  <- pars[grep('^',ss,'dTc',names(pars))]
                    logkep.prior <- logkep.prior+dnorm(dTc,0,eTc[k]/alpha.dTc,log=TRUE)
                }else if(basis=='linear1'){
                    indP <- grep(paste0('^',ss,'per',k),names(pars))
                    indTc <- grep(paste0('^',ss,'Tc',k),names(pars))
                    if(k<=Ntr){
                        P  <-  exp(pars[indP])
                        Tc1  <-  pars[indTc]
                        logkep.prior <- logkep.prior+ dnorm(P,mean=Ptransit[k],sd=ePtransit[k],log=TRUE)
                        logkep.prior <- logkep.prior+ dnorm(Tc1,mean=Tc[k],sd=eTc[k],log=TRUE)
                    }else{
                        logkep.prior <- logkep.prior+ log(1/(par.max[indP]-par.min[indP]))
                        logkep.prior <- logkep.prior+ log(1/(par.max[indTc]-par.min[indTc]))
                    }
                }
                indK <- grep(paste0('^',ss,'lnK',k),names(pars))
                indsw  <-  grep(paste0('^',ss,'sqresinw',k),names(pars))
                indcw  <-  grep(paste0('^',ss,'sqrecosw',k),names(pars))
                logkep.prior <- logkep.prior+ log(1/(par.max[indK]-par.min[indK]))
                logkep.prior <- logkep.prior+ log(1/(par.max[indsw]-par.min[indsw]))
                logkep.prior <- logkep.prior+ log(1/(par.max[indcw]-par.min[indcw]))
                if(length(indI)>0){
                    if(coplanar){
                        if(k==1){
                            logkep.prior <- logkep.prior + log(sin(pars[indI])/2)
                        }else{
                            logp <- dnorm(pars[indI],pars['Inc1'],0.01,log=TRUE)
                            logkep.prior <- logkep.prior + logp
                        }
                    }else{
                        logkep.prior <- logkep.prior + log(sin(pars[indI])/2)
                    }
                }
            }
        }
    }
#    if(target=='CPD-632495'){
    if(FALSE){
        m1 <- pars['Mstar']
        m2 <- k2m(pars['K1'],exp(pars['per1']),pars['e1'],m1,Inc=pars['Inc1'])$ms
        ls2au <- 0.0020039888
        mtot <- m1+m2
        asini <- ((exp(pars['per1'])/365.25)^2*mtot)^(1/3)*sin(pars['Inc1'])*m1/(m1+m2)
#        if(FALSE){
        if(TRUE){
        logkep.prior <- dnorm(exp(pars['per1']),1236.724526,6e-6,log=TRUE)#P prior
        logkep.prior <- logkep.prior+dnorm(pars['e1'],0.8698797,6e-8,log=TRUE)#e prior
        logkep.prior <- logkep.prior+dnorm(pars['omega1'],5.56175365,1.9e-6,log=TRUE)#omega prior
        logkep.prior <- logkep.prior+dnorm(pars['Omega1'],3.2986,0.0349,log=TRUE)#Omega prior
#        logkep.prior <- logkep.prior+dnorm(pars['Inc1'],2.6889,0.0537,log=TRUE)#Inc prior
        logkep.prior <- logkep.prior+dnorm(asini,1296.27448*ls2au,0.00014*ls2au,log=TRUE)#Inc prior
        logkep.prior <- logkep.prior+dnorm(pars['Mo1'],1.13754224683,1.2e-7,log=TRUE)#Mo prior
        }else{
        logkep.prior <- dnorm(exp(pars['per1']),1236.724526,6e-3,log=TRUE)#P prior
        logkep.prior <- logkep.prior+dnorm(pars['e1'],0.8698797,6e-3,log=TRUE)#e prior
        logkep.prior <- logkep.prior+dnorm(pars['omega1'],5.56175365,1.9e-2,log=TRUE)#omega prior
        logkep.prior <- logkep.prior+dnorm(pars['Omega1'],3.2986,0.0349,log=TRUE)#Omega prior
        logkep.prior <- logkep.prior+dnorm(asini,1296.27448*ls2au,0.014*ls2au,log=TRUE)#Inc prior
#        logkep.prior <- logkep.prior+dnorm(pars['Inc1'],2.6889,0.0537,log=TRUE)#Inc prior
        logkep.prior <- logkep.prior+dnorm(pars['Mo1'],1.13754224683,1.2e-2,log=TRUE)#Mo prior
        }
    }

#    if(target=='HR8799' & stability)
    if(stability & Np>1){
        pp <- extract.par(pars,bases=bases)
###The Stability of Multi-Planet Systems: https://ui.adsabs.harvard.edu/abs/1996Icar..119..261C/abstract
        if(any(names(pars)=='Mstar')){
            m1 <- pars['Mstar']
        }else{
            m1 <- out$Mstar
        }
        ind <- sort(pp$P,index.return=TRUE)$ix#decreasing=TRUE,
        m2s <- k2m(pp$K[ind],pp$P[ind],pp$e[ind],m1,Inc=pp$Inc[ind])$ms
        mt <- cumsum(c(m1,m2s))[-1]
        a1 <- ((pp$P[ind][-Np]/365.25)^2*mt[-Np])^(1/3)
        a2 <- ((pp$P[ind][-1]/365.25)^2*mt[-1])^(1/3)
        K <- ((m2s[-Np]+m2s[-1])/3)^(1/3)
        Rh <- K*((a1+a2)/2)
        Delta <- (a2-a1)/Rh
        e2 <- pp$e[-Np]^2+pp$e[-1]^2-2*pp$e[-Np]*pp$e[-1]*cos(pp$omega[-Np]-pp$omega[-1])
        i2 <- pp$Inc[-Np]^2+pp$Inc[-1]^2-2*pp$Inc[-Np]*pp$Inc[-1]*cos(pp$Omega[-Np]-pp$Omega[-1])
        Delta02 <- Delta^2-(4/3*K^2)*(e2+i2)
        sta.prior <- 0
#        ind <- which(Delta02<12)
#        if(length(ind)>0) sta.prior <- sta.prior+sum(dnorm(Delta02[ind],12,0.1,log=TRUE))
        ind <- which(Delta02<100)
        if(length(ind)>0) sta.prior <- sta.prior+sum(dnorm(Delta02[ind],100,1,log=TRUE))
        logkep.prior <- logkep.prior + sta.prior
    }
    logprior <- logprior + logkep.prior
####noise model
    if(out$Nrv>0){
        for(i1 in 1:length(out$ins.rv)){
            if(is.null(nqp.all)){
                nqp <- out[[out$ins.rv[i1]]]$noise$nqp
            }else{
                nqp <- nqp.all[[out$ins.rv[i1]]]
            }
            n <- nqp[1];q <- nqp[2];p <- nqp[3]
            phiprior = wprior = 1/(phi.max-phi.min)
            if(p>0){
                alpha.prior <- 1/(alpha.max-alpha.min)
                logprior <- logprior +p*log(phiprior)+log(alpha.prior)
            }
            if(q>0){
                beta.prior <- 1/(beta.max-beta.min)
                logprior <- logprior +q*log(wprior)+log(beta.prior)
            }
        }
    }
    lp <- 0
    if(priorf>0 & target=='UCAC4569-026385'){
#        lp <- lp+dnorm(out$astrometry[out$iref,'parallax']-pars['dplx'],0.68,0.05,log=TRUE)
        plx <- out$astrometry[out$iref,'parallax']-pars['dplx']
        if(priorf==1) lp <- lp+dnorm(plx,0.58,0.05,log=TRUE)
        if(priorf==2) lp <- lp+dnorm(plx,0.62,0.05,log=TRUE)
        if(priorf==3) lp <- lp+dnorm(plx,0.6,0.1,log=TRUE)
        out$Mstar <- pmfun.spec(plx)
        out$eMstar <- pemfun.spec(plx)
    }
    if(priorf>0 & target=='GaiaBH2'){
#        lp <- lp+dnorm(out$astrometry[out$iref,'parallax']-pars['dplx'],0.68,0.05,log=TRUE)
        plx <- out$astrometry[out$iref,'parallax']-pars['dplx']
        if(priorf==1) lp <- lp+dnorm(plx,0.81,0.1,log=TRUE)
        if(priorf==2) lp <- lp+dnorm(plx,0.89,0.05,log=TRUE)
        if(priorf==3) lp <- lp+dnorm(plx,0.85,0.1,log=TRUE)
    }
    if(target=='HD253754'){
        lp <- lp+dnorm(out$astrometry[out$iref,'parallax']-pars['dplx'],1.1,0.3,log=TRUE)
    }
    if(target=='CPD-632495'){
        bary <- astrometry.bary(pars.kep=pars,tt=2400000.5+55000,bases=bases)
        astro.pt <- hdms2deg(13,02,47.6426,-63,50,08.665)
###without constrained parallax
        dastro <- AstroDiff(c(ra=astro.pt[1],dec=astro.pt[2],parallax=3,pmra=-6.6,pmdec=-4.4,radial_velocity=0),bary[1,])
        era <- 11e-3/3600/24*360*3.6e6*cos(bary[,'dec']/180*pi)
        edec <- 8e-3*1e3#mas
###assign large uncertainty to parallax and RV
        lp <- lp+sum(dnorm(dastro,rep(0,6),c(era,edec,10,1.8,1.4,10),log=TRUE))
#        cat('lp=',lp,'\n')
    }
####Jitter prior
    if(length(out$rel.ins)>0){
        js <- unlist(sapply(out$rel.ins, function(i) pars[grep(paste0('J_.+',i),names(pars))]))
        lp <- lp + sum(dnorm(js,1,0.1,log=TRUE))
#        cat('prior for ',out$rel.ins,':',dnorm(js,1,0.1,log=TRUE),'\n')
    }
    if(length(out$astro.ins)>0){
        js <- unlist(sapply(out$astro.ins, function(i) pars[grep(paste0('J_',i),names(pars))]))
        lp <- lp + sum(dnorm(js,1,0.1,log=TRUE))
#        cat('prior for ',out$astro.ins,':',dnorm(js,1,0.1,log=TRUE),'\n')
    }
####total
    logprior <- logprior + lp
####for special cases
    if(target=='HD41004'){
        k <- m2K(m=0.4,P=exp(pars['per2']),e=pars['e2'],Ms=pars['Mstar'],Inc=pars['Inc2'],unit='ms')
        lp <- dnorm(pars['K2'],mean=k,sd=10,log=TRUE)+dnorm(pars['Mstar'],mean=0.7,sd=0.1,log=TRUE)
#        cat('lp=',lp,'\n')
        logprior <- logprior+lp
    }
    return(logprior)
}

#posterior distribution
posterior <- function(param,tem=1,bases='natural',RVonly=FALSE,verbose=FALSE){
    llike <- loglikelihood(param,bases=bases,RVonly=RVonly,verbose=verbose)
    pr <- prior.func(param,bases=bases)
    post <- llike*tem+pr
    return(list(loglike=llike,logprior=pr,post=post))
}

proposalfunction.simple <- function(param,cov.adapt,par.min,par.max){
    Ntt <- 1e5
    n <- 10
    param.new <- param
    for(k in 1:Ntt){
        par1 <- try(mvrnorm(n=1,mu=param,Sigma=cov.adapt,tol=tol1),TRUE)
        if(class(par1)[1]=='try-error'){
#            cat('ok1\n')
            par1 <- try(mvrnorm(n=1,mu=param,Sigma=nearPD(cov.adapt)$mat,tol=tol1),TRUE)
#        }else{
#            cat('ok2\n')
        }
        param.new <- par1
        param.new[grep('^e\\d',names(param.new))] <- abs(param.new[grep('^e\\d',names(param.new))])
#        cat('e=',param.new[grep('^e\\d',names(param.new))],'\n')
        ind1 <- grep('Mo\\d',names(param.new))
        ind2 <- grep('omega\\d',names(param.new))
        ind3 <- grep('Omega\\d',names(param.new))
        ind4 <- grep('Inc\\d',names(param.new))
        if(length(ind1)>0) param.new[ind1] <- param.new[ind1]%%(2*pi)
        if(length(ind2)>0) param.new[ind2] <- param.new[ind2]%%(2*pi)
        if(length(ind3)>0) param.new[ind3] <- param.new[ind3]%%(2*pi)
        if(length(ind4)>0) param.new[ind4] <- param.new[ind4]%%pi
#        if(length(ind4)>0) param.new[,ind4] <- param.new[,ind4]%%pi
####don't mod the transit epoch by period
        ind1 <- grep('sqresinw',names(param.new))
        ind2 <- grep('sqrecosw',names(param.new))

        sim.stop <- FALSE
        if(length(ind1)>0){
            ind <- c(ind1,ind2)
            if(all(param.new[-ind]>par.min[-ind] & param.new[-ind]<par.max[-ind]) & all((param.new[ind1]^2+param.new[ind2]^2)<1)){
                break()
            }
        }else{
            if(all(param.new>par.min & param.new<par.max)){
                break()
            }else if(FALSE){
#            }else{
                ind <- which(param.new<par.min | param.new>par.max)
                cat(names(param.new)[ind],'\n')
                cat('param.new=',param.new[ind],'\n')
                cat('par.min=',par.min[ind],'\n')
                cat('par.max=',par.max[ind],'\n')
            }
        }
    }
    e <- param.new[grep('^e\\d',names(param.new))]
    if(any(e<0) | k==Ntt){
        cat('Proposing time exceeding the maximum ',Ntt,' times or negative eccentricity!!\n')
#        save(list=ls(all=TRUE),file='test0.Robj')
        stop()
    }
    if(target=='CPD-632495'){
        m1 <- param.new['Mstar']
        ls2au <- 0.0020039888
        a2sini <- 1296.27448*ls2au
        Pyr <- exp(param.new['per1'])/365.25
#    a1 <- (K/sinI)^2/(4*pi^2)*(1-e^2)*P^(2/3)
##https://www.aanda.org/articles/aa/full_html/2011/01/aa15427-10/aa15427-10.html
        K2 <- 2*pi*a2sini/Pyr/sqrt(1-e^2)*4.74047e3
        q <- param.new['K1']/K2
        mtot <- m1*(1+q)
        a2 <- (Pyr^2*mtot)^(1/3)/(1+q)
        if(a2<a2sini){
###non-physical
            param.new <- param
        }else{
            i <- asin(a2sini/a2)
            if(i<pi/2) i <- pi-i
            param.new['Inc1'] <- i
        }
    }
    if(target=='TYC3588-11669-1'){
###from known parameters (m1,K1,a2sini,P) determine the unknown parameters including m2 and inclination
        m1 <- param.new['Mstar']
        ls2au <- 0.0020039888
        a2sini <- 856.26465*ls2au
        Pyr <- exp(param.new['per1'])/365.25
        K1 <- param.new['K1']

#    a1 <- (K/sinI)^2/(4*pi^2)*(1-e^2)*P^(2/3)
##https://www.aanda.org/articles/aa/full_html/2011/01/aa15427-10/aa15427-10.html
        K2 <- 2*pi*a2sini/Pyr/sqrt(1-e^2)*4.74047e3
        q <- K1/K2
        mtot <- m1*(1+q)
        a2 <- (Pyr^2*mtot)^(1/3)/(1+q)
        if(a2<a2sini){
            param.new <- param
        }else{
            i <- asin(a2sini/a2)
#            if(i<pi/2) i <- pi-i
            param.new['Inc1'] <- i
        }
    }
#    cat('param.new-param=',param.new-param,'\n')
    return(param.new)
}

#Metropolis algorithm####
#####Initial proposal function
proposalfunction <- function(param,cov.adapt){
    Ntt <- 1e4
    param.new <- param
    index <- 1:length(param)
    if(!is.null(par.fix)){
        inds <- match(par.fix,names(param))
        param <- param[-inds]
        index <- index[-inds]
        cov.adapt <- cov.adapt[-inds,][,-inds]
    }
    for(k in 1:Ntt){
        par1 <- try(mvrnorm(n=1,mu=param,Sigma=cov.adapt,tol=tol1,empirical=FALSE))#this could cause some non-zero values for fix value, but it is too small to be accounted for.
        if(class(param.new)=='try-error'){
            par1 <- try(mvrnorm(n=1,mu=param,Sigma=nearPD(cov.adapt),tol=tol1),TRUE)
        }
        param.new[index] <- par1
        if(Np>0){
            M0.sim <- param.new[grep('Mo([[:digit:]]{1})',names(param.new))]
            param.new[grep('Mo([[:digit:]]{1})',names(param.new))] <- M0.sim%%(2*pi)
            omega.sim <- param.new[grep('omega([[:digit:]]{1})',names(param.new))]
            param.new[grep('omega([[:digit:]]{1})',names(param.new))] <- omega.sim%%(2*pi)
            Omega.sim <- param.new[grep('Omega([[:digit:]]{1})',names(param.new))]
            param.new[grep('Omega([[:digit:]]{1})',names(param.new))] <- Omega.sim%%(2*pi)
        }
        logplast <- param.new[grep(paste0('per',Np),names(param.new))]
        if(!quantify & Np>0){
            logplast <- param.new[grep(paste0('per',Np),names(param.new))]
            if(any(logplast<logPmaxs & logplast>logPmins)){
                logic.per <- TRUE
            }else{
                logic.per <- FALSE
            }
            if(all(param.new>par.min & param.new<par.max) & all(logic.per)) break()
        }else{
            if(all(param.new>par.min & param.new<par.max)) break()
        }
    }
#    if(k==Ntt & FALSE){
    if(k==Ntt){
	cat('param.new=',param.new,'\n')
	cat('par.min=',par.min,'\n')
	cat('par.max=',par.max,'\n')
        cat('The times of generating the proposed parameters reach the maximum value!\n')
    }
    names(param.new) <- names(startvalue)
    return(param.new)
}

#####decide whether to accept or adjust the new period values
newper <- function(pars0,pars,pers.low,pers.up){
    ind.last <- which(names(pars)==paste0('per',Np))
    if(Np>1){
        inds = grep('per([[:digit:]]{1})',names(pars))
        ind.former <- inds[inds!=ind.last]
        per.former <- pars[ind.former]
        per0.former <- pars0[ind.former]
    }
    per.last <- pars[ind.last]
    tmp <- par.sec(per.last,pers.low,pers.up)
    logic.all <- logic1 <- tmp$logic
    pars[ind.last] <- tmp$par.new
###former period parameters
    if(length(pers.up)>1 & !fixP){
        pf.low <- pers.up[-length(pers.up)]
        pf.up <- pers.low[-1]
        if(Np>1){
            logic2 <- c()
            for(k1 in 1:length(per.former)){
                ind <- which(per0.former[k1]>pf.low & per0.former[k1]<pf.up)
                if(per.former[k1]>pf.low & per.former[k1]<pf.up){
                    logic2 <- c(logic2,TRUE)
                }else{
                    logic2 <- c(logic2,FALSE)
                }
            }
            logic.all <- all(logic1,logic2)
        }
    }
    return(list(pars=pars,logic=logic.all))
}
####
par.sec <- function(par,pars.low,pars.up){
    if(all(par>pars.low & par<pars.up)){
        log1 <- TRUE
        par.new <- par
    }else if(par< max(pars.up) & par>min(pars.low)){
        dpar.low <- min(abs(par-pars.low))
        dpar.up <- min(abs(par-pars.up))
        if(dpar.low<dpar.up){
            ind.min <- which.min(abs(par-pars.low))
            par.new <- min(pars.low[ind.min]+dpar.low,pars.up[ind.min])
        }else{
            ind.min <- which.min(abs(par-pars.up))
            par.new <- max(pars.up[ind.min]-dpar.up,pars.low[ind.min])
        }
        log1 <- TRUE
    }else{
        log1 <- FALSE
        par.new <- par
    }
    return(list(par.new=par.new,logic=log1))
}
####adaptive proposal function
###calculate covariance matrix
covariance.n0 <- function(mat,eps=tol1){
    cov.par <- Sd*cov(mat)+Sd*eps*diag(ncol(mat))
    return(cov.par)
}

covariance.rep <- function(parms,cov1,mu1,n,Nupdate=1,eps=tol1){
    Npar <- length(parms)
    if(n%%Nupdate==0){
        Sd <- 2.4^2/Npar
        mu2 <- (mu1*(n-1)+parms)/n#mu1 and mu2 are the mean vectors of parameters for n-1 and n iterations
        N <- n-1
        cov2 <- (N-1)*cov1/N+Sd*(N*mu1%*%t(mu1)-(N+1)*mu2%*%t(mu2)+parms%*%t(parms)+eps*diag(length(parms)))/N
    }else{
        cov2 <- cov1
    }
    return(cov2)
}
###run MH algorithm
run.metropolis.MCMC <- function(startvalue,cov.start,iterations,n0=10,verbose=FALSE,tem=1,bases=rep('natural',10),RVonly=FALSE){
    n0 <- max(n0,min(1000,iterations/10))
    Npar <- length(startvalue)
    chain  <-  array(dim=c(iterations+1,Npar))
    logpost = loglike = rep(NA,iterations+1)
    colnames(chain) <- names(startvalue)
    chain[1,]<- startvalue
    mu1 <- chain[1,]
    logpost.out = posterior(chain[1,],tem=tem,bases=bases,RVonly=RVonly)
    logpost[1] <- logpost.pre <- logpost.out$post
    loglike[1] <- loglike.pre <- logpost.out$loglike
    logprior.pre <- logpost.out$logprior
    dt0 <- 0
    cov.adapt <- array(data=0,dim=c(Npar,Npar))
    t.start <- proc.time()
    ind2 <- 1:length(startvalue)
    if(!is.null(par.fix)){
        ind1 <- match(par.fix,names(startvalue))
        ind2 <- ind2[-ind1]
    }
    for(i in 1:iterations){
        if(i == n0){
            cov.adapt[ind2,ind2] <- covariance.n0(mat=chain[1:n0,ind2])
        }else if(i > n0){
            cov.adapt[ind2,ind2] <- covariance.rep(chain[i,ind2],cov.adapt[ind2,ind2],mu1[ind2],i)
        }else{
            cov.adapt <- cov.start
        }
#        if(i%%(iterations/100)==0) cat('diag(cov.adapt)=',diag(cov.adapt),'\n')
        proposal  <- chain[i,]
        proposal[ind2] <- proposalfunction.simple(chain[i,ind2],cov.adapt[ind2,ind2],par.min=par.min[ind2],par.max=par.max[ind2])
#        verbose <- FALSE
#        if(i%%round(iterations/2)==0){
#            verbose <- TRUE
#        }
#        proprop <- posterior(proposal,tem=tem,bases=bases,RVonly=RVonly,verbose=verbose)
        proprop <- posterior(proposal,tem=tem,bases=bases,RVonly=RVonly,verbose=FALSE)
        logpost.prop <- proprop$post
	logprior.prop <- proprop$logprior
        loglike.prop <- proprop$loglike

        logpost.cur <- logpost.pre
        loglike.cur <- loglike.pre
	logprior.cur <- logprior.pre
        if(is.na(logpost.prop)){
            probab = 0
        }else{
            probab = exp(logpost.prop-logpost.cur)
        }
        if(runif(1)< probab){
            chain[i+1,]=proposal
###values for next proposal estimation
            logpost.pre <- logpost.prop
            loglike.pre <- loglike.prop
	    logprior.pre <- logprior.prop
        }else{
            chain[i+1,]=chain[i,]
            logpost.pre <- logpost.cur
            loglike.pre <- loglike.cur
            logprior.pre <- logprior.cur
        }
##save values
        logpost[i+1]<- logprior.pre + loglike.pre
        loglike[i+1] <- loglike.pre
        mu1 <- (mu1*i+chain[i+1,])/(i+1)
##moniter parameters
        if(verbose & i%%round(iterations/10)==0){
#            cat('i=',i,';tem=',tem,';diag(cov.adapt)=',diag(cov.adapt),';')
            cat('i=',i,';tem=',tem,';')
            cat('acceptance percentage:',100*(1-mean(duplicated(chain[1:(i+1),1:Npar]))))
            if(Np>=1){
                pars.now <- chain[i+1,]
                pp <- extract.par(pars.now,bases=bases)
                Ps <- pp$P
                Ks  <-  pp$K
                es  <-  pp$e
                cat('; period= ',Ps,'; K=', Ks)
                if(prior.type!='e0'){
                    cat('; e=',es,'\n')
                }
            }
            cat('; max(logpost)=',max(logpost[1:(i+1)]),'; maximum likelihood:',max(loglike[1:(i+1)]),'\n')
            t.start <- proc.time()
        }
        conv <- Rhat <- NULL
        if(i==iterations){
#            cat('i=',i,'\n')
###check convergence of chains
            Nsub <- 5
            mcmc <- chain[,!grepl('Mo|omega|logpost|loglike',colnames(chain))]
            chain.len <- nrow(mcmc)
            subchain.len <- floor(chain.len/Nsub)
            var.sub <- array(data=NA,dim=c(Nsub,ncol(mcmc)))#s
            mean.sub <- array(data=NA,dim=c(Nsub,ncol(mcmc)))#theta.b
            meanvar <- rep(NA,ncol(mcmc))#W
            mean.all <- rep(NA,ncol(mcmc))#theta.bb
            var.single.est <- rep(NA,ncol(mcmc))#B
            for(j1 in 1:ncol(mcmc)){
                for(i1 in 1:Nsub){
                    var.sub[i1,j1] <- var(mcmc[((i1-1)*subchain.len+1):(i1*subchain.len),j1])
                    mean.sub[i1,j1] <- mean(mcmc[((i1-1)*subchain.len+1):(i1*subchain.len),j1])
                }
                meanvar[j1] <- mean(var.sub[,j1])
                mean.all[j1] <- mean(mean.sub[,j1])
                var.single.est[j1] <- subchain.len/(Nsub-1)*sum((mean.sub[,j1]-mean.all[j1])^2)
            }
            var.est <- (subchain.len-1)/subchain.len*meanvar+1/subchain.len*var.single.est#Var.hat
            Rhat <- sqrt(var.est/meanvar)
            acceptance <- 100*(1-mean(duplicated(chain[1:i,])))
            if(any(is.na(Rhat))) Rhat[is.na(Rhat)] <- 1
            if(any(Rhat>1.1)){
                conv <- FALSE
            }else{
                conv <- TRUE
#                chain <- chain[1:i,]
#                if(i>1e5){
#                    break()
#                }
            }
        }
    }
    val <- list(out=cbind(chain,logpost,loglike),Rhat=Rhat,conv=conv,acc=acceptance)
#    cat('boundary of period=',exp(par.min[1]),exp(par.max[1]),'d\n')
#    cat('range of period=',exp(range(chain[,(Np-1)*Nkeppar+1])),'d\n')
    return(val)
}
bin.simple <- function(data,Nbin){
    x <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    xmid <- seq(min(x),max(x),length.out=Nbin+1)
    x1 <- xmid[-1]
    x2 <- xmid[-length(xmid)]
#    xbin <- (x1+x2)/2
    xbin <- sapply(1:Nbin,function(j) sum(x[x<x1[j] & x>x2[j]]/dy[x<x1[j] & x>x2[j]]^2)/sum(1/dy[x<x1[j] & x>x2[j]]^2))
    ybin <- sapply(1:Nbin,function(j) sum(y[x<x1[j] & x>x2[j]]/dy[x<x1[j] & x>x2[j]]^2)/sum(1/dy[x<x1[j] & x>x2[j]]^2))
    dybin <- sapply(1:Nbin,function(j) 1/sqrt(sum(1/dy[x<x1[j] & x>x2[j]]^2)))
    tmp <- cbind(xbin,ybin,dybin)
    tmp[which(!is.na(tmp[,1])),,drop=FALSE]
}

binning.post <- function(par.val,post.val,like.val){
    par.min<- min(par.val)
    par.max <- max(par.val)
    bin <- (par.max-par.min)/Nbins
    post.bin <- c()
    like.bin <- c()
    par.bin <- c()
    ind.na <- c()
    for(i in 1:Nbins){
        ind <- which((par.val>=par.min+(i-1)*bin) & (par.val<par.min+i*bin))
        if(length(ind)==0){
            post.bin <- c(post.bin,min(post.val))
	    like.bin <- c(like.bin,min(like.val))
            ind.na <- c(ind.na,i)
        }else{
#            post.bin <- c(post.bin,max(post.val[ind]))
#            like.bin <- c(like.bin,max(like.val[ind]))
            post.bin <- c(post.bin,NA)
            like.bin <- c(like.bin,NA)
        }
        par.bin <- c(par.bin,par.min+(2*i-1)*bin/2)
    }
    return(list(likes=like.bin,posts=post.bin,pars=par.bin,ind.na=ind.na))
}
###a more efficient way to binning parameters and log posteriors/likelihoods
binning.post2 <- function(par.val,post.val,like.val,Nbins){
    ind <- sort(par.val,index.return=TRUE)$ix
    par.sort <- par.val[ind]
    post.sort <- post.val[ind]
    like.sort <- like.val[ind]
    p1 <- hist(par.sort,breaks=seq(min(par.sort),max(par.sort),length.out=Nbins+1),plot=FALSE)
    index <- c(0,cumsum(p1$counts))
    post.max <- rep(NA,length(p1$mids))
    like.max <- rep(NA,length(p1$mids))
    ind.na <- c()
    for(j in 1:(length(index)-1)){
        if(index[j+1]>index[j]){
            post.max[j] <- max(post.sort[(index[j]+1):index[j+1]])
            like.max[j] <- max(like.sort[(index[j]+1):index[j+1]])
        }
    }
    ind.na <- which(is.na(post.max))
    return(list(likes=like.max,posts=post.max,pars=p1$mids,ind.na=ind.na))
}
kepler.tm <- function(Kp,wp,ep,Ea){
    Kp*sqrt(1-ep^2)*(cos(wp)*sqrt(1-ep^2)*cos(Ea)-sin(wp)*sin(Ea)/(1-ep*cos(Ea)))
}

kepler.ford <- function(Kp,wp,ep,Ea){
    Tp <- 2*atan(sqrt((1+ep)/(1-ep))*tan(Ea/2))
    rv <- Kp*(cos(wp+Tp)+ep*cos(wp))
    return(rv)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

tcol <- function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1,], rgb.val[2,], rgb.val[3,],
                 max = 255,
                 alpha = (100-percent)*255/100,
                 names = name)
    invisible(t.col)
}

####ref. ExoFAST: https://arxiv.org/pdf/1907.09480.pdf
calTc <- function(Tp,P,e,omega){
    theta <- 0.5*pi-omega
    E <- 2*atan(sqrt((1-e)/(1+e))*tan(theta/2))
    M <- E-e*sin(E)
    Tc <- (M%%(2*pi))/(2*pi)*P+Tp
    Tc
}

M02Tp <- function(M0,T0,P){
    T0-(M0%%(2*pi))*P/(2*pi)
}

Tp2M0 <- function(Tp,T0,P){
    (2*pi*(T0-Tp)/P)%%(2*pi)
}

Tc2M0.circular <- function(Tc,T0,P,omega){
    Tp <- Tc-(((0.5*pi-omega)/(2*pi))%%1)*P
    M0 <- (((T0-Tp)/P)%%1)*2*pi
    M0
}

data.distr <- function(x,lp=NULL,xlab,ylab,main='',oneside=FALSE,plotf=TRUE){
    xs <- seq(min(x),max(x),length.out=1e3)
    fitnorm <- fitdistr(x,"normal")
    p <- hist(x,plot=FALSE)
    xfit <- length(x)*mean(diff(p$mids))*dnorm(xs,fitnorm$estimate[1],fitnorm$estimate[2])
    ylim <- range(xfit,p$counts)
    if(plotf){
        plot(p,xlab=xlab,ylab=ylab,main=main,ylim=ylim)
        lines(xs,xfit,col='red')
    }
    xopt <- NA
    if(!is.null(lp)){
        xopt <- x[which.max(lp)]
    }
    x1=Mode(x)
    x2=mean(x)
    x3=sd(x)
    x4=skewness(x)
    x5=kurtosis(x)
    med=median(x)
    xs = sort(x)
    onesig <- 0.682
    confidence <- c(0.99,0.95,0.90,0.682)
    qs <- c((1-confidence)/2,(1+confidence)/2)
    tmp <- as.numeric(quantile(xs,qs))
    x99low = tmp[1]
    x95low = tmp[2]
    x90low = tmp[3]
    xminus.1sig =xlow = tmp[4]
    x99up = tmp[5]
    x95up = tmp[6]
    x90up = tmp[7]
    xplus.1sig = xup = tmp[8]
#    abline(v=c(x1per,x99per),col='blue')
    if(plotf){
        if(!oneside){
            legend('topleft',legend=c(as.expression(bquote('MAP ='~.(format(xopt,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3))))),bty='n')
#            legend('topright',legend=c(as.expression(bquote('MAP ='~.(format(xopt,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
            legend('topright',legend=c(as.expression(bquote('med ='~.(format(med,digit=3)))),as.expression(bquote(q[16]~'='~.(format(xlow,digit=3)))),as.expression(bquote(q[84]~'='~.(format(xup,digit=3))))),bty='n')
        }else{
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
        }
    }
    tmp <- c(xopt=xopt,x99low=x99low,x99up=x99up,x95low=x95low,x95up=x95up,x90low=x90low,x90up=x90up,xlow=xlow,xup=xup,xminus.1sig=xlow, xplus.1sig=xup, mode=x1,med=med,mean=x2,sd=x3,skewness=x4,kurtosis=x5)
    return(tmp)
}

matrix.distr <- function(y,lp=NULL,xlab,ylab,main='',oneside=FALSE,plotf=TRUE){
    tmp <- c()
    for(kk in 1:ncol(y)){
        x <- y[,kk]
        xs <- seq(min(x),max(x),length.out=1e3)
        fitnorm <- fitdistr(x,"normal")
        p <- hist(x,plot=FALSE)
        xfit <- length(x)*mean(diff(p$mids))*dnorm(xs,fitnorm$estimate[1],fitnorm$estimate[2])
        ylim <- range(xfit,p$counts)
        if(plotf){
            plot(p,xlab=xlab,ylab=ylab,main=main,ylim=ylim)
            lines(xs,xfit,col='red')
        }
        xopt <- NA
        if(!is.null(lp)){
            xopt <- x[which.max(lp)]
        }
        x1=Mode(x)
        x2=mean(x)
        x3=sd(x)
        x4=skewness(x)
        x5=kurtosis(x)
	med = median(x)
        xs = sort(x)
        x1per = max(min(xs),xs[floor(length(xs)*0.01)])
        x99per = min(xs[ceiling(length(xs)*0.99)],max(xs))
        x10per = max(min(xs),xs[floor(length(xs)*0.1)])
        x90per = min(xs[ceiling(length(xs)*0.9)],max(xs))
        xminus.1sig = max(min(xs),xs[floor(length(xs)*0.15865)])
        xplus.1sig = min(xs[ceiling(length(xs)*(1-0.15865))],max(xs))
                                        #    abline(v=c(x1per,x99per),col='blue')
        if(plotf){
            if(!oneside){
                legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3))))),bty='n')
                                        #            legend('topright',legend=c(as.expression(bquote('MAP ='~.(format(xopt,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
                legend('topright',legend=c(as.expression(bquote('MAP ='~.(format(xopt,digit=3)))),as.expression(bquote(q[10]~'='~.(format(x10per,digit=3)))),as.expression(bquote(q[90]~'='~.(format(x90per,digit=3))))),bty='n')
            }else{
                legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
            }
        }
        tmp <- cbind(tmp,c(xopt=xopt,x1per=x1per,x99per=x99per,x10per=x10per,x90per=x90per,xminus.1sig=xminus.1sig,xplus.1sig=xplus.1sig,mode=x1,med=med,mean=x2,sd=x3,skewness=x4,kurtosis=x5))
    }
    colnames(tmp) <- colnames(y)
    return(tmp)
}
####function to select output files from MCMC chain
select.file <- function(sub.out){
   fout <- NA
   if(length(sub.out)!=0){
      ind1 <- which(as.logical(sub.out[,5])==TRUE & as.numeric(sub.out[,3])>10 & as.numeric(sub.out[,3])<35)
      if(length(ind1)!=0){
        out2 <- sub.out[ind1,]
        if(length(ind1)==1){
           tmp <- out2
        }else{
#          ind2 <- which.max(as.numeric(out2[,4]))
           ind2 <- which.min(as.numeric(out2[,3]))
           tmp <- out2[ind2,]
        }
        fout <- as.character(tmp[1])
      }
  }
  if(!is.na(fout)){
    cat('Choose this mcmc results for further investigations: ', fout,'\n')
  }else{
    cat('No qualified mcmc chain for further investigations!\n')
  }
   return(fout)
}

###a function to reset -Inf element in an evidences
inf.rm <- function(E,B,flag=NA){
    ind.m <- which(!is.na(E) & E!=Inf & E!=-Inf,arr.ind=T)
    Emin <- min(E[ind.m])
    ind.b <- which(!is.na(B) & B!=Inf & B!=-Inf,arr.ind=T)
    Bmin <- min(B[ind.b])
    Edim <- dim(E)
    cat('Emin=',Emin,'\n')
    for(i in 1:Edim[3]){
        for(j in 1:Edim[2]){
            if(any(E[,j,i]==-Inf) | any(E[,j,i]==Inf) | any(is.na(E[,j,i])) | B[j,i]==Inf | B[j,i]==-Inf | is.na(B[j,i])){
                if(flag!=NA){
                    E[,j,i] <- Emin
                    B[j,i] <- 0
                }else{
                    E[,j,i] <- NA
                    B[j,i] <- NA
                }
            }
        }
    }
    return(list(E=E,B=B))
}
stemPlot <- function(x,y,pch=16,linecol=1,clinecol=1,add=FALSE,pair=FALSE,...){
    if(!add){
        plot(x,y,pch=pch,...)
    }else{
        points(x,y,pch=pch)
    }
    if(!pair){
        for (i in 1:length(x)){
            lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
        }
    }else{
        ind.pair <- which((1:as.integer(length(x)/2))%%2==1)
        ind.solid <- c(2*ind.pair-1,2*ind.pair)
        ind.dashed <- c(1:length(x))[-ind.solid]
        for (i in 1:length(x)){
            if(any(i==ind.solid)){
                lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
            }else{
                lines(c(x[i],x[i]), c(0,y[i]),col=linecol,lty=2)
            }
        }
    }
 #   lines(c(x[1]-2,x[length(x)]+2), c(0,0),col=clinecol,lty=3)
    abline(h=0,lty=2)
}
###functions
tpm <- function(priors, likes, lambda, h){
    prod.pre <- c(rep(0,h),likes[-((length(likes)-h+1):length(likes))])
    Ptpm <- sum(priors*likes/((1-lambda)*likes*priors+lambda*prod.pre))/sum(priors/((1-lambda)*likes*priors+lambda*prod.pre))
    return(log(Ptpm))
}
calc.alpha <- function(param,cov.adapt,post0,order='01'){
#        par.prop <- rbind(par.prop,par.tmp)
    post.prop <- posterior(param,tem=tem)$post
    if(order=='01'){
        if(is.na(post.prop)){
            probab = 0
        }else{
            probab = exp(post.prop-post0)
        }
    }else{
        if(is.na(post.prop)){
            probab = 1
        }else{
            probab = exp(post0-post.prop)
        }
    }
    min(1,probab)#since q1=q0, probab*q0/q1=probab; ref Chib and Jeliazkov 2001
}
####parallel computing
AMH <- function(nburn,Ns){
    if(nburn==0) nburn <- 1
    out.amh <- foreach(n = 1:Ncores, .combine='comb',.multicombine=TRUE) %dopar% {
        tmp <- run.metropolis.MCMC(startvalue,cov.start,Ns)
        list(out=cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)]),cov=tmp$cov)
#        cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)])
    }
#    return(out.amh)
    return(list(out=out.amh$out,cov=out.amh$cov))
}
####parallel chain returning covariance
AMH2 <- function(nburn,Ns){
    if(nburn==0) nburn <- 1
    out.amh <- foreach(n = 1:Ncores, .combine='comb',.multicombine=TRUE) %dopar% {
        tmp <- run.metropolis.MCMC(startvalue,cov.start,Ns)
        list(out=cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)]),cov=tmp$cov)
    }
    return(list(out=out.amh$out,cov=out.amh$cov))
}
####parallel chain with different temperture
AMH3 <- function(nburn,Ns){
    soft <- TRUE
    if(!soft){
        zeta <- tem.min
        Nlow <- floor(Ncores/2)+1
        Nup <- ceiling(Ncores/2)
        eta <- log(tem/tem.min)/log(Nlow)
        tems <- tempering(1:Nlow,zeta,eta)
        zeta <- tem
        eta <- log(1/tem)/log(Nup)
        tems <- c(tems,sort(tempering(Nup:1,zeta,eta))[-1])
    }else{
        zeta <- tem.min
        eta <- log(tem/tem.min)/log(Ncores)
        tems <- tempering(1:Ncores,zeta,eta)
    }
    tems <- sort(tems,decreasing=TRUE)
    if(nburn==0) nburn <- 1
    out.amh <- foreach(n = 1:Ncores, .combine='comb',.multicombine=TRUE) %dopar% {
        tem <- tems[n]
#        if(n==1){
#        if(n<=Ncores/2){
#        if(n<=3*Ncores/4){
#        if(n<=length(ps) & n<=Ncores/2){
        par.tmp <- pars
#        par.tmp <- par.end
        if(n<=(length(par.tmp)/Npar) & n<=(Ncores/2)){
            if(length(par.tmp)==Npar){
                startvalue <- par.tmp
            }else{
                startvalue <- par.tmp[n,]
            }
#            indP <- (Np-1)*Nkeppar+1
#            if(cov.start[indP,indP]<tol1){
#                cov.start[indP,indP] <- startvalue[indP]*1e-3
#            }

            tem <- tems0[Ntrace[n]]
#            tem <- 1
        }else{
#            tem <- tems[1]
#            tem <- tems0[1]
            tem <- 1
            source('prepare_par.R',local=TRUE)
        }
#        source('mcmc_func.R',local=TRUE)
        tmp <- run.metropolis.MCMC(startvalue,cov.start,Ns)
        list(out=cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)]),cov=tmp$cov)
    }
    return(list(out=out.amh$out,cov=out.amh$cov))
}
###tempering function; power function
tempering <- function(x,a,b){
    a*x^b
}
####MCMC stuck steps calculation
stuck <- function(arr,post,Nstuck.min){
#    index <- which(duplicated(arr[,((Np-1)*Nkeppar+1):Npar]))
    index <- which(duplicated(arr[,(Np-1)*Nkeppar+1]))
    Ns <- 1
    Nstuck <- c()
    ind <- c()
    ind.stuck <- c()
    for(k in 2:length(index)){
        if((index[k]-index[k-1])==1){
            Ns <- Ns+1
            ind <- c(ind,index[k])
        }else if(length(ind)>0){
            Nstuck <- c(Nstuck,Ns)
            ind.stuck <- c(ind.stuck,ind[length(ind)])
            Ns <- 1
            ind <- c()
        }else{
	    Ns <- 1
            ind <- c()
	}
    }
    index <- which(Nstuck>Nstuck.min)
    if(length(index)>0){
        nstuck <- Nstuck[index]
        indstuck <- ind.stuck[index]
        index1 <- which.max(post[indstuck])
        indstuck <- indstuck[index1]
        nstuck <- nstuck[index1]
    }else{
        index1 <- which.max(Nstuck)
        indstuck <- ind.stuck[index1]
        nstuck <- Nstuck[index1]
    }
    ind.max <- which.max(post)
    per.stuck <- arr[indstuck,(Np-1)*Nkeppar+1]
    per.max <- arr[ind.max,(Np-1)*Nkeppar+1]
    #shift the to the local maxima
    if(abs(per.stuck-per.max)<(0.01*per.stuck)){
        indstuck <- ind.max
    }
#else{
#        indstuck <- 1
#        nstuck <- 1
#    }
    return(list(Nstuck=nstuck,ind.stuck=indstuck))
}
# Print the hostname for each cluster member
GetClusterInfo <- function() {
  info <- Sys.info()[c("nodename", "machine")]
  paste("Node:", info[1], "with CPU type", info[2])
}

###new Plows and Pups for divided period space
period.division <- function(plows,pups,Pmax,Pmin){
    plows.new <- c()
    pups.new <- c()
    if(any(c(plows,pups)>Pmax) | any(c(plows,pups)<Pmin)){
        for(i0 in 1:length(plows)){
            if((pups[i0]<Pmin & plows[i0]<Pmin) | (plows[i0]>Pmax & pups[i0]>Pmax)){
                plows.new <- plows.new
                pups.new <- pups.new
            }else if(plows[i0]<=Pmin & pups[i0]>=Pmin & pups[i0]<=Pmax){
                plows.new <- c(plows.new,Pmin)
                pups.new <- c(pups.new,pups[i0])
            }else if(plows[i0]<=Pmin & pups[i0]>Pmax){
                plows.new <- c(plows.new,Pmin)
                pups.new <- c(pups.new,Pmax)
            }else if(plows[i0]>=Pmin & pups[i0]<=Pmax){
                plows.new <- c(plows.new,plows[i0])
                pups.new <- c(pups.new,pups[i0])
            }else if(plows[i0]>=Pmin & plows[i0]<=Pmax & pups[i0]>Pmax){
                plows.new <- c(plows.new,plows[i0])
                pups.new <- c(pups.new,Pmax)
            }
        }
    }else{
        plows.new <- plows
        pups.new <- pups
    }
    return(list(plow=plows.new,pup=pups.new))
}

###tell foreach how to combine output
comb <- function(x, ...) {
      mapply(rbind,x,...,SIMPLIFY=FALSE)
}
####calculate timie
time.calc <- function(t1){
    dur <- as.numeric(proc.time()[3]-t1[3])
    time.consumed <- paste(floor(dur/3600),'h',floor((dur%%3600)/60),'m',dur%%60,'s',sep='')
    cat('Time consumed: ',time.consumed,'\n')
}

wtb.disperse <- function(t,x,ex,dt=1,sj=0){
    t0 <- t[1]#the start time point
###note: weight w=1/(ex^2+sj^2)
    ts <- t[1]
    xs <- x[1]
    exs <- ex[1]
    tnew <- c()
    xnew <- c()
    exnew <- c()
    for(i1 in 2:length(t)){
        if((t[i1]-t0)<dt & i1<length(t)){
            ts <- c(ts,t[i1])
            xs <- c(xs,x[i1])
            exs <- c(exs,ex[i1])
        }else{
            tnew <- c(tnew,as.numeric(ts%*%(1/(exs^2+sj^2))/sum(1/(exs^2+sj^2))))
            x0 <- as.numeric(xs%*%(1/(exs^2+sj^2))/sum(1/(exs^2+sj^2)))
            xnew <- c(xnew,x0)
#            xnew <- c(xnew,mean(xs))
            ivw <- sqrt(1/sum(1/(exs^2+sj^2)))#inverse variance weighted uncertainty
            chi2r <- sum((xs-x0)^2/(exs^2+sj^2))/(length(exs)-1)#reduced chi2
            exnew <- c(exnew,sqrt(chi2r*ivw))
            ts <- t[i1]
            xs <- x[i1]
            exs <- ex[i1]
            t0 <- t[i1]
        }
    }
    return(cbind(tnew,xnew,exnew))
}

####weighted time binning
wtb <- function(t,x,ex,dt=1,sj=0){
#    t <- sort(t)
    t0 <- t[1]#the start time point
###note: weight w=1/(ex^2+sj^2)
    ts <- t[1]
    xs <- x[1]
    exs <- ex[1]
    tnew <- c()
    xnew <- c()
    exnew <- c()
    for(i1 in 2:length(t)){
        if((t[i1]-t0)<dt & i1<length(t)){
            ts <- c(ts,t[i1])
            xs <- c(xs,x[i1])
            exs <- c(exs,ex[i1])
        }else{
            tnew <- c(tnew,mean(ts))
            xnew <- c(xnew,xs%*%(1/(exs^2+sj^2))/sum(1/(exs^2+sj^2)))
#            xnew <- c(xnew,mean(xs))
            exnew <- c(exnew,sqrt(1/sum(1/(exs^2+sj^2))))
            ts <- t[i1]
            xs <- x[i1]
            exs <- ex[i1]
            t0 <- t[i1]
        }
    }
    return(cbind(tnew,xnew,exnew))
}
####weighted time binning with scattering as bin error
wtb.sigma <- function(t,x,ex,dt=1,sj=0){
#    t <- sort(t)
    t0 <- t[1]#the start time point
###note: weight w=1/(ex^2+sj^2)
    ts <- t[1]
    xs <- x[1]
    exs <- ex[1]
    tnew <- c()
    xnew <- c()
    exnew <- c()
    for(i1 in 2:length(t)){
        if((t[i1]-t0)<dt & i1<length(t)){
            ts <- c(ts,t[i1])
            xs <- c(xs,x[i1])
            exs <- c(exs,ex[i1])
        }else{
            tnew <- c(tnew,mean(ts))
            xnew <- c(xnew,xs%*%(1/(exs^2+sj^2))/sum(1/(exs^2+sj^2)))
#            xnew <- c(xnew,mean(xs))
            exnew <- c(exnew,sd(xs))
            ts <- t[i1]
            xs <- x[i1]
            exs <- ex[i1]
            t0 <- t[i1]
        }
    }
    return(cbind(tnew,xnew,exnew))
}
wtb.number <- function(t,x,ex,N=10){
###bin a fixed number of data points
    Nbin <- ceiling(length(t)/N)
    exnew <- xnew <- tnew <- rep(NA,Nbin)
    for(j in 1:Nbin){
        inds <- ((j-1)*N+1):min(j*N,length(t))
        tnew[j] <- t[inds]%*%(1/ex[inds]^2)/sum(1/(ex[inds]^2))
        xnew[j] <- x[inds]%*%(1/ex[inds]^2)/sum(1/(ex[inds]^2))
        exnew[j] <- sqrt(1/sum(1/ex[inds]^2))
    }
    return(cbind(tnew,xnew,exnew))
}

wtb.simple <- function(t,x,ex,dt=1,sj=0,N=NULL){
#    t <- sort(t)
    t0 <- t[1]#the start time point
###note: weight w=1/(ex^2+sj^2)
    tmin <- min(t)
    tmax <- max(t)
    if(is.null(N)){
        ts <- seq(tmin,tmax,by=dt)
    }else{
        ts <- seq(tmin,tmax,length.out=N+1)
    }
    out <- c()
    for(j in 1:(length(ts)-1)){
        if(j==1){
            inds <- which(t>=ts[j] & t<=ts[j+1])
        }else{
            inds <- which(t>ts[j] & t<=ts[j+1])
        }
        if(length(inds)>0){
            tmean <- mean(t[inds])
            xmean <- x[inds]%*%(1/(ex[inds]^2+sj^2))/sum(1/(ex[inds]^2+sj^2))
            exmean <- sqrt(1/sum(1/(ex[inds]^2+sj^2)))
            out <- rbind(out,c(tmean,xmean,exmean))
        }
    }
    return(out)
}

#####extract variable values from file names
extract.Nsamp <- function(f){
    Ns <- c()
    for(i in 1:length(f)){
        f1 <- gsub('.*Nsamp','',f[i])
        Ns <- c(Ns,as.integer(gsub('_.+','',f1)))
    }
    return(Ns)
}
extract.tem <- function(f){
    f1 <- gsub('.*tem','',f)
    return(as.numeric(gsub('_acc\\d+','',f1)))
}
extract.Ndata <- function(f){
    if(grepl('Ndata',f)){
        f1 <- gsub('.*Ndata','',f)
    }else{
        f1 <- gsub('.*w0_N','',f)
    }
    return(as.integer(gsub('_.+','',f1)))
}
extract.Lmax <- function(f){
    f1 <- gsub('.+negLmax','',f)
    f2 <- gsub('-.+','',f1)
    f3 <- gsub('.pdf','',f2)
    return(-as.numeric(f3))
}
extract.acc <- function(f){
    f1 <- gsub('.*acc','',f)
    return(as.numeric(f1))
}
combine.index <- function(ind1,ind2,norm=FALSE){
    if(length(ind1)!=0){
        if(!is.matrix(ind1)) ind1 <- matrix(ind1,ncol=1)
        if(!is.matrix(ind2)) ind2 <- matrix(ind2,ncol=1)
        if(nrow(ind2)<nrow(ind1)){
            ind <- cbind(ind1,c(ind2,rep(NA,nrow(ind1)-nrow(ind2))))
        }else if(nrow(ind2)>nrow(ind1)){
            ind1 <- rbind(ind1,matrix(NA,nrow=nrow(ind2)-nrow(ind1),ncol=ncol(ind1)))
            ind <- cbind(ind1,ind2)
        }else{
            ind <- cbind(ind1,ind2)
        }
    }else{
        ind <- ind2
    }
    if(norm){
        ind <- scale(ind)
    }
    if(!is.matrix(ind)){
        ind <- matrix(ind,ncol=1)
    }
    return(ind)
}
####period par transformation
period.transformation <- function(pers,period.type=period.par){
    if(period.type=='logP'){
        Pers <- exp(pers)
    }else if(period.type=='nu'){
        Pers <- 1/pers
    }else if(period.type=='P'){
        Pers <- pers
    }else{
        Pers <- NA
    }
    return(Pers)
}
period.trans2 <- function(ps,period.type=period.par){
    if(period.type=='logP'){
        Pers <- log(ps)
    }else if(period.type=='nu'){
        Pers <- 1/ps
    }else if(period.type=='P'){
        Pers <- ps
    }else{
        Pers <- NA
    }
    return(Pers)
}


msini2K <- function(msini,P,e,Ms,unit='mj'){
####P in unit of days
    G <- 4*pi^2
    if(unit=='me'){
        msini <- msini*Me2j
    }else if(unit=='ms'){
        msini <- msini/Mj2s
    }
    28.4329*msini*(msini*Mj2s+Ms)^(-2/3)*(P/365.25)^(-1/3)/sqrt(1-e^2)
}

m2K <- function(m,P,e,Ms,Inc,unit='mj'){
####P in unit of days
##https://www.researchgate.net/publication/253789798_Radial_Velocity_Techniques_for_Exoplanets
    G <- 4*pi^2
    m2 <- m
    if(unit=='me'){
        m2 <- m*Me2j
    }else if(unit=='ms'){
        m2 <- m/Mj2s
    }
    msini <- m2*sin(Inc)
    28.4329*msini*(m2*Mj2s+Ms)^(-2/3)*(P/365.25)^(-1/3)/sqrt(1-e^2)
}

calcmass <- function(K,P,e,Mstar,Inc){
    K <- K/1e3/4.74047#from m/s to au/yr
    P <- P/365.25
    a1 <- (K/sin(Inc))^2/(4*pi^2)*(1-e^2)*P^(2/3)
    mp0 <- mp <- sqrt(a1*Mstar^(4/3))
    for(j in 1:100){
        mp <- sqrt(a1*(Mstar+mp)^(4/3))
        if(all(abs(mp-mp0)/mp0<0.0001)) break()
        mp0 <- mp
    }
    return(mp)
}
k2m <- function(K,P,e,Ms,Inc=NULL,Niter=100,tol=1e-6,more=FALSE,type=2){
###If Inc is given, k2m will determine absolute mass
###If Inc is not given, k2m will approximately give msini if m is small
    Me2s <- 3.003e-6#Earth mass in solar unit
    Mj2s <- 1/1048
    if(is.null(Inc)){
        sinI <- 1
    }else{
        sinI <- sin(Inc)
    }
    K <- K/1e3/4.74047#from m/s to au/yr
    P <- P/365.25#yr
    ##https://ui.adsabs.harvard.edu/abs/2000A%26AS..145..161P/abstract
    ##http://exoplanets.astro.yale.edu/workshop/EPRV/Bibliography_files/Radial_Velocity.pdf
    type <- 2
    if(type==1){
        a1 <- (K/sinI)^2/(4*pi^2)*(1-e^2)*P^(2/3)#; previous version, problematic
        mp0 <- mp <- sqrt(a1*Ms^(4/3))#; previous version, problematic
    }else if(type==2){
        a1sini <- P*K*sqrt(1-e^2)/(2*pi)
        a1 <- a1sini/sinI#semi-major axis of the primary star around mass center
        ##assume mp is small and get initial value
        mp0 <- mp <- a1*P^(-2/3)*Ms^(2/3)
    }else if(type==3){
        a1sini <- P*K*sqrt(1-e^2)/(2*pi)
        a1 <- a1sini/sinI#semi-major axis of the primary star around mass center
        mp0 <- mp <- (a1sini^3/P^2*Ms^2)^(1/3)/sinI
    }
    for(j in 1:Niter){
        if(type==1)  mp <- sqrt(a1*(Ms+mp)^(4/3))
        if(type==2)  mp <- a1*P^(-2/3)*(Ms+mp)^(2/3)
        if(type==3)  mp <- (a1sini^3/P^2*(Ms+mp)^2)^(1/3)/sinI
        if(all(abs(mp-mp0)/mp0<tol)) break()
        mp0 <- mp
    }
    Mpj <- mp/Mj2s
    Mpe <- Mpj*Mj2s/Me2s
    if(more){
        a <- (P^2*(Ms+mp))^(1/3)#au                                                                                          as <- a*mp/(Ms+mp)
        ap <- a*Ms/(Ms+mp)
        return(list(mj=Mpj,me=Mpe,ms=mp,ap=ap,a=a,as=as))
    }else{
        return(list(mj=Mpj,me=Mpe,ms=mp))
    }
}

fix <- function(x,par){
    if(length(par)>0){
        for(j in 1:length(par)){
            ind <- grep(par[j],names(startvalue))
            if(is.matrix(x)){
                x[ind,] <- 0
                x[,ind] <- 0
            }else{
#                if(grepl('per',par[j])){
                    x[ind] <- startvalue[ind]
#                }else{
#                    x[ind] <- 0
#                }
            }
        }
    }
    return(x)
}
lag.sta <- function(Ms,m1,a1,e1,m2,e2){
    mu1 <- m1/Ms
    mu2 <- m2/Ms
    alpha <- mu1+mu2
    gamma1 <- (1-e1^2)^0.5
    gamma2 <- (1-e2^2)^0.5
    f <- function(delta){
        alpha^-3*(mu1+mu2/delta^2)*(mu1*gamma1+mu2*gamma2*delta)^2-1-3^(4/3)*mu1*mu2/alpha^(4/3)
    }
    d <- uniroot.all(f,c(0,10))
    return(d^2*a1)
}
show.peaks <- function(ps,powers,levels){
    ind <- which(powers==max(powers) | (powers>(max(powers)-log(100)) & powers>levels[3]))
    pmax <- ps[ind]
    ppmax <- powers[ind]
    j0 <- 1
    p0 <- pmax[1]
    pp0 <- ppmax[1]
    pms <- p0
    pos <- pp0
    if(length(pmax)>1){
        for(j in 2:length(pmax)){
            if(abs(pmax[j]-p0) < 0.1*p0){
                if(ppmax[j]>pp0){
                    j0 <- j
                    p0 <- pmax[j]
                    pp0 <- ppmax[j0]
                    pms[length(pms)] <- p0
                    pos[length(pos)] <- pp0
                }
		    }else{
                j0 <- j
                p0 <- pmax[j]
                pp0 <- ppmax[j0]
                pms <- c(pms,p0)
                pos <- c(pos,pp0)
            }
        }
    }else{
        pms <- pmax
        pos <- ppmax
    }
    return(cbind(pms,pos))
}
check.window <- function(t,Dt,Nbin){
    n <- Nbin-1
    dt <- (max(t)-min(t)-Dt)/n
    tstart <- min(t)+(0:n)*dt
    tend <- min(t)+(0:n)*dt+Dt
    Ns <- c()
    for(j in 1:n){
        ind <- which(t>tstart[j] & t<tend[j])
        Ns <- c(Ns,length(ind))
    }
    return(Ns)
}

combine.list <- function(sets){
    N <- length(sets)
    data <- c()
    for(j in 1:N){
        data <- c(data,sets[[j]])
    }
    return(data)
}

###optimize angular parameters
optimize.ang <- function(mcmc){
    omega = omega1 = omega2 = omega3 = mcmc.out[,grep('omega([[:digit:]]{1})',colnames(mcmc.out))]
    if(any(grepl('Mo',colnames(mcmc.out)))){
        Mo = Mo1 = Mo2 = Mo3 = mcmc.out[,grep('Mo([[:digit:]]{1})',colnames(mcmc.out))]
    }
    if(Np==1){
        omega = omega1 = omega2 = omega3 = matrix(omega1)
        if(any(grepl('Mo',colnames(mcmc.out)))){
            Mo = Mo1 = Mo2 = Mo3 = matrix(Mo1)
        }
    }
    for(j in 1:Np){
        ind.o2 <- which(omega1[,j]>pi)
        ind.o3 <- which(omega1[,j]<pi)
        omega2[ind.o2,j] <- omega2[ind.o2,j]-2*pi
        omega3[ind.o3,j] <- omega3[ind.o3,j]+2*pi
        k <- which.min(c(sd(omega1[,j]),sd(omega2[,j]),sd(omega3[,j])))
        if(k==1){
            omega[,j] <- omega1[,j]
        }else if(k==2){
            omega[,j] <- omega2[,j]
        }else{
            omega[,j] <- omega3[,j]
        }
        if(any(grepl('Mo',colnames(mcmc)))){
            ind.o2 <- which(Mo1[,j]>pi)
            ind.o3 <- which(Mo1[,j]<pi)
            Mo2[ind.o2,j] <- Mo2[ind.o2,j]-2*pi
            Mo3[ind.o3,j] <- Mo3[ind.o3,j]+2*pi
            k <- which.min(c(sd(Mo1[,j]),sd(Mo2[,j]),sd(Mo3[,j])))
            if(k==1){
                Mo[,j] <- Mo1[,j]
            }else if(k==2){
                Mo[,j] <- Mo2[,j]
            }else{
                Mo[,j] <- Mo3[,j]
            }
            mcmc[,grep('Mo([[:digit:]]{1})',colnames(mcmc))] <- Mo
        }
    }
    mcmc[,grep('omega([[:digit:]]{1})',colnames(mcmc.out))] <- omega
    mcmc
}
error.ellipse <- function(x,y,covar,percent=95){
#http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
#http://www.r-tutor.com/elementary-statistics/probability-distributions/chi-squared-distribution
    tmp <- eigen(covar)
    lambda <- abs(tmp$values)
    v1 <- tmp$vectors[,1]
    angle <- atan(v1[2]/v1[1])
    if(angle < 0) angle <- angle+2*pi
    s <- qchisq(percent/100, df=2)
    a <- sqrt(s*lambda[1])
    b <- sqrt(s*lambda[2])
    t <- seq(0, 2*pi, by=pi/1000)
    xc <- x
    yc <- y
    xt <- a*cos(t)
    yt <- b*sin(t)
    rot <- array(c(cos(angle),sin(angle),-sin(angle),cos(angle)),dim=c(2,2))
    xy <- rot%*%rbind(xt,yt)
    xt <- xy[1,]+xc
    yt <- xy[2,]+yc
    return(cbind(xt,yt))
}
ellipse.edge <- function(covar,xs,ys,ic){
                                        # covariance of x,y
                                        # xs and ys are the coordinates of the center of the ellipse and two nearby points on the reference orbit to determine the slope
                                        # ... are any arguments that can be passed to function lines
    tmp <- eigen(covar)
    lambda <- abs(tmp$values)
    v1 <- tmp$vectors[,1]
    angle <- atan(v1[2]/v1[1])
    if(angle < 0) angle <- angle+2*pi
    a <- sqrt(5.991*lambda[1])
    b <- sqrt(5.991*lambda[2])
    t <- seq(0, 2*pi, by=pi/1000)
    xc <- xs[ic]
    yc <- ys[ic]
    xt <- a*cos(t)
    yt <- b*sin(t)
    rot <- array(c(cos(angle),sin(angle),-sin(angle),cos(angle)),dim=c(2,2))
    xy <- rot%*%rbind(xt,yt)
    xt <- xy[1,]+xc
    yt <- xy[2,]+yc
    s1 <- mean(diff(ys))/mean(diff(xs))
    ind.up <- which(yt>yc)
    ind.low <- which(yt<=yc)
    s2up <- diff(yt[ind.up])/diff(xt[ind.up])
    s2low <- diff(yt[ind.low])/diff(xt[ind.low])
    ind1 <- ind.up[which.min(abs(s2up-s1))]
    ind2 <- ind.low[which.min(abs(s2low-s1))]
    xup <- xt[ind1]
    xlow <- xt[ind2]
    yup <- yt[ind1]
    ylow <- yt[ind2]
###determine semi-ellipse divided by the orthogonal line
    slope <- (yup-ylow)/(xup-xlow)
    intercept <- (ylow*xup-xlow*yup)/(xup-xlow)
    ind.semi.up <- which(yt>=(intercept+slope*xt))
    ind.semi.low <- which(yt<=(intercept+slope*xt))
    xy <- cbind(xt[c(ind1,ind2)],yt[c(ind1,ind2)])
###sort the index
    t2 <- (t+angle)%%(2*pi)
    ind.semi.up <- ind.semi.up[sort(t2[ind.semi.up],index.return=TRUE)$ix]
    ind.semi.low <- ind.semi.low[sort(t2[ind.semi.low],index.return=TRUE)$ix]
    xs.up <- xt[ind.semi.up]
    xs.low <- xt[ind.semi.low]
    ys.up <- yt[ind.semi.up]
    ys.low <- yt[ind.semi.low]
##output
    return(list(edge=xy,xt=xt,yt=yt,xup=xs.up,yup=ys.up,xlow=xs.low,ylow=ys.low))
#    lines(xt,yt,...)
#    points(xt[c(ind1,ind2)],yt[c(ind1,ind2)],...)
}
rm.out <- function(x,y,alpha=3){
    r <- sqrt(x^2+y^2)
    ind.out <- which(abs(diff(r))>abs(median(diff(r)))+alpha*sd(diff(r)))
    ind.in <- (1:length(r))[-ind.out]
    for(k in ind.out){
        ind.rep <- ind.in[which.min(abs(ind.in-k))]
        x[k] <- x[ind.rep]
        y[k] <- y[ind.rep]
    }
    return(cbind(x,y))
}
