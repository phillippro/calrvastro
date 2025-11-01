args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
   star <- args[1]
   pmra <- as.numeric(args[2])
   pmdec <- as.numeric(args[3])
   plx <- as.numeric(args[4])
   vr <- as.numeric(args[5])
   dT <- as.numeric(args[6])
}else{
   star <- 'HD209100'
   pmra <- 3960.930
   pmdec <- -2539.230
   plx <- 276.0600
   vr <- -40.50448
   dT <- 3
}
mu <- sqrt(pmra^2+pmdec^2)
#Astrometric determination of WD radial velocities with Gaia?
#https://www.researchgate.net/publication/232244606_Astrometric_determination_of_white_dwarf_radial_velocities_with_Gaia
#https://www.aanda.org/articles/aa/pdf/2021/08/aa41344-21.pdf
#https://articles.adsabs.harvard.edu/pdf/1999A%26A...348.1040D
#Table1.2.2: https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf/99adf6e3-6893-4824-8fc2-8d3c9cbba2b5
#**important paper: https://arxiv.org/pdf/1208.3048.pdf

A <- 9.7779222168e8
mudot <- function(mu,plx,vr){
##mu: mas; plx: mas; vr: km/s
    -2*mu*plx*vr/A#mas/yr^2
}
plxdot <- function(mu,plx,vr){
##
    -plx^2*vr/A#mas/yr
}
rvdot <- function(mu,plx,vr){
##https://ui.adsabs.harvard.edu/abs/1977VA.....21..289V/abstract
    2.3e-8*mu^2/plx
}
mdot <- function(plx,vr){
##https://ui.adsabs.harvard.edu/abs/1977VA.....21..289V/abstract
    2.17e-9*vr*plx
}
##https://ui.adsabs.harvard.edu/abs/1977VA.....21..289V/abstract

cat(star,'PA effect in proper motion over',dT,'years:',mudot(mu,plx,vr)*dT,'mas/yr\n')
cat(star,'PA effect in position over',dT,'years:',0.5*mudot(mu,plx,vr)*dT^2,'mas\n')
cat(star,'PA effect in RV over ',dT,'years:',rvdot(mu,plx,vr)*dT*1e3,'m/s\n')
cat(star,'PA effect in parallax over',dT,'years:',plxdot(mu,plx,vr)*dT,'mas\n')
cat(star,'PA effect in apparent magnitude over',dT,'years:',log(10)*2.5*mdot(plx,vr)*dT*1e6,'ppm\n')

