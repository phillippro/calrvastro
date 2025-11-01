library(pracma)
###read ccd error derived from al error; ref. Holl et al. 2023
sigma.ccd.gdr3 <- read.csv('data/gdr3_ccd_al.csv')
gmags <- c(0,sigma.ccd.gdr3[,1])
sc3s <- c(sigma.ccd.gdr3[1,2],sigma.ccd.gdr3[,2])#extrapolate to brighter stars
sc3 <- approxfun(gmags,sc3s)#ccd error as a function of gmag

sigma.al <- read.csv('data/gdr2_al.csv')
gmags <- c(0,sigma.al[,1])
sa2s <- c(sigma.al[1,2],sigma.al[,2])#extrapolate to brighter stars
sa2 <- approxfun(gmags,sa2s)#ccd error as a function of gmag

tab.u0.gdr2.2d  <- read.csv('data/table_u0_g_col.txt')
g22s <- unique(tab.u0.gdr2.2d[,1])
c22s <- unique(tab.u0.gdr2.2d[,2])
u22s <- array(NA,dim=c(length(c22s),length(g22s)))
for(i in 1:length(g22s)){
    ind  <- which(tab.u0.gdr2.2d[,1]==g22s[i])
    u22s[,i] <- tab.u0.gdr2.2d[ind,3]
}
u0.gc22 <- function(xp,yp) interp2(g22s,c22s, u22s, xp, yp, method = "linear")

Nbin <- 8#8 ccd transits per FOV transit on average

ind2 <- out$cat.ind[[which(out$cats=='GDR2')]]
ind3 <- out$cat.ind[[which(out$cats=='GDR3')]]
if(!is.na(astrometry[igdr2,'astrometric_matched_observations'])){
    Nfov2 <- astrometry[igdr2,'astrometric_matched_observations']
    Nccd2 <- Nfov2*Nbin
}
if(!is.na(astrometry[igdr3,'astrometric_matched_transits'])){
    Nfov3 <- astrometry[igdr3,'astrometric_matched_transits']
    Nccd3 <- Nfov3*Nbin
}
br <- astrometry[igdr2,'bp_rp']
u02 <- u0.gc22(gmag,br)
sal2 <- sa2(gmag)
sccd3 <- sc3(gmag)#ccd error
sfov3 <- sccd3/sqrt(Nbin)#fov erro
sccd2 <- sqrt((sal2*u02)^2+0.18^2)
sfov2 <- sccd2/sqrt(Nbin)#fov error
ruwe3 <- astrometry[igdr3,'ruwe']
ruwe2 <- astrometry[igdr2,'ruwe']
