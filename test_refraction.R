###calculate dRV from delay
###assume straight line delay zero
###only defection matters
###assume the height of troposphere
Htropo <- 5e4#50km
###assume 50 arcmin change of refraction during one night
dR <- 10*60/206265/12/3600#rad/s
###refraction Doppler shift
zref <- Htropo*dR/3e8#s/s
##refraction RV
rv.ref <- zref*3e8
cat('RV due to change of refraction=',zref,';dRV=',rv.ref,'m/s\n')

###derived from dt/dtau=z
###A: upper height limit of troposhere
###zref <- vOA*R
vGO <- 2*pi*6000e3/24/3600#m/s
vGA <- 2*pi*(6000+50)*1e3/24/3600#m/s
vOA <- vGA-vGO#m/s
###R=1 arcmin
R <- 60/206265
rv.ref1 <- vOA*R
cat('constant refraction RV=',rv.ref1,'m/s\n')

###conclusion the
