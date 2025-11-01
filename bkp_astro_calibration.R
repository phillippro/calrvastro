##DR2 zero-point parallax and frame rotation from Lindegren et al. 2020a; https://www.aanda.org/articles/aa/pdf/2020/01/aa36161-19.pdf
##calculate the zero-point parallax as a function of color and magnitude (https://pypi.org/project/gaiadr3-zeropoint/)
##GDR3 zero-point parallax
##-0.029mas from Groenewegen 2021 (ui.adsabs.harvard.edu/abs/2021A&A...654A..20G/abstract)
##-0.068mas from Fabricius+2021 (https://www.aanda.org/articles/aa/pdf/2021/05/aa39834-20.pdf)
##Kervella et al. 2019; https://doi.org/10.1051/0004-6361/201834371 and Brandt et al. 2018 have similar frame rotation values but are derived from third proper motion and thus not used
##GDR2 rotation: wx=-0.077,wy=-0.096,wz=-0.002
##EDR3/DR3 frame rotation; Fabricius et al. 2021
if(FALSE){
rot2 <- c(wx=-0.077,wy=-0.096,wz=-0.002)
rot3 <- c(wx=-0.120,wy=0.173,wz=0.09)
##transform from Hipparcos frame to Gaia DR3 frame
#out$astrometry[out$ihip,] <- frame.rotation(out$astrometry[out$ihip,],rot=-rot3,dt=2016-1991.25)
##transform GDR3 from GDR3 frame to Hipparcos frame
out$astrometry[out$igdr3,] <- frame.rotation(out$astrometry[out$igdr3,],rot=rot3,dt=2016-1991.25)
##transform from Gaia DR2 frame to Gaia DR3 frame and then to Hipparcos frame
out$astrometry[out$igdr2,] <- frame.rotation(out$astrometry[out$igdr2,],rot=rot2-rot3,dt=2016-2015.5)
out$astrometry[out$igdr2,] <- frame.rotation(out$astrometry[out$igdr2,],rot=rot3,dt=2016-1991.25)
##correct for parallax zero point
out$astrometry[out$igdr2,'parallax'] <- out$astrometry[out$igdr2,'parallax']+0.05
out$astrometry[out$igdr3,'parallax'] <- out$astrometry[out$igdr3,'parallax']+0.068
}else{
    par.cal <- c(0.081167303,  0.026200736,  0.014835092,  0.001022954,  0.064735816, -0.012307563, -0.019974004)
Calibrate2GDR3()
    out$astrometry[out$igdr2,] <- GaiaCalibrate(out$astrometry[out$igdr2,],par.cal,dt=2015.5-2016)
}
