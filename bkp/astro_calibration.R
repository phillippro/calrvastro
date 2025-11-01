###
#• For straightforward calibration of GDR2 data, employ SS-based values for sources with G <10.5mag and NST-based values for sources with G >10.5 mag.
#• For precise calibration of GDR2 data, utilize the Python scripts provided in the appendix to address bias on a case-by-case basis.
#• Calibrate the 2-parameter GDR1 solutions using parameters determined with NSTP2.
#• For GDR1 targets with 5-parameter solutions, use the SSP5 values to calibrate G < 10.5 mag stars and use the NSTP2 values to calibrate G > 10.5 mag stars.
#• Use TYC positions without correction, in conjunction with GDR2 and GDR3, to constrain orbits.
#• Calibrate HIP1 and HIP2 data using either non-zero ⃗ε or non-zero ω⃗ given by Lindegren et al. (2016) but not both, the ICRF2 could be transformed to GDR3 frame using the calibration parameters given by Charlot et al. (2020); Lunz et al. (2023).
#• Quadratically add an astrometric jitter of 2.16 mas (or mas/yr) to HIP2 astrometry.
if(!is.null(out$astrometry) & calibrate){
    cat('calibration!\n')
    if(any(colnames(out$astrometry)=='catalog')){
    ih <- which(out$astrometry[,'catalog']=='HIP')
    it <- which(out$astrometry[,'catalog']=='TYC')
    i1 <- which(out$astrometry[,'catalog']=='GDR1')
    i2 <- which(out$astrometry[,'catalog']=='GDR2')
    i3 <- which(out$astrometry[,'catalog']=='GDR3')
    }
    ih <- out$ihip
    it <- out$ihip
    i1 <- out$igdr1
    i2 <- out$igdr2
    i3 <- out$igdr3
###the values of calibration parameters are from Feng+2024
    if(length(ih)>0){
        par.HIPtoVLBI2015 <- c(0,0,0,0.126,-0.185,-0.076,0.089)
        par.VLBI2015toVLBI2020 <- c(0.008,0.015,0,0,0,0,0)
        par.VLBI2020toGDR3 <- c(0.226,0.327,0.168,0.022,0.065,-0.016,0)
        out$astrometry[ih,] <- Calibrate2GDR3(out$astrometry[ih,],par.HIPtoVLBI2015,dt=-23.75,gamma=-1)
        out$astrometry[ih,] <- Calibrate2GDR3(out$astrometry[ih,],par.VLBI2015toVLBI2020,dt=-5.015,gamma=-1)
        out$astrometry[ih,] <- Calibrate2GDR3(out$astrometry[ih,],par.VLBI2020toGDR3,dt=-4.015,gamma=-1)
    }
    if(length(i1)>0){
        if(is.na(out$astrometry[i1,'parallax'])){
            par.cal <- c(0.00,  -0.13,  -0.01,  0,0,0,0)
        }else{
            if(gmag>10.5) par.cal <- c(0.15,  -0.45,  -0.05,  -0.05,-0.35,-0.05,0.01)
            if(gmag<=10.5) par.cal <- c(0.39,-0.17,0.12,0.02,-0.03,0.02,0.00)
        }
        out$astrometry[i1,] <- Calibrate2GDR3(out$astrometry[i1,],par.cal,dt=2015-2016,gamma=-1)
    }
    if(length(i2)>0){
        par.cal <- c(-0.08,  -0.02,  -0.01,  0.00,  -0.07, 0.01, 0.02)
        out$astrometry[i2,] <- Calibrate2GDR3(out$astrometry[i2,],par.cal,dt=2015.5-2016,gamma=-1)
    }
}
