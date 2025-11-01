mass2lum <- function(m){
    if(m<0.43){
        l <- 0.23*m^2.3
    }else if(m<2){
        l <- m^4
    }else if(m<55){
        l <- 1.4*m^3.5
    }else{
        l <- 32000*m
    }
    l
}
m2l.eker15 <- function(m){
#https://www.aanda.org/articles/aa/pdf/2018/11/aa34109-18.pdf
#https://iopscience.iop.org/article/10.1088/0004-6256/149/4/131/pdf
    if(m>0.38 & m<1.05){
        ll <- 4.841*log10(m)-0.026
    }else if(m>1.05 & m<2.4){
        ll <- 4.328*log10(m)-0.002
    }else if(m>2.4 & m<7){
        ll <- 3.962*log10(m)+0.120
    }else if(m>7 & m<32){
        ll <- 2.726*log10(m)+1.237
    }else if(m<0.38){
        ll <- 2.440*log10(m)-1.03494
    }
    10^ll
}
#https://www.aanda.org/articles/aa/pdf/2018/11/aa34109-18.pdf

mlr <- function(L,dL){
#http://iopscience.iop.org/article/10.1088/0004-6256/149/4/131/pdf
    logL <- log10(L)
    logL12 <- 4.841*log10(c(0.38,1.05))-0.026
    logL23 <- 4.328*log10(c(1.05,2.4))-0.002
    logL34 <- 3.962*log10(c(2.4,7))+0.120
    logL45 <- 2.726*log10(c(7,32))+1.237
    sigmas <- c(0.121,0.108,0.165,0.158)
    alphas <- c(4.841,4.328,3.962,2.726)
    Mlim <- c(0.38,1.05,2.4,7,7,32)
#cat('L1=',L1,';L2=',L2,';L3=',L3,';L4=',L4,';L=',L,'\n')
    dMs <- Ms <- c()
    for(j in 1:length(L)){
        if(logL[j]>=logL12[1] & logL[j]<logL12[2]){
#            cat('Low mass\n')
            M <- 10^((log10(L[j])+0.026)/4.841)
            s <- sigmas[1]
            alpha <- alphas[1]
        }else if(logL[j]>=logL12[2] & logL[j]<logL23[2]){
#            cat('Intermediate mass\n')
            M <- 10^((log10(L[j])+0.002)/4.328)
            s <- sigmas[2]
            alpha <- alphas[2]
        }else if(logL[j]>=logL23[2] & logL[j]<logL34[2]){
#            cat('High mass\n')
            M <- 10^((log10(L[j])-0.120)/3.962)
            s <- sigmas[3]
            alpha <- alphas[3]
        }else if(logL[j]>=logL34[2] & logL[j]<logL45[2]){
#            cat('Very high mass\n')
            M <- 10^((log10(L[j])-1.237)/2.726)
            s <- sigmas[4]
            alpha <- alphas[4]
        }
        dL <- sqrt(dL^2+((s/0.4343)*L)^2)#add fitting residual
        dM <- dL/L*M/alpha
        Ms <- c(Ms,M)
        dMs <- c(dMs,dM)
    }
    return(cbind(Ms,dMs))
}
mag2mass <- function(M,dM){
#absolute magnitude to mass in solar mass
#http://iopscience.iop.org/article/10.3847/0004-6256/152/5/141/meta
    if(M>9){
        cs <- c(0.19226,-0.050737,0.010137, -0.00075399,-1.9858e-5)
        dcs <- c(0.000424, 0.000582, 0.00021, 4.57e-5,1.09e-5)
        x0 <- 13
        u <- M-x0
        du <- dM
        mass <- cs[1]+cs[2]*u+cs[3]*u^2+cs[4]*u^3+cs[5]*u^4
        dmass <- dcs[1]+cs[2]*du+dcs[2]*u+dcs[3]*u^2+2*u*cs[3]*du+3*u^2*cs[4]*du+u^3*dcs[4]+dcs[5]*u^4+4*u^3*du*cs[5]
        dmass <- sqrt(0.023^2+dmass^2)
    }else{
                                        #Mass–luminosity relation of intermediate-mass stars by O. Yu. Malkov 2007
        mass <- exp(0.525 - 0.147*M + 0.00737*M^2)
        dmass1 <- -0.147*dM+0.00737*2*M*dM
        dmass <- sqrt(0.05^2+dmass1^2)
    }
    return(cbind(mass,dmass))
}
jy2jd.decimal <- function(year){
###2010.0 = jD 2455197.5 from Gaia
###https://www.cosmos.esa.int/documents/29201/1645651/GDR2_DataModel_draft.pdf
    (year-2010.0)*365.25+2455197.5
}
blackbody <- function(lambda,pars){
    c <- 299792458#m/s
    nu <- c/(lambda*1e-6)
    logA <- pars['logA']
    T <- pars['T']
    beta <- pars['beta']
    tau0 <- pars['tau0']
    h <- 6.62607004e-34
    k <- 1.38064852e-23
    nu0 <- 1e6*c/850
    tau <- tau0*(nu/nu0)^beta
    S <- 1e26*exp(logA)*(1-exp(-tau))*2*h*nu^3/c^2/(exp(h*nu/(k*T))-1)
    if(logFit) S <- log(S)
    return(S)
}
