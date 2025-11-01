u0.table <- read.csv('table_u0_g_col.txt')
###number of transits
calc.chi2 <- function(dabs,eabs){
    sum(dabs^2/eabs^2)
}

calc.gof <- function(chi2,nu){
    sqrt(9*nu/2)*((chi2/nu)^(1/3)+2/(9*nu)-1)
}

calc.jitter <- function(dabs,eabs,nu){
    p <- sum(dabs^2)/nu-eabs^2
    if(p>0){
        sqrt(p)
    }else{
        0
    }
}

calc.ruwe <- function(chi2,nu){
    sqrt(chi2/nu)
}

###fig. 9 of Lindegren+2018
gs <- 3:21
ss2 <- c(2  ,  2,1.8,0.5,0.3,0.3, 0.3, 0.3, 0.3,0.25,0.25, 0.2, 0.3,0.4,0.7,1.2,2.4,4.5,9)
ss3 <- c(0.5,0.5,0.4,0.2,0.2,0.2,0.16,0.16,0.18,0.16,0.17,0.19,0.25,0.4,0.7,1.2,2.4,4.5,9)
g2s.dr2 <- approxfun(gs,ss2)
g2s.dr3 <- approxfun(gs,ss3)
sigma.mag <- function(g){
    dloge <- log10(7)-log10(0.05)
    dmag <- 21-12
    alpha <- dloge/dmag
    sigma <- rep(0.05,length(g))
    ind <- which(g>12)
    if(length(ind)>0) sigma[ind] <- 0.05+10^(alpha*(g[ind]-12))
    sigma
}
duplicate <- TRUE
#duplicate <- FALSE
fitting <- TRUE
#fitting <- FALSE
additive <- TRUE
#additive <- FALSE
ff <- '/Users/ffeng/Documents/projects/dwarfs/data/hg/'
dr2 <- read.csv(paste0(ff,target,'_gaia2_hip.csv'))
dr3 <- read.csv(paste0(ff,target,'_gaia3_hip.csv'))
#ns <- c('astrometric_matched_observations','astrometric_gof_al','astrometric_chi2_al','astrometric_excess_noise','astrometric_excess_noise_sig','ruwe')#'astrometric_weight_al',
ns <- c('astrometric_matched_observations','astrometric_excess_noise','ruwe')
#sigma.fov <- g2s(dr3['phot_g_mean_mag'])#mas
sigma.fov <- sigma.mag(dr3['phot_g_mean_mag'])
sigma2.res <- g2s.dr2(dr2['phot_g_mean_mag'])
sigma3.res <- g2s.dr3(dr3['phot_g_mean_mag'])
dt2 <- (out$gost[out$cat.ind[[1]],1]-gdr3.epoch[1])/365.25
dt3 <- (out$gost[out$cat.ind[[2]],1]-gdr3.epoch[1])/365.25
if(fitting){
    dra2 <- lm(reflex.gost[out$cat.ind[[1]],'dra']~dt2)$residuals
    ddec2 <- lm(reflex.gost[out$cat.ind[[1]],'ddec']~dt2)$residuals
    dra3 <- lm(reflex.gost[out$cat.ind[[2]],'dra']~dt3)$residuals
    ddec3 <- lm(reflex.gost[out$cat.ind[[2]],'ddec']~dt3)$residuals
}else{
    dra2 <- reflex.gost[out$cat.ind[[1]],'dra']
    dra3 <- reflex.gost[out$cat.ind[[2]],'dra']
    ddec2 <- reflex.gost[out$cat.ind[[1]],'ddec']
    ddec3 <- reflex.gost[out$cat.ind[[2]],'ddec']
}
ruwe2 <- ruwe3 <- chi2.dr2 <- chi2.dr3 <- jitter2 <- jitter3 <- c()
for(j in 1:1000){
    d2 <- dra2*sin(out$gost[out$cat.ind[[1]],'psi'])+ddec2*cos(out$gost[out$cat.ind[[1]],'psi'])
    d3 <- dra3*sin(out$gost[out$cat.ind[[2]],'psi'])+ddec3*cos(out$gost[out$cat.ind[[2]],'psi'])
    if(additive){
        d2 <- d2+rnorm(length(dra2),0,sigma2.res)
        d3 <- d3+rnorm(length(dra3),0,sigma3.res)
    }
    if(duplicate){
        dabs2 <- dabs3 <- c()
        ts2 <- ts3 <- c()
        for(d in d2){
            dabs2 <- c(dabs2,rnorm(8,d,sigma.fov))
        }
        for(d in d3){
            dabs3 <- c(dabs3,rnorm(8,d,sigma.fov))
        }
        for(t in dt2){
            ts2 <- c(ts2,rep(t,8))
        }
        for(t in dt3){
            ts3 <- c(ts3,rep(t,8))
        }
        nu2 <- 8*length(out$cat.ind[[1]])-5
        nu3 <- 8*length(out$cat.ind[[2]])-5
    }else{
        ts2 <- dt2
        ts3 <- dt3
        dabs2 <- d2
        dabs3 <- d3
        nu2 <- length(out$cat.ind[[1]])-5
        nu3 <- length(out$cat.ind[[2]])-5
    }
    chi2.dr20 <- calc.chi2(dabs2,sigma.res)
    chi2.dr30 <- calc.chi2(dabs3,sigma.res)
    chi2.dr2 <- c(chi2.dr2,chi2.dr20)
    chi2.dr3 <- c(chi2.dr3,chi2.dr30)
    jitter2 <- c(jitter2,calc.jitter(dabs2,sigma.fov,nu2))
    jitter3 <- c(jitter3,calc.jitter(dabs3,sigma.fov,nu3))
    ruwe2 <- c(ruwe2,calc.ruwe(chi2.dr20,nu2))
    ruwe3 <- c(ruwe3,calc.ruwe(chi2.dr30,nu3))
}
###add norminal
for(n in ns){
    if(n=='astrometric_matched_observations'){
        q2 <- length(out$cat.ind[[1]])
        q3 <- length(out$cat.ind[[2]])
    }
    if(n=='astrometric_chi2_al'){
        q2 <- mean(chi2.dr2)
        q3 <- mean(chi2.dr3)
        if(additive){
            q3 <- sqrt(q3^2 + 1.4^2*nu2^2)
        }
    }
    if(n=='astrometric_gof_al'){
        q2 <- calc.gof(mean(chi2.dr2),nu2)
        q3 <- calc.gof(mean(chi2.dr3),nu3)
        if(additive){
#            q3 <- q3 + 1.4^2*nu2^2
            q3 <- q3
        }
    }
    if(n=='astrometric_excess_noise'){
        q2 <- mean(jitter2)
        q3 <- mean(jitter3)
    }
    if(n=='ruwe'){
        q2 <- mean(ruwe2)
        q3 <- mean(ruwe3)
#        if(additive) q2 <- sqrt(q2^2+1.4^2)
#        if(additive) q3 <- sqrt(q3^2+1.4^2)
    }
    if(n=='astrometric_excess_noise_sig'){
        q2 <- q3 <- NA
    }
    if(n=='ruwe'){
###Lindegren+2018
        ind <- which.min(sqrt((dr2['phot_g_mean_mag']-u0.table[,'g_mag'])^2 +(dr2[,'bp_rp']-u0.table[,'bp_rp'])^2))
        u0 <- u0.table[ind,3]
        p2 <- sqrt(dr2[,'astrometric_chi2_al']/(dr2[,'astrometric_n_good_obs_al']-5))/u0
#        ind <- which.min(sqrt((dr3['phot_g_mean_mag']-u0.table[,'g_mag'])^2 +(dr3[,'bp_rp']-u0.table[,'bp_rp'])^2))
#        u0 <- u0.table[ind,3]
#        p3 <- sqrt(dr3[,'astrometric_chi2_al']/(dr3[,'astrometric_n_good_obs_al']-5))/u0
        p3 <- dr3[,n]
    }else if(n=='astrometric_matched_observations'){
        p2 <- dr2[,'astrometric_matched_observations']
        p3 <- dr3[,'astrometric_matched_transits']
    }else{
        p2 <- dr2[,n]
        p3 <- dr3[,n]
    }
    cat('GDR2',n,':',p2,';Predicted:',q2,'\n')
    cat('GDR3',n,':',p3,';Predicted:',q3,'\n')
}
