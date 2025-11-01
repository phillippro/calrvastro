options(scipen = 3)
Args <- commandArgs(trailingOnly=TRUE)
cat('Args=',Args,'\n')
if(length(Args)!=0){
    ids <- as.character(Args[1])
    Niter <- as.integer(Args[2])
    Nbin.per <- as.numeric(Args[3])#frac.burn <- as.numeric(Args[3])
    nbin.per <- as.numeric(Args[4])#frac.adapt <- as.numeric(Args[4])
    tem <- as.numeric(Args[5])
    inicov <- as.numeric(Args[6])
    Pini <- as.integer(Args[7])
    noise.models <- as.character(Args[8])
    period.par <- as.character(Args[9])
    Ncores <- as.integer(Args[10])
    Np <- as.integer(Args[11])
    mode <- as.integer(Args[12])
    Nw <- as.integer(Args[13])
    prior.type <- as.character(Args[14])
    calibrations <- as.numeric(Args[15])
    if(Nw>1){
        ids <- c(ids,Args[16:(15+(Nw-1))])
        noise.models <- c(noise.models,Args[(16+Nw-1):(15+2*(Nw-1))])
    }
}else{
    ids <- 'HD128621_season1_53413to54400_ind0'
    Niter <- 1.0e3
    Nbin.per <- 1
    nbin.per <- 1
    tem <- 0.1
    inicov <- 1e-3
    Pini <- 200#day
    noise.models <- 'ARMA01'#noise.model: white, GP(R), ARMA, TJ(Ntj=1,noise vary with RHK or SA index, the third column of HARPS data), TJ(Ntj=3,vary with FWHM, BIS, RHK), TARMA, TGP, ARMATJ(ARMA+TJ), GPTJ,TJAR(model the RV contributed by index as a AR(p)-like model), ARMATJAR(ARMA+TJAR), PSID (previous-subsequent index dependent, this model is similar to TJAR but without time-varying/index-dependent jitter), ARMAPSID(ARMA+PSID)
    period.par <- 'logP'
    Ncores <- 4
    Np <- 0
    mode <- 11#data,sim
    Nw <- 7#fit models to multiple wavelength data sets simultaneously
    prior.type <- 'mt'
    calibrations <- 0
    if(Nw>1){
#	 ids <- c(ids,'HD128621_season2_54400to54760_ind0.6.3.5.2.20.17.1.18','HD128621_season3_54760to55120_ind0.3.6.2.1.5.7.10.14.4','HD128621_season4_55120to55480_ind0.6.5.2.3.13.15.21.1.14','HD128621_season5_55480to55840_ind0.6.3.11.18.8.7.21.20','HD128621_season6_55840to56200_ind0.6.21.19.13.15.12.3.5.9','HD128621_season7_56200to56560_ind0.9.1.7.14.8.6.2')
	 ids <- c(ids,'HD128621_season2_54400to54760_ind0','HD128621_season3_54760to55120_ind0','HD128621_season4_55120to55480_ind0','HD128621_season5_55480to55840_ind0','HD128621_season6_55840to56200_ind0','HD128621_season7_56200to56560_ind0')
#        noise.models <- c(noise.models,'ARMA02','ARMA06','ARMA06','ARMA05','ARMA05','ARMA07')
        noise.models <- c(noise.models,'ARMA01','ARMA01','ARMA01','ARMA01','ARMA01','ARMA01')
    }
}
#polynomial orders
if(exists('Npoly0') & exists('Npoly.sub0')){
   Npoly <- Npoly0
   Npoly.sub <- Npoly.sub0
}else{
   Npoly <- floor(mode/10)
   Npoly.sub <- mode%%10
}
detrend <- FALSE
#detrend.type <- 'qp'
#detrend.type <- 'sin'
detrend.type <- 'poly'
detrend.plot <- FALSE
##########################################
##############part I: global parameters
##########################################
##############check if the parameters can be passed from previous runs
if(exists('calibrations0')){
    calibrations <- calibrations0
    rm(calibrations0)
}
if(exists('noise.models0')){
    noise.models <- noise.models0
    rm(noise.models0)
}
cat('calibrations=',calibrations,'\n')
if(exists('Niter1')){
    Niter <- Niter1
    rm(Niter1)
}
if(exists('ids0')){
    ids <- ids0
    rm(ids0)
}
if(exists('Nw0')){
    Nw <- Nw0
    rm(Nw0)
}
if(exists('Ncores0')){
    Ncores <- Ncores0
    rm(Ncores0)
}
if(exists('Np1')){
    Np <- Np1
}
if(exists('noise.models0')){
    noise.models <- noise.models0
    rm(noise.model0)
}
time.unit <- 365.24
frac.burn <- 0.5
frac.adapt <- 1
index.more <- FALSE
dRV.mode <- 'force'#center; seq
#independence <- TRUE#assuming that the Keplerian signal, noise components in different aperture data sets are independent
independence <- FALSE
if(!exists('quantify')){
    quantify <- FALSE
}
#if(Np==0) quantify <- TRUE
workflow <- TRUE
tree <- FALSE
tree.burn <- FALSE
#chain.type <- 'parallel'#'section', 'adapt', 'parallel','normal'
sim.data <- FALSE
data.mode <- 'single'#'combined'#combined different data sets for one target
scaled <- TRUE
Nupdate <-  1#adaptive frequency
Nburn = round(Niter/Ncores*(frac.burn))#burning fraction
n0 = max(round(Niter/Ncores*(1-frac.adapt)+2),10)#where the adaptive cov. used#at least start from the second interation
#n0 = Niter
eps = 1e-8##optional covariance parameter to avoid negative determinate
Ntry = 1000#the times limit to try 
if(exists('prior.type0')){
    prior.type <- prior.type0
    rm(prior.type0)
}
cat('prior.type=',prior.type,'\n')
verbose <- FALSE
#verbose <- TRUE
Nverbose <- round(Niter/10)#frequency of showing
#the parameters are stored in the following list object
par.data <- list()
par.data <- lapply(1:Nw,function(i) par.data[[ids[i]]] <- list())
cat('ids=',ids,'\n')
cat('length(par.data)=',length(par.data),'\n')
names(par.data) <- ids
#############################################################################
###the loop to fill parameters or data into the list object
#############################################################################
ins <- c('SOPHIE','HARPS','KECK','ELODIE','HIRES','HIRES','HARPN','PFS','APF','UVES','LICK','PRD')
for(jj in 1:Nw){
    id <- ids[jj]
    noise.model <- noise.models[jj]
##########
    id2 <- id
    if(grepl('wbin',id)) id2 <- gsub('wbin1hr_jitter0_','',id)
    star.name <- target <- gsub('_.+','',id2)
    if(grepl('KD2_3',id2)) star.name <- target <- 'K2-3'
    if(grepl('MP_Mus',id)){
        star.name <- target <- paste0(target,'_Mus')
    }
    ID <- gsub('_ind.+','',id)
    ID <- gsub('\\d+ep.+','',ID)
    ID <- gsub('\\d+ep\\d+','',ID)
    ID <- gsub('\\d+ep','',ID)
    ID <- gsub('ro\\d','',ID)
    ID <- gsub('ro','',ID)
    if(grepl('bin',id)){
        ID <- gsub('_\\dbin.+','',ID)
    }
####instrument
    instrument <- 'HARPS'
    if(grepl('_',ID)) instrument <- gsub('.+_','',ID)
#    if(grepl(paste(ins,collapse='|'), ID)){
#        for(j in 1:length(ins)){
#            if(grepl(ins[j],ID)){
#                instrument <- ins[j]
#            }
#        }
#        target <- paste0(target,'_',instrument)
#        ID <- gsub(paste0('_',instrument,'.*'),'',ID)
#    }
    cat('instrument=',instrument,'\n')
####Inds
    f1 <- gsub('.+ind','',id)
    f2 <- gsub('_.+','',f1)
    cat('f2=',f2,'\n')
    if(grepl('\\.',f2)){
        Inds <- as.integer(unlist(strsplit(f2,'\\.')))
    }else{
        Inds <- as.integer(unlist(strsplit(f2,'')))
    }
    Inds <- Inds[which(Inds>0)]
    cat('Inds=',Inds,'\n')
####divide the data into sevaral epoch ranges (or chunks) and then apply independent noise models to each. The independent parameters are determined by ep.par. 
#ep.all <- c('a','b','s','ma','c1','c2','c3','d2-1','d3-2','d4-3','d5-4','d6-5','dc3-2','dc2-1')
    if(all(Inds==0)){
        ep.all <- c('s','ma')
    }else{
        ep.all <- c('s','ma',paste0('c',1:length(Inds)))
    }
    if(grepl('ep',id)){
        f1 <- gsub('.+_','',id)
        par.epo <- as.integer(gsub('ep.*','',f1))
        if(par.epo>50){
            epoch.type <- 'gap'
            dt.tvn <- par.epo
        }else{
            epoch.type <- 'duration'
            Nepoch <- par.epo
        }
        f1 <- gsub('.+ep','',id)
        f1 <- gsub('_.+','',f1)
        nepoch <- as.integer(gsub('par.*','',f1))
        if(nepoch==0){
            nepoch <- 1:Nepoch
        }
###### if length(nepoch)==1, the epoch groups are analyzed independently.
###### Otherise, they are analyzed in combination, but with independent noise parameters. 
        if(grepl('eppar|ep\\dpar',id)){
            f1 <- gsub('_con|_ccf','',id)
            ep.par <- c(gsub('.+par','',f1))
            cat('ep.par=',ep.par,'\n')
            if(ep.par==''){
                ep.par <- ep.all
            }
        }else if(grepl('\\dep',id) & !grepl('ep\\d',id)){
            ep.par <- ep.all
        }else{
            ep.par <- ''
        }
    }else{
        ep.par <- ''
        epoch.type <- 'duration'
        nepoch <- Nepoch <- 1
    }
    cat('id=',id,'\n')
    if(grepl('GP',noise.model) & nchar(noise.model)>2){
	gp.type <- substring(noise.model,3,nchar(noise.model))
	noise.model <- substring(noise.model,1,2)
    }
    p <- 0
    q <- 0
    cat('noise.model=',noise.model,'\n')
    if(grepl('ARMA',noise.model)){
        if(grepl('\\d',noise.model)){
	    pq <- as.integer(gsub('ARMA','',noise.model))
	    if(pq>=10 & grepl('ARMA0',noise.model)){
	        p <- floor(pq/100)
		q <- pq%%100
	    }else{
	        p <- floor(pq/10)
		q <- pq%%10
	    }
#            p <- as.integer(substring(noise.model,nchar(noise.model)-1,nchar(noise.model)-1))
#            q <- as.integer(substring(noise.model,nchar(noise.model),nchar(noise.model)))
            noise.model <- sub('p\\dq\\d+','',noise.model)
        }
        if(nchar(noise.model)>6){
            arma.type <- substring(noise.model,6,nchar(noise.model))
        }
	noise.model <- gsub('\\d+','',noise.model)#substring(noise.model,1,4)
    }
    if(p==0 & q==0 & grepl('ARMA',noise.model)){
        noise.model <- 'white'
    }

    data.type <- ''
    cat('ID=',ID,'\n')
    if(!exists('AB')){
        AB <- FALSE#whether or not model the RV(Aindex) and RV(Bindex)
    }
    if(!exists('TJAR.mode')){
        if(floor(Niter/100000)==19){
            TJAR.mode <- 'jitter'#jitter or RV; model the jitter/RV dependence on activity index
        }else if(floor(Niter/100000)==18){
            TJAR.mode <- 'RV'
        }else{
            TJAR.mode <- 'RV'
        }
    }
    if(noise.model=='PSID' | noise.model=='ARMAPSID'){
        TJAR.mode <- 'RV'
    }
    TJAR.AB <- PSID.AB <- FALSE
    if(!exists('kep.type0')){
        kep.type <- 'pure'#pure or AK(apodized Keplerian) or QP (quasi-periodic) or AKQP (both apodized and quasi-periodic)
    }else{
        kep.type <- kep.type0
        rm(kep.type0)
    }
    TJAR.sym <- 'sym'#past or sym or future or asym, which define how the RV component depend on the Sindex.
    fix.par <- ''
    if(!exists('fixP')){
        fix.par <- 'NA'
    	if(any(grepl('per',fix.par))){
	     fixP <- TRUE
        }else{
             fixP <- FALSE
        }
    }
    arma.type <- 'abs'#squared, vary, abs, quasi-periodic 
    if(!exists('gp.type')){
	gp.type <- 'abs'#sq, abs, qp
    }
    if(noise.model=='TARMA' ){
        arma.type <- 'index'
                                        #    arma.type <- 'index.exponent'
    }
    Par <- 0
    if(noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='PSID' | noise.model=='ARMAPSID'){
        Par <- 1
    }
    ep.name <- ''
    if(grepl('ep',id)){
        m <- regexpr('\\d+ep\\d+',id)
        ep.name <- regmatches(id,m)
        cat('ep.name=',ep.name,'\n')
    }
    cat('q=',q,'\n\n')
####print input parameters
############################################
###part II: read data 
############################################
    folder <- paste0('../data/aperture/',star.name,'/')
    if(!file.exists(folder))     folder <- paste0('../data/aperture/',target,'/')
    f <- paste0(folder,ID,'.dat')
    if(!file.exists(f)){
       if(instrument!='KECK'){
          f <- paste0(folder,target,'_',instrument,data.type,'.dat')
       }else{
          f <- paste0(folder,target,'_',instrument,data.type,'.vels')
       }
    }
    cat('f=',f,'\n')
    if(!file.exists(f)){
#       folder <- paste0('../data/aperture/',ID,'/')
       f <- paste0(folder,ID,'.dat')
    }
    cat('folder=',folder,'\n')
    cat('f=',f,'\n')
    data <- read.table(f,nrows=1)
    if(class(data[1,1])=='factor'){
        data <- read.table(f,header=TRUE,check.names=FALSE)
    }else{
        data <- read.table(f)
    }
#    data <- read.table(f,header=TRUE,check.names=FALSE)
    ferror <- gsub('.dat$','_error.dat',f)
    if(ncol(data)>3 & grepl('HARPS',id) & file.exists(ferror)){
        eproxy <- read.table(ferror,header=TRUE,check.names=FALSE)
    }else{
        eproxy <- c()
    }
    times0 <- trv <- data[,1]%%2400000
    ind <- sort(trv,index.return=TRUE)$ix
    trv <- trv[ind]
    RV0 <- data[ind,2]
    if(Npoly.sub>1){
       val <- lm(RV0 ~ poly((trv-min(trv))/time.unit,Npoly.sub))
       RV <- residuals(val)
    }else{
       RV <- RV0
    }
    eRV <- data[ind,3]
    if(!all(Inds==0)){
        proxies <- data[,3+Inds]
	if(!is.null(eproxy)){
            eproxies <- eproxy[,Inds]
	}else{
            eproxies <- proxies
	    if(is.null(dim(proxies))){
                 eproxies <- rep(1,length(proxies))
	    }else{
                 eproxies[1:nrow(proxies),1:ncol(proxies)] <- 1
            }
	}
        proxy.names <- gsub('.proxy.+','',colnames(data)[3+Inds])
        if(length(Inds)==1){
            proxies <- matrix(proxies,ncol=1)
            eproxies <- matrix(eproxies,ncol=1)
        }
        for(j in ncol(proxies)){
###normalization
            eproxies[,j] <- eproxies[,j]/sd(proxies[,j])
            proxies[,j] <- scale(proxies[,j])
        }
        if(is.null(dim(proxies))){
            proxies <- matrix(proxies,ncol=1)
        }
	if(detrend){
            qp <- function(t,a,b,T,tau,c,d){
                a*sin(2*pi*t/T+b*sin(2*pi*t/tau)+c)+d
            }
            sinusoid <- function(t,a,T,c,d){
                a*sin(2*pi*t/T+c)+d
            }
            for(j in 1:length(Inds)){
                proxy <- proxies[,j]
                if(detrend.type=='qp'){
                    val <- nls(proxy~qp(trv,a,b,T,tau,c),start=list(a=sd(proxy),b=0.1,T=max(trv)-min(trv),tau=max(trv)-min(trv), c=0,d=0),control=nls.lm.control(maxiter=500))
                }else if(detrend.type=='sin'){
                    val <- nls(proxy~sinusoid(trv,a,T,c,d),start=list(a=sd(proxy),T=max(trv)-min(trv),c=0,d=0),control=nls.lm.control(maxiter=500))
                }else{
                    x <- trv-min(trv)
                    val <- lm(proxy ~ poly(x,degree=5,raw=TRUE))
                }
                co <- coef(val) 
                if(detrend.plot){
                    plot(trv,proxy,xlab='time',ylab='new proxy',ylim=1.1*range(proxy))
                    legend('topright',legend=paste0('r0=',format(cor(RV,proxy),digit=3)),bty='n')
                    if(detrend.type=='qp'){
                        curve(qp(x,a=co['a'],b=co['b'],T=co['T'],tau=co['tau'],c=co['c'],d=co['d']),add=TRUE ,lwd=2, col="steelblue")
                    }else if(detrend.type=='sin'){
                        curve(sinusoid(x,a=co['a'],T=co['T'],c=co['c'],d=co['d']),add=TRUE ,lwd=2, col="steelblue")
                    }else{
                        xx <- seq(0,max(trv)-min(trv),length.out=1000)
                        lines(xx+min(trv), predict(val, data.frame(x=xx)), col="steelblue")
                    }
                    points(trv,residuals(val),col='green')
                    legend('topleft',legend=paste0('r1=',format(cor(RV,residuals(val)),digit=3)),bty='n')
                    Sys.sleep(3)
                }
                proxy <- residuals(val)
                proxies[,j] <- proxy
            }
	}
    }else{
        proxies <- NA
        eproxies <- NA
        proxy.names <- NA
    }
###epochs
    epochs <- c()#the epochs used to divide the data into chunks
    Nes <- c()
    ind.remain <- c()
    if(epoch.type=='gap'){
            tmp <- trv
            ind1 <- which(diff(tmp)>dt.tvn)
            ind2 <- ind1+1
            val <- c(min(tmp)-1,(trv[ind1,1]+trv[ind2,1])/2,max(tmp)+1)
            epochs <- combine.index(epochs,val)
            Nes <- c(Nes,length(val)-1)
            if(length(nepoch)<Nepoch & Nepoch>1){
                ind.keep <- which(trv>val[nepoch] & trv<val[nepoch+1])
                ind.remain <- combine.index(ind.remain,ind.keep)	   
                trv[-ind.keep] <- NA
                RV[-ind.keep] <- NA
                eRV[-ind.keep] <- NA
                if(!all(is.na(proxies))){
                    proxies[-ind.keep,] <- NA
                    eproxies[-ind.keep,] <- NA
                }
            }
    }else{
            tmp <- trv
            dur <- (max(tmp)+10-min(tmp))/Nepoch#avoid the maximum time point
            val <- min(tmp)-1
            for(k in 1:Nepoch){
                val <- c(val,min(tmp)+k*dur)
            }
            if(length(nepoch)==1){
#                epochs <- matrix(val,ncol=1)
                epochs <- val
                Nes <- Nepoch
            }else{
                epochs <- combine.index(epochs,val)
                Nes <- c(Nes,Nepoch)
            }            
            if(length(nepoch)<Nepoch & Nepoch>1 & length(nepoch)==1){
                ind.keep <- which(trv>val[nepoch] & trv<val[nepoch+1])
                ind.remain <- ind.keep
                trv <- trv[ind.keep]
                times0 <- times0[ind.keep]
                RV <- RV[ind.keep]
                eRV <- eRV[ind.keep]
                if(!all(is.na(proxies))){
                    proxies <- proxies[ind.kepp,]
                    eproxies <- eproxies[ind.kepp,]
                }
            }
    }
    if(is.null(dim(proxies)) & !all(is.na(proxies))){
        proxies <- matrix(proxies,ncol=1)
        eproxies <- matrix(proxies,ncol=1)
    }
#    tnorm <- poly((trv-min(trv))/time.unit,degree=Npoly)[,1:Npoly,drop=FALSE]
####save all data and parameters for one data set
    vars <- unique(c('ID','trv','t0','times0','RV','eRV','mom.name','mom','epoch.type','epochs','Nes','Nepoch','nepoch','ep.par','ep.name','proxies','proxy.names','eproxies','Inds','noise.model','p','q','data.type','Nd','Ntj','gp.type','folder','arma.type','instrument','target','star.name','data','cmax','cmin','cini'))
    for(k in 1:length(vars)){
        if(exists(vars[k])){
            par.data[[id]][[vars[k]]] <- eval(parse(text = vars[k]))
####remove variables of previous run
            rm(list=vars[k])
        }
    }
}
###some parameters derived from the whole data sets
###time range
trv.all <- c()
RV.all <- c()
eRV.all <- c()
qs <- c()
ps <- c()
for(i3 in 1:Nw){
    var <- names(par.data[[ids[i3]]])
    for(k in 1:length(var)){
        assign(var[k],par.data[[ids[i3]]][[var[k]]])
    }
    trv.all <- c(trv.all,trv)
    RV.all <- c(RV.all,RV)
    eRV.all <- c(eRV.all,eRV)
    qs <- c(qs,q)
    ps <- c(ps,p)
}
Ndata <- length(trv.all)
rmin <- min(RV.all)
rmax <- max(RV.all)
tmin <- min(trv.all)
tmax <- max(trv.all)
#####for very significant signals
err.rel <- mean(eRV.all)/(rmax-rmin)
if(err.rel<1e-2){
    for(i3 in 1:Nw){
        par.data[[i3]]$erv <- par.data[[i3]]$eRV/(err.rel*10)
    }
    small.err <- TRUE
}else{
    small.err <- FALSE
}
######talk type for gen_mcmc.R
#talk.type <- 'like'
talk.type <- 'post'
if(Ndata>500 & median(diff(trv.all))<1){
   talk.type <- 'post'
}
####cmax
for(j in 1:Nw){
    var <- names(par.data[[ids[j]]])
    for(k in 1:length(var)){
        assign(var[k],par.data[[ids[j]]][[var[k]]])
    }
    Kmax = 2*(rmax-rmin)
    Kmin = 0#-Kmax
    cini <- cmax <- cmin <- c()
    if(length(Inds)>0){
        for(i in 1:ncol(proxies)){
            tmp <- Kmax/(max(proxies[,i])-min(proxies[,i]))
            cmax <- c(cmax,tmp)
            cmin <- c(cmin,-tmp)
            cini <- c(cini,0)
        }
    }
    par.data[[ids[j]]]$cmax <- cmax
    par.data[[ids[j]]]$cmin <- cmin
    par.data[[ids[j]]]$cini <- cini
}
