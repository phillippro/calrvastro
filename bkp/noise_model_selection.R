source('mcmc_func.R')
library(xtable)
options(scipen=9)
folder <- '/car-data/ffeng/dwarfs/output/'
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    targets <- as.character(args[1])
    data.type <- as.character(args[2])
}else{
#    targets <- 'HD128621_season1_53413to54400'
#    targets <- 'HD128621_season2_54400to54760'
#    targets <- 'HD128621_season3_54760to55120'
#    targets <- 'HD128621_season4_55120to55480'
#    targets <- 'HD128621_season5_55480to55840'
#    targets <- 'HD128621_season6_55840to56200'
    targets <- 'HD128621_season7_56200to56560'
#    targets <- 'list'
    data.type <- ''
#    data.type <- 'con'
#    data.type <- 'ccf'
}
if(!file.exists(folder)){
    folder <- 'output/'
}
if(targets=='list'){
    targets <- list.dirs(path = folder, full.names = FALSE, recursive = TRUE)
    if(any(targets=='')){
        targets <- targets[targets!='']
    }
}
verbose <- TRUE
#data.type <- 'ccf'
#data.type <- 'NA'
qs <- 0:8
#aps <- c(0,3,6,9,18)
if(grepl('season',targets)){
    f1 <- gsub('.+season','',targets)
    Nseason <- as.integer(gsub('_.+','',f1))
}else{
    Nseason <- 0
}
if(Nseason==1){
    Inds <- c(0,1,3,5,2,6,15,7,4,9,18)
}else if(Nseason==2){
    Inds <- c(0,6,3,5,2,20,17,1,18,15,11)
}else if(Nseason==3){
    Inds <- c(0,3,6,2,1,5,7,10,14,4,16)
}else if(Nseason==4){
    Inds <- c(0,6,5,2,3,13,15,21,1,14,20)
}else if(Nseason==5){
    Inds <- c(0,6,3,11,18,8,7,21,20,9,22)
}else if(Nseason==6){
    Inds <- c(0,6,21,19,13,15,12,3,5,9,17)
}else if(Nseason==7){
    Inds <- c(0,9,1,7,14,8,6,2,10,12,5)
}
ind.name <- c()
for(j in 1:length(Inds)){
    ind.name <- c(ind.name,paste0('Ind',paste(Inds[1:j],collapse='.')))
}
logBF <- maxlogL <- array(-1e6,dim=c(length(Inds),length(qs),length(targets)),dimnames=list(ind.name, paste0('MA(',qs,')'), targets))
opt.noise <- array(data=0,dim=c(length(targets),5))
#e.g.
#keppure_priormt_poly10_Ndata3735_quantifyFALSE_1per1_Nw1_HD128621_season5_55480to55840_ind0_0planet_ARMA08_Nsamp4000000_tem1_acc20_pretemPd_negLmax11531
opt.noise[,1] <- targets
for(i1 in 1:length(targets)){
    target <- targets[i1]
    cat('\n',target,'\n')
    ID <- gsub('_.+','',target)
    folder.target <- paste0(folder,ID)
    flist <- list.files(path=folder.target,pattern=paste0('poly10.+quantifyTRUE.+',target,'.+0planet.+pdf'))
    ind <- which(grepl('keppure_priormt',flist) & grepl('_ind0',flist))
    flist <- flist[ind]
    BF.opt <- 0
    for(i2 in 1:length(qs)){
        q <- qs[i2]
        ind <- grep(paste0('ARMA0',q),flist)
        if(length(ind)>0){
        fs <- flist[ind]
        for(i3 in 1:length(Inds)){
            Ind <- Inds[1:i3]
            index <- paste(Inds[1:i3],collapse='.')
            ind <- grep(paste0('ind',index,'_'),fs)
            f <- fs[ind]
#            cat('f=',f,'\n')
            ind <- grep('c2p',f)
            if(length(ind)>0){
                f <- f[ind]
            }
            if(length(f)>0){
            Ns <- extract.Nsamp(f)
####extract maximum likelihood
            lmax <- extract.Lmax(f)
            ind <- which(Ns==max(Ns))
            if(length(ind)>1){
                ind2 <- which.max(lmax[ind])
                ind <- ind[ind2]
            }
            f <- f[ind[1]]
            if(verbose){
                cat(paste0(q,'MA','; Ind',index,'\n'))
                cat(paste0(folder.target,'/',f),'\n\n')
            }
            maxlogL[i3,i2,i1] <- lmax[ind]
            Ndata <- extract.Ndata(f)
#            Ndata <- 330
####BIC-estimated BFs
            Npar <- 0
            if(q==0){
                Npar <- length(Ind)-1
            }else{
                Npar <- q+length(Ind)
            }
            BIC <- -2*maxlogL[i3,i2,i1]+Npar*log(Ndata)
            BIC0 <- -2*maxlogL[1,1,i1]
            if(i2==1 & i3==1){
                logBF[i3,i2,i1] <- 0
            }else{
                logBF[i3,i2,i1] <- -(BIC-BIC0)/2
            }
	    thap <- 150
	    thma <- 150
            if(i2>1 & i3>1){
                if((logBF[i3,i2,i1]-max(logBF[1:(i3-1),1:i2,i1]))>log(thap) & (logBF[i3,i2,i1]-max(logBF[1:i3,1:(i2-1),i1]))>log(thma) & logBF[i3,i2,i1]>BF.opt){
                    opt.noise[i1,2:5] <- c(ind.name[i3],q,maxlogL[i3,i2,i1],Ndata)
                    BF.opt <- logBF[i3,i2,i1]
                }
            }else if(i2>1){
                if((logBF[i3,i2,i1]-max(logBF[1:i3,1:(i2-1),i1]))>log(thma) & logBF[i3,i2,i1]>BF.opt){
                    opt.noise[i1,2:5] <- c(ind.name[i3],q,maxlogL[i3,i2,i1],Ndata)                    
                    BF.opt <- logBF[i3,i2,i1]
                }
            }else if(i3>1){
                if((logBF[i3,i2,i1]-max(logBF[1:(i3-1),1:i2,i1]))>log(thap) & logBF[i3,i2,i1]>BF.opt){
                    BF.opt <- logBF[i3,i2,i1]
                    opt.noise[i1,2:5] <- c(ind.name[i3],q,maxlogL[i3,i2,i1],Ndata)                    
                }
            }else{
                opt.noise[i1,2:5] <- c(ind.name[1],0,maxlogL[i3,i2,i1],Ndata)       
            }
        }
        }
    }
    }
}
save(maxlogL,logBF,opt.noise,file='select_noise_model.Robj')
targets.none <- c()
for(i in 1:dim(maxlogL)[3]){
    for(j in 1:dim(maxlogL)[2]){
        for(k in 1:dim(maxlogL)[1]){
            if(abs(maxlogL[k,j,i])>1e5) targets.none <- rbind(targets.none,c(dimnames(maxlogL)[[3]][i],gsub('MA','',dimnames(maxlogL)[[2]][j]),gsub('Ind','',dimnames(maxlogL)[[1]][k])))
        }
    }
}
colnames(opt.noise)=c('name','Ind','q','logLmax','Ndata')
write.table(opt.noise,file='select_noise_model.txt',quote=FALSE,append=FALSE,row.names=FALSE)
write.table(targets.none,file='targets_none.txt',quote=FALSE,append=FALSE,col.names=FALSE,row.names=FALSE)
cat('save data:\n')
cat('select_noise_model.Robj\n')
cat('select_noise_model.txt\n')
cat('\nResults:\n')
cat(colnames(opt.noise),'\n')
cat(opt.noise,'\n')
#print(xtable(formatC(logBF[,,1],format='f',digit=2)),quote=FALSE)
