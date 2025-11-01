##################################################
####periodograms for activity indicators and plot
##################################################
for(i in ins){
    acts <- names(out[[i]][['activity']])
    for(act in acts){
        indice <- out[[i]][['activity']][[act]]
        if(!grepl('AP',act)){
            if(!grepl('window',act)){
                cat('calulate BFP for',act,'for',i,'set of',target,'\n')
                per <- BFP(indice[,1],indice[,2],indice[,3],Nma=1,Nar=0,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=FALSE,noise.only=FALSE,Nsamp=Nsamp)
            }else{
                cat('calulate LS for',act,'for',i,'RV set of',target,'\n')
                per <- lsp(times=indice[,1],x=indice[,2],ofac=ofac,from=fmin,to=fmax,alpha=c(0.1,0.01,0.001))
            }
            out[[i]][['BFP']][[act]] <- per
        }
    }
}

###determine rotation period
GPf <- FALSE
if(any(grepl('photometry',names(out)))){
#    per <- BFP(out[['photometry']]$data[,1],out[['photometry']]$data[,2],out[['photometry']]$data[,3],Nma=1,Nar=1,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=TRUE,noise.only=FALSE,Nsamp=Nsamp,gp.par=rep(NA,3))
 per <- BFP(out[['photometry']]$data[,1],out[['photometry']]$data[,2],out[['photometry']]$data[,3],Nma=0,Nar=0,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=TRUE,gp.par=rep(NA,3),noise.only=TRUE,Nsamp=Nsamp,sampling=sampling)
    out[['photometry']][['BFP']] <- list(GP=per)
    GPf <- TRUE
}
#source('signal_likelihood.R')

##################################################
####periodograms for different noise models for RVs
##################################################
for(j in 1:Nset){
    cat('\n')
    for(k in 1:length(noise.types)){
        data <- out[[ins[j]]]$RV
        for(i in 1:min(Nmax,4)){
            Nma <- Nar <- 0
            GP <- FALSE
            gp.par <- rep(NA,3)
            if(noise.types[k]=='MA'){
                Nma <- 1
            }
            if(noise.types[k]=='AR'){
                Nar <- 1
            }
            if(grepl('GP',noise.types[k]) & GPf){
                GP <- TRUE
                par.opt <- out['photometry']$BFP$GP$par.opt
            }
            cat('calulate BFP for',i,'signal for',noise.types[k],'for',ins[j],'RV set of',target,'\n')
            per <- BFP(data[,1],data[,2],data[,3],Nma=Nma,Nar=Nar,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=GP,noise.only=FALSE,Nsamp=Nsamp, gp.par=gp.par,renew=TRUE)
            out[[ins[j]]][['BFP']][[paste0('sig',i)]][[noise.types[k]]] <- per
            data[,2] <- per$res.s
        }
    }
}
