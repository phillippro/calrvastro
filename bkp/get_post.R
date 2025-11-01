args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    f <- gsub('.pdf','.Robj',args)
    if(!grepl('results',f)) f <- paste0('results/',f)
    cat(f,'\n')
    load(f)
    mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
    lnl <- mc[,'loglike']
    lnp <- mc[,'logpost']
    cat('maximum posterior:',lnp[which.max(lnp)],'\n')
    cat('maximum likelihood:',lnl[which.max(lnl)],'\n')
}