replacePar <- FALSE
pars <- read.table(fpar)
paropt <- as.numeric(pars[,2])
ns0 <- ns <- as.character(pars[,1])
if(any(grepl('^b\\d',ns)) | any(grepl('^s\\d',ns))){
    ns[grep('^b\\d',ns)] <- paste0('b_',out$ins.rv)
    ns[grep('^s\\d',ns)] <- paste0('J_',out$ins.rv)
    for(j in 1:length(out$ins.rv)){
        i <- out$ins.rv[j]
        if(any(grepl(paste0('^w',j,'\\d'),ns0))){
            ii <- grep(paste0('^w',j,'\\d'),ns0)
            jj <- grep(paste0('^beta',j),ns0)
            ns[ii] <- paste0(gsub(paste0('w',j),'w',ns0[ii]),'_',i)
            ns[jj] <- paste0('beta_',i)
        }
    }
    ns[grep('^d.+0$',ns)] <- gsub('0','',ns[grep('^d.+0$',ns)])
}
names(paropt) <- ns
ind <- grep('jitter',ns)
if(length(ind)>0){
    ns[ind] <- gsub('jitter','logJ',ns[ind])
    paropt[ind] <- log(paropt[ind])
    names(paropt) <- ns
}
replacePar <- TRUE
nn <- intersect(names(startvalue),ns)
cat('replacing parameters by ',paste0('pars/',target,'.par:'),nn,'\n')
par.replace <- paropt[nn]
par.replace[grep('omega',nn)] <- par.replace[grep('omega',nn)]%%(2*pi)
par.replace[grep('Omega',nn)] <- par.replace[grep('Omega',nn)]%%(2*pi)
par.replace[grep('Inc',nn)] <- par.replace[grep('Inc',nn)]%%(2*pi)
par.replace[grep('Mo',nn)] <- par.replace[grep('Mo',nn)]%%(2*pi)
startvalue[nn] <- par.replace
startvalue1 <- startvalue
ll <- loglikelihood(startvalue,prediction=TRUE)
if(out$Nrv>0){
    for(jj in 1:length(out$ins.rv)){
        i <- out$ins.rv[jj]
        startvalue1[paste0('b_',i)] <- startvalue1[paste0('b_',i)]+mean(ll$res[[i]]$res2)
    }
    if(length(out$data.binary)>0){
        for(kk in 1:length(out$ins.binary)){
            i <- out$ins.binary[kk]
            inds <- out$ind.binary[[i]]
            startvalue1[paste0('b_',i)] <- startvalue1[paste0('b_',i)]+mean(ll$res[[i]]$res2)
        }
    }
}
if(any('dra'==names(startvalue))){
    startvalue1[c('dra','ddec','dplx','dpmra','dpmdec')] <- startvalue1[c('dra','ddec','dplx','dpmra','dpmdec')]-ll$res$GDR3
}
llike <- loglikelihood(startvalue)
llike1 <- loglikelihood(startvalue1)
if(llike1>llike){
    cat('update offsets!\n')
    startvalue <- startvalue1
}else{
    cat('keep original offsets!\n')
}
