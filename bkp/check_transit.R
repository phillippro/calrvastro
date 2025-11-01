Prefs <- rep(0,10)
#if(any(grepl('per1',par.fix))){
    c22 <- read.table('../data/combined/eph_exoclock_2022.txt',header=TRUE,row.names=NULL)
    stars.c22 <- substr(c22[,1],1,nchar(c22[,1])-1)
    ps.c22 <- substr(c22[,1],nchar(c22[,1]),nchar(c22[,1]))
    ivshina <- read.table('../data/combined/eph_ivshina22.csv',sep='|',header=TRUE)
    wang <- read.table('../data/combined/eph_wang23.txt',header=TRUE,row.names=NULL)
    target1 <- target
    if(target=='HD39091') target1 <- 'pimen'
    toi <- read.table('toi_info.txt',header=TRUE)
    ftim <- list.files(paste0('../data/combined/',target),pattern='transit.+\\.tim\\d$')
    nps <- as.integer(substr(ftim,nchar(ftim),nchar(ftim)))
    if(length(nps)>0) dP <- nps
    ps <- letters[seq( from = 2, to = 10 )][nps]
    for(j in 1:length(nps)){
        np <- nps[j]
        p <- ps[j]
        ind0 <- which((toupper(gsub('-','',target))==toupper(gsub('-','',stars.c22)) | toupper(gsub('-','',target1))==toupper(gsub('-','',stars.c22))) & p==ps.c22)
        ind1 <- which((toupper(gsub('-','',target))==toupper(gsub('-','',wang[,'Type'])) | toupper(gsub('-','',target1))==toupper(gsub('-','',wang[,'Type']))) & p==wang[,'ID'])
        ind2 <- which(toupper(gsub('-','',target))==toupper(gsub('-','',ivshina[,'Sys'])) & toupper(gsub('-','',target))==toupper(gsub('-','',ivshina[,'Sys'])))
        ind3 <- which(toupper(gsub('-','',target))==toupper(gsub('-','',toi[,'Name'])) & toupper(gsub('-','',target))==toupper(gsub('-','',toi[,'Name'])))
        if(length(ind0)>0){
            Prefs[np] <- c22[ind0[1],4]
        }else if(length(ind1)>0){
            Prefs[np] <- ivshina[ind1[1],'P0']
        }else if(length(ind2)>0){
            Prefs[np] <- wang[ind2[1],'P0']
        }else if(length(ind3)>0){
            Prefs[np]  <- toi[ind3[1],'P']
        }
    }
#}
