tab <- read.csv('mother_cross.csv',sep='|',header=TRUE)
poss <- unique(tab[,1])
out <- c()
for(pos in poss){
    inds <- which(tab[,1]==pos & !is.na(tab[,'RV']))
    if(length(inds)>0){
        ind <- inds[which.min(tab[inds,2])]
    }else{
        inds <- which(tab[,1]==pos)
        ind <- inds[which.min(tab[inds,2])]
    }
    out <- rbind(out,tab[ind,c('target_position','RV','e_RV','angle_arcsec','RA_ICRS','e_RA_ICRS','DE_ICRS','e_DE_ICRS','Source','Gmag')])
}
write.table(out,file='mother_RVs_opt.txt',quote=FALSE,row.names=FALSE)
