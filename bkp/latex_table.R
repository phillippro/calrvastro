#targets <- gsub('\\/.+','',fs)
#targets <- stars
#ind <- match(exception,targets)
#if(reload){
###if(FALSE){
fout <- paste0('phase',version,'.Robj')
cat(fout,'\n')
#save(list=ls(all=TRUE),file=fout)
###}
cat(fout,'\n')
ns <- c('\\tablehead{\\colhead{Name}','\\colhead{Other Name}','$P$ [yr]','$K$ [m/s]','e','$\\omega$ [deg]','$M_0$ [deg]','$I$ [deg]','$\\Omega$ [deg]','$\\ln J_{\\rm hip}$','$\\ln J_{\\rm gaia}$','$\\Delta \\alpha$ [mas]','$\\Delta \\delta$ [mas]','$\\Delta \\mu_\\alpha$ [mas/yr]','$\\Delta \\mu_\\delta$ [mas/yr]','$T_p$ [JD]','$m_p$ [$M_{\\rm Jup}$]','$a$ [au]')
nss <- c()
for(jj in 1:length(ns)){
    if(jj<3){
        nss <- c(nss,ns[jj])
     }else{
        nss <- c(nss,ns[jj],'')
     }
}
ns <- nss
Nstar <- length(stars)
#table.out0 <- table.out
#table.out <- table.out0
table1 <- table.out
table2 <- table.out
table1 <- cbind(ns,table1[,-ncol(table1)],paste0(table1[,ncol(table1)],'\\\\'))
table2 <- cbind(ns,table2[,-ncol(table2)],paste0(table2[,ncol(table2)],'\\\\'))
tab.all <- cbind(ns,table.out[,-ncol(table.out)],paste0(table.out[,ncol(table.out)],'\\\\'))
tab.all[2,ncol(tab.all)] <- gsub('\\\\\\\\','}',tab.all[2,ncol(tab.all)])
tab.all[3,1] <- paste0('\\startdata ',tab.all[3,1])
tab.all[nrow(tab.all),ncol(tab.all)] <- paste0(tab.all[nrow(tab.all),ncol(tab.all)],' \\enddata')
if(grepl('_astro_',fs[1])) solution <- 'dr2'
if(grepl('_hg3_',fs[1])) solution <- 'hg3'
if(grepl('_hg2--_',fs[1])) solution <- 'hgca'
#solution <- astrotype
fout <- paste0(solution,'_table_planet.txt')
fout0 <- paste0(solution,'_table_all.txt')
fout1 <- paste0(solution,'_table_planet1.txt')
fout2 <- paste0(solution,'_table_planet2.txt')
#cat(fout,'\n')
#cat(fout0,'\n')
#cat(fout1,'\n')
#cat(fout2,'\n')
NN <- nrow(table1)
#index <- c(1:3,NN-(2:0),4:(NN-3))
index <- 1:NN
dat1 <- table1[index,]
dat1[2,ncol(dat1)] <- paste0(dat1[2,ncol(dat1)],'\\hline')
dat1[nrow(dat1),ncol(dat1)] <- paste0(dat1[nrow(dat1),ncol(dat1)],'\\hline\\hline')

dat2 <- table2[index,]
dat2[2,ncol(dat2)] <- paste0(dat2[2,ncol(dat2)],'\\hline')
dat <- rbind(dat1,dat2)
#write.table(dat,file=fout,quote=FALSE,sep='&',row.names=FALSE,col.names=FALSE)
#write.table(tab.all,file=fout0,quote=FALSE,sep='&',row.names=FALSE,col.names=FALSE)
#write.table(dat1,file=fout1,quote=FALSE,sep='&',row.names=FALSE,col.names=FALSE)
#write.table(dat2,file=fout2,quote=FALSE,sep='&',row.names=FALSE,col.names=FALSE)


