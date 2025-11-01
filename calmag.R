source('timing_function.R')
source('mass_luminosity.R')
library(magicaxis)
dir <- '/Users/ffeng/Documents/projects/dwarfs/rvastro/'
file.age <- list.files(paste0(dir,'ATMO_2020_cooling_tracks_age/'),pattern='txt$',full.name=TRUE)
file.age2 <- list.files(paste0(dir,'ATMO_2020_cooling_tracks_age_v2/'),pattern='txt$',full.name=TRUE)
file.mass <- list.files(paste0(dir,'ATMO_2020_cooling_tracks_mass/'),pattern='txt$',full.name=TRUE)
host <- read.table('host_par.txt',header=TRUE)
info <- read.csv('../data/code/SigClassify/common_withSig1.csv',sep=',')

#starage <- read.csv('starage_0.1_to_1asec.csv')
#tab <- read.table('bd_list.txt',header=TRUE)
tab <- read.table('companion_m5_120.txt',header=TRUE)
tmp <- c()
cc <- colnames(tab)
epoch <- sum(time_Yr2jd(2023))
ss <- data.frame()
#ss <- c()
for(j in 1:nrow(tab)){
    ii <- which(tab[j,c('')]>5 & tab[j,indm]<75 & tab[j,c('p','e_p','E_p')]>1000/365.25)
    if(length(ii)>0){
        star <- tab[j,1]
        ind1 <- which(info[,1]==star)
###get the separation at 2023 and the closest and widest separations for each companion
        for(i in ii){
            if(length(ind1)>0){
                lum <- as.numeric(info[j,'lum_val'])
                if(is.na(lum)){
                    cat('mlr!\n')
                    lum <- m2l.eker15(tab[j,'mstar'])
                }
            }else{
                cat('mlr!\n')
                lum <- m2l.eker15(tab[j,'mstar'])
            }
            ss <- rbind(ss,data.frame(tab[j,1],tab[j,'mstar'],tab[j,'emstar'],as.numeric(lum),as.numeric(tab[j,'plx'])))
            alphas <- tab[j,c('a','e_a','E_a')]*tab[j,'plx']#mas
            tmp <- rbind(tmp,c(as.numeric(tab[j,indm[i]+0:2]),as.numeric(tab[j,c('p','e_p','E_p')]),as.numeric(alphas)))
        }
    }
}
tmp <- data.frame(ss,tmp)
cn <- c('star','mstar','emstar','lum','plx','mpJ.opt','mpJ.lower','mpJ.upper','per.opt','per.lower','per.upper','rho.opt','rho.lower','rho.upper')
colnames(tmp) <- cn
ii <- which(tmp[,'rho.opt']>100)
tmp <- tmp[ii,]

star <- tmp[,'star']
mp <- tmp[,'mpJ.opt']
mstar <- tmp[,'mstar']
rho <- tmp[,'rho.opt']
lum <- tmp[,'lum']
#ii <- which(is.na(lum))
#lum[ii] <- unlist(sapply(ii,function(i)
ages <- c(0.1,1,10)
edr3 <- read.csv('~/Documents/projects/dwarfs/data/combined/star_EDR3_info.csv',sep='|')
###load and add missing entries
simbad <- read.csv('~/Documents/projects/dwarfs/data/combined/star_simbad_info.csv',sep='|')
simbad <- rbind(simbad[1,],simbad)
simbad[1,c(2:3,5:10,13:19)] <- c('WASP-8','WASP-8','359.9002966060205', '-35.0313676342542','109.747', '7.614','11.0873','-1.94','10.49','9.87','9.91','9.6125','8.501','8.218','8.086')
##
sn <- tmp[,'star']
sn[sn=='HD115404A'] <- 'HD115404'
ii <- match(sn,gsub(' ','',simbad[,2]))
simbad <- simbad[ii,]
contrastJ <- contrastH <- contrastK <- contrastL <- c()
for(age in ages){
    if(age==1){
        dat <- read.table(paste0(dir,'ATMO_2020_cooling_tracks_age/1.0Gyr_ATMO_CEQ_vega.txt'),header=TRUE)
    }else if(age==5){
        dat <- read.table(paste0(dir,'ATMO_2020_cooling_tracks_age/5.0Gyr_ATMO_CEQ_vega.txt'),header=TRUE)
    }else if(age==10){
        dat <- read.table(paste0(dir,'ATMO_2020_cooling_tracks_age/10.0Gyr_ATMO_CEQ_vega.txt'),header=TRUE)
    }else if(age==0.1){
        dat <- read.table(paste0(dir,'ATMO_2020_cooling_tracks_age_v2/0.1Gyr_ATMO_CEQ_vega.txt'),header=TRUE)
    }
    mm <- dat[,'Mass']*1047.348644
    teff0 <- dat[,'Teff']
    loglum0 <- dat[,'Luminosity']#base10
    Jmag0 <- dat[,'MKO_J']
    Hmag0 <- dat[,'MKO_H']
    Kmag0 <- dat[,'MKO_K']
    if(min(dat[,'Mass'])*1047>min(mp)){
        mm <- c(mm,min(mp))
        teff0 <- c(teff0,teff0[1])
        loglum0 <- c(loglum0,loglum0[1])
        Jmag0 <- c(Jmag0,Jmag0[1])
        Hmag0 <- c(Hmag0,Hmag0[1])
        Kmag0 <- c(Kmag0,Kmag0[1])
    }
    Tfun <- approxfun(mm,teff0)
    Lfun <- approxfun(mm,loglum0)
    Jfun <- approxfun(mm,Jmag0)
    Hfun <- approxfun(mm,Hmag0)
    Kfun <- approxfun(mm,Kmag0)
    Jpmag <- Jfun(mp)
    Hpmag <- Hfun(mp)
    Kpmag <- Kfun(mp)
    logLp <- Lfun(mp)
    d <- 1000/tmp[,'plx']#pc
    Jsmag <- as.numeric(simbad[,'Mag.J'])-2.5*log10((d/10)^2)#log10
    Hsmag <- as.numeric(simbad[,'Mag.H'])-2.5*log10((d/10)^2)
    Ksmag <- as.numeric(simbad[,'Mag.K'])-2.5*log10((d/10)^2)
    ind <- which(is.na(Jsmag) | is.na(Hsmag) | is.na(Ksmag))
    if(length(ind)>0){
        mags <- c(Jsmag[ind],Hsmag[ind],Ksmag[ind])
        i3 <- (1:length(star))[-ind]
        i2 <- i3[which(abs(tmp[i3,'mstar']-tmp[ind,'mstar'])<0.1 & abs(tmp[i3,'lum']-tmp[ind,'lum'])<0.1 & !is.na(Jsmag[i3]) & !is.na(Hsmag[i3]) & !is.na(Ksmag[i3]))]
        i2 <- i2[which.min(abs(tmp[i2,'lum']-tmp[ind,'lum']))]
        if(is.na(Jsmag[ind])) Jsmag[ind] <- Jsmag[i2]
        if(is.na(Hsmag[ind])) Hsmag[ind] <- Hsmag[i2]
        if(is.na(Ksmag[ind])) Ksmag[ind] <- Ksmag[i2]
    }
    ii <- which(star=='HD74014')
    cat('Jpmag-Jsmag=',Jpmag[ii]-Jsmag[ii],'\n')
    contrastJ <- cbind(contrastJ,10^((Jsmag-Jpmag)/2.5))
    contrastH <- cbind(contrastH,10^((Hsmag-Hpmag)/2.5))
    contrastK <- cbind(contrastK,10^((Ksmag-Kpmag)/2.5))
    contrastL <- cbind(contrastL,10^(logLp-log10(lum)))
}
out <- data.frame(tmp,rho,contrastJ,contrastH,contrastK,contrastL)
colnames(out) <- c(colnames(tmp),'rho.mas',paste0('ContrastJ',ages),paste0('ContrastH',ages),paste0('ContrastK',ages),paste0('ContrastL',ages))
#index <- which(tab[,'Mcomp.Mjup']<100 & tab[,'Mcomp.Mjup']<100)
index <- 1:nrow(out)
fout1 <- 'targets_for_Subaru.txt'
cat(fout1,'\n')
write.table(out[index,],file=fout1,quote=FALSE,row.names=FALSE)

#fout2 <- 'targets_for_Subaru_allinfo.txt'
fout2 <- 'test.txt'
cat(fout2,'\n')
out.all <- cbind(out,tmp)[index,]
write.table(out.all,file=fout2,quote=FALSE,row.names=FALSE)

xlab <- expression(rho*' [mas]')
cc <- c('black','red','blue','green')
fout <- paste0('contrast_sep_companion.pdf')
cat(fout,'\n')
pdf(fout,13,4)
size <- 1.2
par(mfrow=c(1,3),oma=c(0,0,0,5),mar=c(4,4,1.5,1),cex=size,cex.lab=size)
for(age in ages){
    bb <- paste0(c('ContrastJ','ContrastH','ContrastK','ContrastL'),age)
    for(i in 1:length(bb)){
        if(i==1){
            plot(out[,'rho.mas'],out[,bb[i]],xlab=xlab,ylab='Contrast',main=paste(age,'Gyr'),col=cc[1],pch=1,ylim=range(out[,bb]),log='xy',yaxt='n')
            magaxis(side=2,majorn=5)
        }else{
            points(out[,'rho.mas'],out[,bb[i]],col=cc[i],pch=1)
        }
    }
    if(age==ages[length(ages)]){
        legend('topright',inset=c(-0.5,0),legend=c('J band','H band','K band','Luminosity'),col=cc,bty='n',xpd=NA,pch=1)
    }
}
dev.off()

fout <- paste0('contrast_mp_companion.pdf')
cat(fout,'\n')
pdf(fout,13,4)
xlab <- expression(m[c]*' ['*m[Jup]*']')
par(mfrow=c(1,3),oma=c(0,0,0,5),mar=c(4,4,1.5,1),cex=size,cex.lab=size)
ind <- which(out[,'rho.mas']>100)
for(age in ages){
    bb <- paste0(c('ContrastJ','ContrastH','ContrastK','ContrastL'),age)
    for(i in 1:length(bb)){
        if(i==1){
            plot(mp[ind],out[ind,bb[i]],xlab=xlab,ylab='Contrast',main=paste(age,'Gyr'),col=cc[1],pch=1,ylim=range(out[,bb]),log='xy',yaxt='n')
            magaxis(side=2,majorn=5)
        }else{
            points(mp[ind],out[ind,bb[i]],col=cc[i],pch=1)
        }
    }
    if(age==ages[length(ages)]){
        legend('topright',inset=c(-0.5,0),legend=c('J band','H band','K band','Luminosity'),col=cc,bty='n',xpd=NA,pch=1)
    }
}
dev.off()

source('contrast_table.R')
