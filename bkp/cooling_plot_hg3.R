file.age <- list.files('ATMO_2020_cooling_tracks_age/',pattern='txt$',full.name=TRUE)
file.age2 <- list.files('ATMO_2020_cooling_tracks_age_v2/',pattern='txt$',full.name=TRUE)
file.mass <- list.files('ATMO_2020_cooling_tracks_mass/',pattern='txt$',full.name=TRUE)
ff <- '../paper/astrometry/v6/'

###paper BDs
#star <- c('HD 131664 B','GJ 777 B','GJ 779 B', 'HD 4747 B')
star <- c('HD 131664 B','HD 190360 b','HD 190406 B', 'HD 4747 B')
m <- c(128,1.8,73,68.3)
em <- c(17.6,0.2,6.1,5.9)
#a <- c(2.32,6.7,2,3.2)
T <- c(1200,149.5,1680,1700)
eTl <- c(450,26.5,170,100)
eTu <- c(800,26.5,170,100)
#L <- c(1600,1600,1700,1700)
# 55.8 2.32Gyr  ?
# 2.6 ~6.7Gyr ?
# 72.2 1-3Gyr 1510-1850K
# 65.5 1.4-5.6Gyr 1700+/-100K

fin <- paste0(ff,'cooling_age_mass_hg3_v1.pdf')
cat(fin,'\n')
pdf(fin,6,6)
#par(mfrow=c(2,2))
size <- 1.5
par(mar=c(5,5,1,1),cex.lab=size,cex.axis=size,cex.main=size)
lum <- c()
teff <- c()
mass <- c()
ages <- c()
for(f in file.age){
    tab <- read.table(f,header=TRUE)
    teff <- cbind(teff,tab[,'Teff'])
    lum <- cbind(lum,tab[,'Luminosity'])
    mass <- cbind(mass,tab[,'Mass'])
    ages <- c(ages,tab[1,'Age'])
}
mass <- mass*1047.348644#Mjup
#plot(unlist(mass),unlist(teff),col='white',xlab=expression('Mass ['*M[Jup]*']'),ylab=expression(T[eff]*' [K]'),xlim=rev(c(55,75)),ylim=c(1550,1750))
plot(unlist(mass),unlist(teff),col='white',xlab=expression('Mass ['*M[Jup]*']'),ylab=expression(T[eff]*' [K]'),xlim=rev(c(1,130)),ylim=c(100,2000))
#plot(unlist(mass),unlist(lum),col='white',xlab=expression('Mass ['*M[Jup]*']'),ylab=expression('log(L)['*L[sun]*']'),xlim=rev(c(20,80)),ylim=c(-6.5,-4))
#age.select <- c(0.1,0.2,0.4,1,2,4,10)
#age.select <- c(1,2,4,6,8,10)
age.select <- c(1,2,4,10)
ts <- seq(1200,800,length.out=length(age.select))
n <- 0
index <- sort(ages,index.return=TRUE,decreasing=TRUE)$ix
for(j in index){
    if(any(ages[j]==age.select)){
        n <- n+1
        ind <- sort(mass[,j],index.return=TRUE)$ix
        lines(mass[ind,j],teff[ind,j],lty=1,col='darkgrey')
                                        #    lines(mass[ind,j],lum[ind,j],lty=1,col='grey')
                                        #    Ntxt <- round(nrow(mass)/2)
        Ntxt <- which.min(abs(teff[,j]-ts[n]))
        text(mass[Ntxt,j],teff[Ntxt,j],pos=3,labels=paste(ages[j],'Gyr'),col='darkgrey')
    }
}
cols <- c('brown','black','green','blue','orange')
ind <- c(1,2,3,4)
points(m[ind],T[ind],col=cols[ind],pch=20)
arrows(m[ind]-em[ind],T[ind],m[ind]+em[ind],T[ind],length=0.1,angle=90,code=3,col=cols[ind],cex=0.8,pch=20)
#arrows(m[ind],T[ind]-eTl[ind],m[ind],T[ind]+eTu[ind],length=0.1,angle=90,code=3,col=cols[ind],cex=0.8,pch=20)
##add Jupiter
points(1,134,pch=15,col=cols[5])
arrows(m[1:2],T[1:2]-eTl[1:2],m[1:2],T[1:2]+eTu[1:2],length=0.1,angle=90,code=3,col=cols[1:2],cex=0.8,pch=20,lty=2)
arrows(m[3:4],T[3:4]-eTl[3:4],m[3:4],T[3:4]+eTu[3:4],length=0.1,angle=90,code=3,col=cols[3:4],cex=0.8,pch=20)
legend('topright',legend=c(star,'Jupiter','ATMO 2020'),col=c(cols,'darkgrey'),pch=c(20,20,20,20,15,NA),lty=c(rep(NA,5),1),bty='n',cex=size)
dev.off()

teff <- lum <- mass <- ages <- c()
for(f in file.age2){
    tab <- read.table(f,header=TRUE)
    teff <- cbind(teff,tab[,'Teff'])
    lum <- cbind(lum,tab[,'Luminosity'])
    mass <- cbind(mass,tab[,'Mass'])
    ages <- c(ages,tab[1,'Age'])
}

fin <- paste0(ff,'cooling_age_mass_hg3_v2.pdf')
cat(fin,'\n')
pdf(fin,6,6)
#par(mfrow=c(2,2))
par(mar=c(5,5,1,1),cex.lab=size,cex.axis=size,cex.main=size)
#par(mar=c(5,5,1,1))
lum <- c()
teff <- c()
mass <- c()
ages <- c()
for(f in file.age2){
    tab <- read.table(f,header=TRUE)
    teff <- cbind(teff,tab[,'Teff'])
    lum <- cbind(lum,tab[,'Luminosity'])
    mass <- cbind(mass,tab[,'Mass'])
    ages <- c(ages,tab[1,'Age'])
}
mass <- mass*1047.348644#Mjup
#plot(unlist(mass),unlist(teff),col='white',xlab=expression('Mass ['*M[Jup]*']'),ylab=expression(T[eff]*' [K]'),xlim=rev(c(55,75)),ylim=c(1550,1750))
plot(unlist(mass),unlist(teff),col='white',xlab=expression('Mass ['*M[Jup]*']'),ylab=expression(T[eff]*' [K]'),xlim=rev(c(1,15)),ylim=c(150,400))
age.select <- c(1,2,4,6,8,10)
ts <- seq(320,280,length.out=length(age.select))
n <- 0
index <- sort(ages,index.return=TRUE,decreasing=TRUE)$ix
for(j in index){
    if(any(ages[j]==age.select)){
        n <- n+1
        ind <- sort(mass[,j],index.return=TRUE)$ix
        lines(mass[ind,j],teff[ind,j],lty=1,col='darkgrey')
                                        #    lines(mass[ind,j],lum[ind,j],lty=1,col='grey')
                                        #    Ntxt <- round(nrow(mass)/2)
        Ntxt <- which.min(abs(teff[,j]-ts[n]))
        text(mass[Ntxt,j],teff[Ntxt,j],pos=3,labels=paste(ages[j],'Gyr'),col='darkgrey')
    }
}
cols <- c('brown','black','green','blue')
points(m,T,col=cols,pch=20)
arrows(m-em,T,m+em,T,length=0.1,angle=90,code=3,col=cols,cex=0.8,pch=20)
arrows(m,T-eTl,m,T+eTu,length=0.1,angle=90,code=3,col=cols,cex=0.8,pch=20)
#legend('topright',legend=star,col=cols,pch=20,bty='n')
legend('topright',legend=star[2],col=cols[2],pch=20,bty='n',cex=size)
dev.off()
