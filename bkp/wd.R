##assume double suns
cs <- rainbow(100)
plxs <- seq(10,1000,by=10)#mas
pdf('test.pdf',16,16)
par(mfrow=c(4,4))
dt <- 0.5#GDR3 and GDR2 time gap
for(k in 1:length(plxs)){
    plx <- plxs[k]
    a <- 10^seq(0,4,by=0.01)
    v <- sqrt(4*pi^2*2/a)#au/yr
    pm <- v*plx#mas/yr
    if(k==1){
        plot(a,pm,xlab='a [au]',ylab='pm [mas/yr]',type='l',col=cs[1],log='x')
    }else{
        lines(a,pm,col=cs[k])
    }
    if(k%%10==0) text(max(a),min(pm),labels=paste(plx,'mas'),xpd=NA,pos=4)
}
plxs <- c(1,10,100,1000)
for(plx in plxs){
    a <- 10^seq(2,5,by=0.1)
    v <- sqrt(4*pi^2*2/a)#au/yr
    pm <- v*plx#mas/yr
    plot(a,v*4.74047,xlab='a [au]',ylab='v [km/s]',type='l',col=cs[1],log='x',main=paste(plx,'mas'))
    plot(a,pm*dt,xlab='a [au]',ylab='pm*dt [mas]',type='l',col=cs[1],log='x',main=paste(plx,'mas'))
    g <- 4*pi^2/a^2#au/yr^2
    ang <- 0.5*g*dt^2*plx
    plot(a,ang,xlab='a [au]',ylab='0.5*g*dt^2 [mas]',type='l',col=cs[1],log='x',main=paste(plx,'mas'))
    plot(a,g*dt,xlab='a [au]',ylab='g*dt [mas/yr]',type='l',col=cs[1],log='x',main=paste(plx,'mas'))
}
dev.off()
