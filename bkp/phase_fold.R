for(i3 in 1:Np){
t3c <- datac[[k1]][,i3]
RV3c <- datac[[k1]][,Np+i3]
t4c <- datab[[k1]][,i3]
RV4c <- datab[[k1]][,Np+i3]
eRV4c <- datab[[k1]][,2*Np+i3]
t5c <- datas[[k1]][,i3]
RV5c <- datas[[k1]][,Np+i3]
pchc <- c(1,2,3)
if(data.combine){
    ylim <- c(rmin,rmax)
#if(plot.single){
#    pdf.file <- paste0('P',format(Popt[i3],digit=3),'d_phase_fold.pdf')
#    cat('pdf: ',pdf.file,'\n')
#}
    for(k2 in 1:Nw){
        t3c <- datac[[k2]][,i3]
        RV3c <- datac[[k2]][,Np+i3]
        t4c <- datab[[k2]][,i3]
        RV4c <- datab[[k2]][,Np+i3]
        eRV4c <- datab[[k2]][,2*Np+i3]
        t5c <- datas[[k2]][,i3]
        RV5c <- datas[[k2]][,Np+i3]
        if(k2==1){
            rspan <- max(RV2.all)-min(RV2.all)
            plot(t3c,RV3c,xlab=paste0('orbital phase [days]'),ylab='RV(data-model) [m/s]',pch=pchc[k2],cex=0.2,col=cols[k2],ylim=c(min(RV5c)-rspan/8,max(RV5c)+rspan/8))#col='black',
        }else{
            points(t3c,RV3c,cex=0.2,col=cols[k2])
        }
        points(t4c,RV4c,col=cols[k2],cex=1.2,pch=20)
        arrows(t4c,RV4c-eRV4c,t4c,RV4c+eRV4c,length=0.05,angle=90,code=3,col=cols[k2],cex=1.2,pch=20)
        points(t5c,RV5c,col="red",pch=20,cex=0.2)
#        legend('topright',legend=paste0(i3,'planet'),bty='n')
        legend('topright',legend=paste0('P=',format(Popt[i3],digit=3),'d'),bty='n')
    }
#    if(plot.single){
#        dev.off()
#    }
}else{
###k1 is the number of data set
#    ylim <- range(RV4c[!is.na(RV4c)]-eRV4c[!is.na(eRV4c)],RV4c[!is.na(RV4c)]+eRV4c[!is.na(eRV4c)],RV5c)
    ylim <- range(RV3c,RV5c)
    plot(t3c,RV3c,xlab=paste0('orbital phase[days]'),ylab='RV(data-model)[m/s]',pch=pchc[k1],ylim=ylim,main='',cex=0.2,col='grey')#col='black',
    points(t4c,RV4c,cex=1.2,pch=20,col=cols[k1])
    arrows(t4c,RV4c-eRV4c,t4c,RV4c+eRV4c,length=0.05,angle=90,code=3,col=cols[k1],cex=1.2,pch=20)
    points(t5c,RV5c,col="red",pch=20,cex=0.2)
    legend('topright',legend=paste0(i3,'planet'),bty='n')
}
}
#if(plot.single){
#    dev.off()
#}
