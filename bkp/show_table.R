show.table <- function(table,alpha=0.3,b=0){
    fonts <- rep(1,length(table))
    plot.new()
    cn <- c('',sapply(1:length(table),function(i) table[[i]]['name']),'')
cat('cn=',cn,'\n')
    rn1 <- names(table[[1]][-1])
    rn2 <- names(table[[length(table)]][-1])
#    legend('topleft',inset=c(0,b-0.1),legend=rn,bty='n',xpd=NA,horiz=TRUE)
#    legend('topleft',inset=c(0,b),legend=cn,bty='n',xpd=NA)
    legend('topleft',inset=c(0,b),legend=c('',rn1),bty='n',xpd=NA,text.font=2)
    for(j in 1:length(table)){
        legend('topleft',inset=c(alpha*j,b),legend=table[[j]],bty='n',xpd=NA,text.font=c(2,rep(fonts[j],length(table[[j]])-1)))
    }
    legend('topleft',inset=c((length(table)+1)*alpha,b),legend=c('',rn2),bty='n',xpd=NA,text.font=2)
}

###diagnostics
Nmodel <- Noverlap <- quality <- rep(NA,Nsig)
for(i in 1:Nsig){
        n1 <- sum(as.numeric(sapply(1:ncol(Psig),function(k) any(abs(Psig[,k]-Pkep[i])/Pkep[i]<0.1))))
        Nmodel[i] <- n1
        n2 <- sum(as.numeric(sapply(1:nrow(Pact),function(k) any(abs(Pact[k,]-Pkep[i])/Pkep[i]<0.1))))
        Noverlap[i] <- n2
        score <- 0
###activity dependency
        if(n2==0){
            score <- score+1
        }
###noise model dependency
        if(n1==length(noise.types)){
            score <- score+1
        }
###chunck dependency
        if(gammas[i]<1 | is.na(gammas[i]) | (gamma>=1 & Pkep[i]>max(dTs)/2)){
            score <- score+1
        }
###statistical significant
        if(logBF5[i]>=5){
            score <- score+1
        }
        if(score==4) quality[i] <- 'A'
        if(score==3) quality[i] <- 'B'
        if(score==2) quality[i] <- 'C'
        if(score==1) quality[i] <- 'D'
        if(score==0) quality[i] <- 'E'
}

#####show derived parameters of each fit
par.out <- list()
for(j in 1:(Nsig+1)){
    if(j<=Nsig){
        par.show <- par.kep[Nkeppar*(j-1)+(1:Nkeppar)]
        if(period.par=='logP'){
            par.show[1] <- exp(par.show[1])
        }
        tmp <- c(paste(j,'signal'),round(logBF3[j],3),round(logBF5[j],3),round(par.show,3),Nmodel[j],Noverlap[j],round(gammas[j],3),quality[j])
        if(basis=='natural'){
        names(tmp) <- c('name','logBF3','logBF5','Period[day]','K[m/s]','e','omega[rad]','Mo[rad]','Nmodel','Noverlap','gamma','quality')
}else{
        names(tmp) <- c('name','logBF3','logBF5','Period[day]','ln(K[m/s])','sqrt(e)*sin(omega)','sqrt(e)*cos(omega)','Tc[JD]','Nmodel','Noverlap','gamma','quality')
}
        par.out[[j]] <- tmp
    }else{
        par.out[[j]] <- c('noise',round(par.kep[(length(par.kep)-Nkeppar*Nsig+1):length(par.kep)],3))
    }
}
#fonts <- rep(1,ncol(par.out))
#inds <- ncol(par.shows)+Ntype*(0:5)+ind.show+1
#fonts[inds] <- 2
show.table(par.out,alpha=0.2)
#fonts=fonts
