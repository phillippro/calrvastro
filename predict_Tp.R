source('timing_function.R')
source('mcmc_func.R')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    star <- args[1]
    M0 <- as.numeric(args[2])
    P <- as.numeric(args[3])
}else{
   star <- 'HD17382'
   M0 <- 2.583087
   P <- 5523.9
}
fs <- list.files(path=paste0('../data/combined/',star),pattern='dat$|vels$|rv$',full.name=TRUE)
fs <- fs[!grepl('photo',fs)]
tmin  <- c()
out <- c()
for(f in fs){
    cat(f,'\n')
    tab <- read.table(f,header=TRUE)
    out <- rbind(out,as.matrix(tab[,1:3]))
	tmin <- c(tmin,min(tab[,1]))
}
T0 <- min(tmin)
jd0 <- sum(time_Cal2JD(cbind(2021,9,1)))
cat(jd0,'\n')
Tp <- M02Tp(M0,T0,P)+(1:200)*P
cat('tail(Tp)=',tail(Tp),'\n')
ind <- which(Tp>jd0 & Tp<jd0+5*365.25)
Tp <- Tp[ind]
cat('head(Tp)=',head(Tp),'\n')
cals <- time_Jd2cal(cbind(Tp,0))
for(j in 1:nrow(cals)){
    cat(cals[j,],'\n')
}
fout <- paste0('results/',star,'/',star,'_Tp.pdf')
cat(fout,'\n')
pdf(fout,6,6)
ind <- which(Tp>min(out[,1]) & Tp<max(out[,1]))
ccs <- c()
for(j in ind){
yr <- cals[j,1]
mon <- cals[j,2]
day <- cals[j,3]+cals[j,4]
ccs <- c(ccs,paste0(yr,'-',mon,'-',day))
}
plot(out[,1],out[,2],xlab='BJD',ylab='RV[m/s]',main=ccs)
abline(v=Tp,col='green')
abline(v=Tp+30,col='red')
abline(v=Tp-30,col='red')
abline(v=Tp+365,col='blue')
abline(v=Tp-365,col='blue')
dev.off()

