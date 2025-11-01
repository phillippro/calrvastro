Ps <- seq(1,1000,by=0.1)
a <- (Ps/365.25)^(2/3)*1e3#mas
theta <- (1000%%Ps)/Ps*(2*pi)#rad
rho <- 2*a*sin(theta/2)#maximum separation
pmp <- rho/1000#assume a baseline of 1000 days
fout <- 'planet_PM.pdf'
cat(fout,'\n')
pdf(fout,4,4)
plot(Ps,pmp,type='l',xlab='Period [day]',ylab='Maximum planet-induced PM [mas/yr]',main='assume 1000 day baseline')
dev.off()


