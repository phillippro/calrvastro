#************************************************************************************
#This code calculates habitable zone 'fluxes' using the expression given in the 
# Kopparapu et al.(2014) paper. The corresponding output file is 'HZs.dat'. 
#! It also generates a file 'HZ_coefficients.dat' that gives the coefficients for 
#! the analytical expression.#
#
#!
#!
#! Ravi kumar Kopparapu April 19 2014
#!************************************************************************************##
#
#implicit none
#real *8 seff(6),seffsun(6),teff,a(6),b(6),c(6),d(6),tstar
#integer i

#!*********#***************************************************************************
#! Output files.#

#open(9,file='HZs.dat')
#open(10,file='HZ_coefficients.dat')

#!************************************************************************************
#! Coeffcients to be used in the analytical expression to calculate habitable zone flux 
#! boundaries


seffsun <- c(1.776,1.107,0.356,0.320,1.188,0.99)

a <- c(2.136e-4,1.332e-4,6.171e-5,5.547e-5,1.433e-4,1.209e-4)
b <- c(2.533e-8,1.580e-8,1.698e-9,1.526e-9,1.707e-8,1.404e-8)
c <- c(-1.332e-11,-8.308e-12,-3.198e-12,-2.874e-12,-8.968e-12,-7.418e-12)
d <- c(-3.097e-15,-1.931e-15,-5.575e-16,-5.011e-16,-2.084e-15,-1.713e-15)

#!************************************************************************************
#! Calculating HZ fluxes for stars with 2600 K < T_eff < 7200 K. The output file is
#! 'HZ_fluxes.dat'
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    teff <- as.numeric(args[1])##effective temperature
    l <- as.numeric(args[2])#luminosity
    mass <- as.numeric(args[3])
}else{
    teff <- 3131#2600.0d0
#    teff <- 5777
    l <- 0.00298#luminosity
#    l <- 1
#    mass <- 7#stellar mass in solar unit
    mass <- 1
}
ms <- c(1,5,0.1)
#ind <- which.min(abs(ms-mass))
tstar <- teff - 5780.0
seff <- rep(NA,length(seffsun))
for(i in 1:6){
     seff[i] = seffsun[i] + a[i]*tstar + b[i]*tstar^2 + c[i]*tstar^3 + d[i]*tstar^4
}
d <- (l/seff)^0.5
#cat('Real planetary mass:',mass,'Me\n')
#cat('approximate planetary mass:',ms[ind],'\n')
ind <- 2
cat('inner HZ (Recent Venus limit)=',d[1],'au; P=',format(365.24*sqrt(d[1]^3/mass),digit=2),'d; stellar flux:',seff[1],'Earth flux\n')
cat('5 Earth mass inner HZ (Runaway Greenhouse limit)=',d[3+ind],'au; P=',format(365.24*sqrt(d[3+ind]^3/mass),digit=2),'d; stellar flux:',seff[3+ind],'\n')
cat('1 Earth mass inner HZ (Runaway Greenhouse limit)=',d[2],'au; P=',format(365.24*sqrt(d[2]^3/mass),digit=2),'d; stellar flux:',seff[2],'\n')
ind <- 3
cat('0.1 Earth mass inner HZ (Runaway Greenhouse limit)=',d[3+ind],'au; P=',format(365.24*sqrt(d[3+ind]^3/mass),digit=2),'d; stellar flux:',seff[3+ind],'\n\n')
cat('outer HZ (Maximum Greenhouse limit)=',d[3],'au; P=',format(365.24*sqrt(d[3]^3/mass),digit=2),'d; stellar flux:',seff[3],'\n')
cat('outer HZ (Early Mars limit)=',d[4],'au; P=',format(365.24*sqrt(d[4]^3/mass),digit=2),'d; stellar flux:',seff[4],'\n')





