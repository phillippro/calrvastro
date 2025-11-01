##light second to au
##parameters for CPD-632495 from the following papers:
##The geometric distance and binary orbit of PSR B1259–63, Miller-Jones+2018
##The kinematics and orbital dynamics of the PSR B1259−63/LS 2883 system from 23 yr of pulsar timing, Shannon+2014
ls2au <- 0.0020039888
P <- 1236.724526
eP <- 6e-6
Pyr <- P/365.25
ePyr <- eP/365.25
Pyrs <- rnorm(1e3,Pyr,ePyr)
Tp  <- 2400000.5+53071.2447290
eTp <- 7e-7
omega2 <- 2.420161
eomega2 <- 1.9e-6
omega <- (omega2+pi)%%(2*pi)
eomega <- 1.9e-6
Omega <- 189/180*pi
eOmega <- 2/180*pi
e <- 0.86987970
ee <- 6e-8
es <- rnorm(1e3,e,ee)
asini <- 1296.27448*ls2au
easini <- 0.00014*ls2au
###the orbital parameters for the visible star
i <- 2.687807
ei <- 0.05235988
Is <- rnorm(1e3,i,ei)
asinis <- rnorm(1e3,asini,ei)
a <- asini/sin(i)
As <- asinis/sin(Is)
plx <- 0.38
eplx <- 0.05
per1 <- log(P)
dper1 <- (log(P+eP)-log(P-eP))/2
tmin <- 2448348.75
Mo <- 2*pi*((tmin-Tp)%%P)/P
Ps <- rnorm(1e3,P,eP)
eMo <- sd(2*pi*((tmin-rnorm(1e3,Tp,eTp))%%Ps)/Ps)
mtot <- a^3/Pyr^2
mtots <- As^3/Pyrs^2
M2 <- 1.4
M1 <- mtots-M2
k <- sqrt(4*pi^2*M2^2/mtot/a/(1-e^2))*sin(i)*4.74047e3#m/s
ks <- sqrt(4*pi^2*M2^2/mtots/As/(1-es^2))*sin(Is)*4.74047e3



