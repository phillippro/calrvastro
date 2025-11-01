m1 <- 11.7
ls2au <- 0.0020039888
a2sini <- 856.26465*ls2au
Pd <- 269.436227
Py <- Pyr <- 269.436227/365.25
K1 <- 8322.365512
#    a1 <- (K/sinI)^2/(4*pi^2)*(1-e^2)*P^(2/3)
K2 <- 2*pi*a2sini/Pyr/sqrt(1-e^2)*4.74047e3
q <- K1/K2
m2 <- m1*q

mtot <- m1*(1+q)
a2 <- (Pyr^2*mtot)^(1/3)/(1+q)
i <- asin(a2sini/a2)
