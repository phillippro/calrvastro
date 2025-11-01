args <- as.numeric(commandArgs(trailingOnly=TRUE))
if(length(args)>0){
    g <- args[1]
    eg <- args[2]
    bp <- args[3]
    rp <- args[4]
}else{
    g <- 8.6414
    eg <- 0.0010
    bp <- 10.0157
    rp <- 7.3236
}
g2h <- function(g,eg,bp,rp){
    hs <- ehs <- c()
    for(j in 1:3){
        if(j==1){
            c0 <- 0.0172
            c1 <- -0.4565
            c2 <- -0.0666
            c3 <- 0.0076
            egh <- 0.04
            gh <- c0+c1*(bp-rp)+c2*(bp-rp)^2+c3*(bp-rp)^3
            h <- g-gh
            eh <- sqrt(eg^2+egh^2)
            hs <- c(hs,h)
            ehs <- c(ehs,eh)
        }
        if(j==2){
            c0 <- 0.0023
            c1 <- 0.1253
            c2 <- -0.0279
            c3 <- -0.0029
            eh <-0.05
            hbp <- c0+c1*(bp-rp)+c2*(bp-rp)^2+c3*(bp-rp)^3
            h <- hbp+bp
            hs <- c(hs,h)
            ehs <- c(ehs,eh)
        }
        if(j==3){
            c0 <- 0.0023
            c1 <- 1.1253
            c2 <- -0.0279
            c3 <- -0.0029
            eh <- 0.05
            hrp <- c0+c1*(bp-rp)+c2*(bp-rp)^2+c3*(bp-rp)^3
            h <- hrp+rp
            hs <- c(hs,h)
            ehs <- c(ehs,eh)
        }
    }
    cbind(hs,ehs)
}
hs <- g2h(g,eg,bp,rp)
cat('Hp mag=',hs[,1],'mag\n')
cat('Error of Hp mag=',hs[,2],'mag\n')

###mag to flux
###10^(dm/2.5) or 2.512^dm
