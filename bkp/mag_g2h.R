args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    g1 <- as.numeric(args[1])
    g2 <- as.numeric(args[2])
    br1 <- as.numeric(args[3])
    br2 <- as.numeric(args[4])
}else{
    g1 <- 8.6432
    g2 <- 12.2781
    br1 <- 1.0301
    br2 <- 2.6280
}
dg <- g1-g2

gh <- function(br){
a0 <- 0.0172
a1 <- -0.4565
a2 <- -0.0666
a3 <- 0.0076
sigma <- 0.04
a0+a1*br+a2*br^2+a3*br^3
}

gh1 <- gh(br1)
gh2 <- gh(br2)
dh <- dg-gh1+gh2
f21 <- 2.512^(dh)
cat("Hipparcos brightness ratio of secondary relative to primary:",f21,'\n')