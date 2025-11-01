###ref. Feng+24
aq <- function(a,q,alpha,type='bright'){
    if(type=='bright'){
        eta <- (1-q^(alpha-1))*(1+q^alpha)^(-1)
        cat('eta=',eta,'\n')
       a*q/(1+q)*eta
    }else{
       a*q/(1+q)
    }
}
