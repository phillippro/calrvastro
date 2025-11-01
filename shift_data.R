ii <- sort(tmins,index.return=TRUE)$ix
v1 <- t1 <- c()
drvs <- c()
for(k in ii){
      i <- ins[k]
      if(k==ii[1]){
          drv <- 0
      }else{
#	ts <- trv.all[-out[[i]]$index]#other times
#        rv <- RV.all[-out[[i]]$index]#other times
          t <- out[[i]]$RV[,1]
          dts <- sapply(t,function(tim) min(abs(tim-t1)))
          ind <- which.min(dts)
          index <- which.min(abs(t[ind]-t1))
          drv <- -out[[i]]$RV[ind,2]+v1[index]#data shift amount
          drvs <- c(drvs,drv)#model offset
      }
      tt <- out[[i]]$RV[,1]
      vv <- out[[i]]$RV[,2]
      t1 <- c(t1,tt)
      v2 <- vv+drv
      v1 <- c(v1,v2)
      out[[i]]$data[,2] <- out[[i]]$RV[,2] <- v2
}

