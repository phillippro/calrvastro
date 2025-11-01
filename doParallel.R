library(foreach)
library(doMC)
library(doParallel)

# use the environment variable SLURM_NTASKS_PER_NODE to set the number of cores
#registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
Ncores <- 4
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}

# Bootstrapping iteration example
x <- iris[which(iris[,5] != "setosa"), c(1,5)]
iterations <- 10000 # Number of iterations to run

# Parallel version of code
# Note the '%dopar%' instruction
parallel_time <- system.time({
  r <- foreach(icount(iterations), .combine=cbind) %dopar% {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
  }
})

# Shows the number of Parallel Workers to be used
getDoParWorkers()

# Prints the total compute time.
parallel_time["elapsed"]
