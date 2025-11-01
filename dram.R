###This routine is used to conduct adaptive MCMCs to draw posteriors of models for parameter inference
prepare_parallel_worker <- function(){
    if(!exists('.dram_worker_ready', envir = .GlobalEnv, inherits = FALSE) ||
       !isTRUE(get('.dram_worker_ready', envir = .GlobalEnv))){
        source('mcmc_func.R')
        assign('.dram_worker_ready', TRUE, envir = .GlobalEnv)
    }
}

###Make sure the current (master) session also has the dependencies ready
prepare_parallel_worker()
Pmin0 <- Pmin
Pmax0 <- Pmax
#Npar0 <- Npar
logdP <- (log(Pmax0)-log(Pmin0))/Ncores
mcmc.all <- list()
mcmc.opt <- list()
kep.up <- kep.low <- kep.opt <- c()
Nsig <- 0
basis0 <- basis
bases <- c(rep(basis,Ntr),rep(basis2,Nmax-Ntr))
ind.transit0 <- ind.transit
ll.prim <- per.prim <- c()
ll0 <- NA
data.astrometry0 <- out$astrometry

#####mcmc loop
backend <- foreach::getDoParName()
workers <- foreach::getDoParWorkers()
cat(sprintf("foreach backend: %s (%d workers)\n", backend, workers))
must_have <- c(
  # loop / structure
  "basis0","basis2","bases","Ntr0","ind.transit0","Pmin0","Pmax0",
  "Nmin","Nmax","Ncores","Niter0","Nkeppar","kep.opt",
  # data container and frequently used fields
  "out","data.astrometry0",
  # IC / hot_chain prerequisites
  "rmax","rmin","Ptransit","ePtransit","Tc","eTc","tmin","tmax",
  "RV.all","dP","alpha.dTc","alpha.dP","offset","jitter","rel.name",
  "ind.binary","rvc.type","par.fix","fpar"
)
export_vars <- must_have[ vapply(must_have, exists, logical(1)) ]
if (!length(export_vars)) export_vars <- "prepare_parallel_worker"  # always export at least one

# helper to print what we export (optional)
cat("Exporting to workers:", paste(export_vars, collapse=", "), "\n")

for (np in Nmin:Nmax) {

  # Build the foreach spec with backend-specific options *inside* foreach()
    fe <- foreach(
        ncore = 1:Ncores,
        .errorhandling = "pass",
        .inorder = FALSE,
        .export = c("prepare_parallel_worker", export_vars),
        .options.snow = if (backend %in% c("doParallelSNOW","doSNOW","doParallel")) list(preschedule = FALSE) else NULL,
        .options.multicore = if (backend %in% c("doParallelMC","doMC")) list(preschedule = FALSE) else NULL
    )
  # -------------------- PARALLEL SECTION (workers) --------------------
  mcmc <- fe %dopar% {
    prepare_parallel_worker()
    tryCatch({
      replacePar <- FALSE

      if (np > 0) {
        if (np > Ntr0) {
          Ntr <- 0
          ind.transit <- 0
        } else {
          Ntr <- Ntr0
          ind.transit <- ind.transit0
        }
        Np <- 1
        basis <- if (np > Ntr0) basis2 else basis0

        cat('\n Find signal', np, 'using hot chain!\n')
        bases0 <- bases
        bases  <- basis

        if (Nmin == Nmax) {
          Np    <- Nmin
          bases <- bases0
        }

        # init
#        save(list=ls(all=TRUE),file='test2.Robj')
        sys.source('initial_condition.R', envir = environment())
        RVonly <- FALSE
        if (exists("fpar") && file.exists(fpar)) sys.source('reload_optpar.R', envir = environment())

        par.hot <- startvalue
        if (!is.null(par.fix) && any(par.fix == 'Mstar')) par.hot['Mstar'] <- out$Mstar
        if (any(startvalue < par.min | startvalue > par.max)) {
          cat('The following initial values exceed the prior boundary:',
              names(startvalue)[startvalue < par.min | startvalue > par.max], '\n')
        }
        if (Niter0 >= 1e3) sys.source('hot_chain.R', envir = environment())
        bases <- bases0
      }

      # cold chain
      cat('\n Constrian signal', np, 'using cold chain!\n')
      Ntr <- Ntr0; ind.transit <- ind.transit0
      Np  <- np
      Pmin <- Pmin0; Pmax <- Pmax0

      if (out$Nrv > 0) {
        for (i in out$ins.rv) out[[i]]$RV[,2] <- out[[i]]$data[,2]
      }
      if (out$Nastro > 0) out$astrometry <- data.astrometry0

      sys.source('initial_condition.R', envir = environment())
      if (Np > 1 && Nmin <= 1) {
        startvalue <- as.numeric(c(kep.opt, par.hot))
        startvalue <- assign.names(startvalue, Np = Np, p = ps, q = qs, n = ns, bases = bases)
      } else {
        startvalue <- par.hot
      }

      tmp <- run.metropolis.MCMC(
        startvalue, cov.start,
        iterations = min(1e7, Niter0 * max(Nsig, 1)),
        out1 = out, tem = 1, bases = bases
      )

      # burn-in
      outmat <- if (Niter0 > 1e5) {
        tmp$out[-(1:floor(nrow(tmp$out) / 2)), , drop = FALSE]
      } else {
        tmp$out
      }

      if (!is.matrix(outmat) || nrow(outmat) == 0) stop("empty chain returned")
      outmat  # success return

    }, error = function(e) {
      # Tag the error so we can filter on master
      structure(list(.error = conditionMessage(e)))
    })
  }
  # -------------------- END PARALLEL SECTION --------------------

  # Split successes vs errors
  is_err <- vapply(mcmc, function(x) is.list(x) && !is.null(x$.error), logical(1))
  if (any(is_err)) {
    cat("Worker errors:\n",
        paste0(which(is_err), ": ",
               vapply(mcmc[is_err], `[[`, "", ".error"), collapse = "\n"),
        "\n")
  }
  mcmc <- mcmc[!is_err]

  # Drop non-matrix / empty results
  is_good <- vapply(mcmc, function(mm) is.matrix(mm) && nrow(mm) > 0 && ncol(mm) > 1, logical(1))
  if (any(!is_good)) {
    bad_idx <- which(!is_good)
    cat("Dropping", length(bad_idx), "chains (empty/non-matrix):",
        paste(bad_idx, collapse = ", "), "\n")
    mcmc <- mcmc[is_good]
  }

  # HARD GUARD: no chains left
  if (!length(mcmc)) {
    stop(sprintf("All parallel chains failed for np=%s. See worker errors above.", np))
  }

  # Save mcmc results
  if (!save.memory) mcmc.all[[paste0('sig', np)]] <- mcmc

  # downstream calculations (safe indexing)
  llmax <- vapply(mcmc, function(mm) max(mm[, ncol(mm)]), numeric(1))
  cat('Maximum likelihood for all chains:', llmax, '\n')

  best_idx <- which.max(llmax)
  mcopt <- mcmc[[best_idx]]
  mcmc.opt[[paste0('sig', np)]] <- mcopt

  if (np == 1) {
    per.prim <- vapply(mcmc, function(mm) mm[ which.max(mm[, ncol(mm) - 1]), 1], numeric(1))
    ll.prim  <- llmax
  }

  Npar <- ncol(mcopt) - 2L
  if (Npar < 1) stop("Unexpected chain shape (need at least 3 columns: params, logpost, loglike).")

  par.opt <- mcopt[ which.max(mcopt[, 'loglike']), 1:Npar, drop = TRUE ]
  if (np > 0) cat('optimal period:', extract.par(par.opt, bases = bases)$P, 'days\n')

  if (np > 0) {
    kep.opt <- par.opt[1:(np * Nkeppar)]
    res <- loglikelihood(par.opt, predict = TRUE)$res
    if (out$Nrv > 0) {
      for (i in out$ins.rv) out[[i]]$RV[, 2] <- as.numeric(res[[i]]$res1)
    }
  }

  ll <- max(mcopt[, ncol(mcopt)])
  if (np > 0) {
    if (is.na(ll0)) ll0 <- ll - 1e2
    lnbf3 <- ll - ll0 - 1.5 * log(length(unlist(res)))
    cat('lnbf3=', lnbf3, '\n')
    Nsig <- Nsig + 1
  }
  ll0 <- ll
}
if(Nmin==Nmax) Nsig <- Nmin
###change the data back into raw data
if(out$Nrv>0){
    for(i in out$ins.rv){
        out[[i]]$RV[,2] <- out[[i]]$data[,2]
    }
}
if(any(names(out)=='astrometry')){
    out$astrometry <- data.astrometry0
}
mcopt <- mcmc.opt[[paste0('sig',Nsig)]]
Npar <- ncol(mcopt)-2
par.opt <- mcopt[which.max(mcopt[,'logpost']),1:Npar]
if(Nsig>0){
    Popt <- extract.par(par.opt,out1=out,bases=bases)$P
}else{
    Popt <- NA
}
if(!save.memory){
    out[['mcmc.all']] <- mcmc.all
}
out[['mcmc.opt']] <- mcmc.opt
out$tauT <- RV.kepler(par.opt,out1=out)$tauT
