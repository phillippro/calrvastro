# ---- Setup ---------------------------------------------------------------
library(doParallel)
library(doRNG)        # reproducible foreach
# library(iterators)  # only if you need advanced iterators

nii <- FALSE
Pmin0 <- Pmin; Pmax0 <- Pmax
logdP <- (log(Pmax0) - log(Pmin0)) / Ncores

mcmc.all <- list()
mcmc.opt  <- list()
kep.up <- kep.low <- kep.opt <- c()
Nsig <- 0
basis0 <- basis
bases  <- c(rep(basis, Ntr), rep(basis2, Nmax - Ntr))
ind.transit0 <- ind.transit
ll.prim <- per.prim <- c()
ll0 <- NA
data.astrometry0 <- out$astrometry

# ---- Cluster: start once, preload once ----------------------------------
cl <- makeCluster(Ncores)
registerDoParallel(cl)

# Legacy alias if some code still reads Niter
Niter <- Niter0

# Export ONLY what workers need for hot_chain and likelihood/MCMC
clusterExport(cl, c(
  # hot_chain.R needs:
  "Niter0","Ncores","cov.start","RVonly",
  # state that hot_chain reads/writes:
  "per.prim","ll.prim",'rmax','rmin','Np','bases','Pmin','Pmax','Ntr','RV.all','dP','tmax','tmin','par.global','trv.all',
  # sometimes referenced by helpers:
  "out"
), envir = environment())

# Preload packages + core funcs once (but DO NOT source hot_chain.R here)
clusterEvalQ(cl, {
  RNGkind("L'Ecuyer-CMRG")
  # Put any library() calls here if needed by your sourced files.
  source('mcmc_func.R',        local = TRUE)  # defines run.metropolis.MCMC, loglikelihood, etc.
  source('initial_condition.R', local = TRUE) # if it defines helpers (safe to load once)
  if (file.exists('reload_optpar.R')) source('reload_optpar.R', local = TRUE)
  # DO NOT source('hot_chain.R') here — it’s a script and needs per-task variables.
  if (!exists("Niter") && exists("Niter0")) Niter <- Niter0
  NULL
})

clusterSetRNGStream(cl, 123)

for (np in Nmin:Nmax) {

  if (!nii) {
    # Prepare immutable inputs for workers: create a compact list to avoid exporting the whole workspace.
    # Anything big (like 'out') should be trimmed to what's actually used inside likelihood.
    out_worker <- out              # <- best is to make a lightweight view: e.g., drop unused fields before sending
    out_worker$astrometry <- data.astrometry0

    # Reproducible parallel foreach over chains
    mcmc <- foreach(
      ncore = 1:Ncores,
      .packages = character(),     # add packages used INSIDE worker if needed
      .export   = c("run.metropolis.MCMC", "assign.names", "extract.par",
                    "loglikelihood", "RV.kepler"),  # functions referenced directly
      .combine  = 'list',
      .multicombine = TRUE,
      .maxcombine = Ncores,
      .errorhandling = 'pass'
    ) %dorng% {

      # --------------------- Worker body (pure function style) ---------------------
      # Build all per-chain state locally to avoid global mutation.
      # NOTE: Don't use <<- or modify outer variables!

      # Select bases / Np / init depending on np
      Ntr <- Ntr0; ind.transit <- ind.transit0
      bases0 <- c(rep(basis0, Ntr0), rep(basis2, max(0, Nmax - Ntr0)))
      bases  <- if (np > Ntr0) basis2 else basis0
      Np     <- if (Nmin == Nmax) Nmin else np

      # Special cases
      if (np > 0) {
        # Initial conditions
        source('initial_condition.R')# should expose a function; if not, wrap its logic in a function once.
        # Here we assume it defines 'startvalue', 'par.min', 'par.max', etc. in the worker env.
        # If 'fpar' is involved, guard its existence:
        if (exists("fpar") && file.exists(fpar)) {
          source('reload_optpar.R')
        }

        # Prepare hot chain start
        par.hot <- startvalue
        if (!is.null(par.fix) && any(par.fix == 'Mstar')) {
          par.hot['Mstar'] <- out_worker$Mstar
        }

        # Optional hot chain
        if (Niter0 >= 1e3) {
            source('hot_chain.R')
          # Assumes hot_chain.R defines a function that updates par.hot in-scope
          # e.g., par.hot <- run_hot_chain(par.hot, ...)
        }
      }

      # Cold chain setup (independent of hot chain state)
      Ntr <- Ntr0; ind.transit <- ind.transit0; Np <- np
      Pmin <- Pmin0; Pmax <- Pmax0

      # Reset data inside the local 'out_worker' copy
      if (out_worker$Nrv > 0) {
        for (ins in out_worker$ins.rv) {
          out_worker[[ins]]$RV[, 2] <- out_worker[[ins]]$data[, 2]
        }
      }
      if (out_worker$Nastro > 0) out_worker$astrometry <- data.astrometry0

      # Recompute initial condition for cold chain
      # (Again, best if this is a function you call with explicit args)
      # initial_condition()

      if (Np > 1 && Nmin <= 1) {
        startvalue <- as.numeric(c(kep.opt, par.hot))
        startvalue <- assign.names(startvalue, Np = Np, p = ps, q = qs, n = ns, bases = bases)
      } else {
        startvalue <- par.hot
      }

      # Run MCMC
      iter_target <- min(1e7, Niter0 * max(Nsig, 1L))
      res <- run.metropolis.MCMC(
        startvalue,
        cov.start,
        iterations = iter_target,
        out1 = out_worker,
        tem = 1,
        bases = bases
      )

      # Burn-in truncation on worker to return less data
      outmat <- res$out
      if (Niter0 > 1e5) {
        outmat <- outmat[-(1:floor(nrow(outmat) / 2)), , drop = FALSE]
      }

      # Return a compact result: best row + (optionally) thin chain
      # To cut memory, keep full chain only if you need it downstream.
      list(
        best_row = outmat[which.max(outmat[, ncol(outmat)]), , drop = FALSE],
        chain    = if (get("save.memory", inherits = FALSE, ifnotfound = TRUE)) NULL else outmat
      )
      # ---------------------------------------------------------------------------
    }

  } else {
    # Non-parallel path
    # initial_condition()
    mcmc <- list()
  }

  # --------------------- Master: combine results ---------------------
  # Filter errors
  bad <- which(vapply(mcmc, function(z) !is.list(z) || is.null(z$best_row), logical(1)))
  if (length(bad)) mcmc <- mcmc[-bad]

  # Save chains if requested
  if (!save.memory) {
    mcmc.all[[paste0('sig', np)]] <- lapply(mcmc, `[[`, "chain")
  }

  # Find best by loglike (assumes last column is loglike)
  llmax  <- vapply(mcmc, function(z) z$best_row[1, ncol(z$best_row)], numeric(1))
  mcopt  <- mcmc[[which.max(llmax)]]$chain
  mcmc.opt[[paste0('sig', np)]] <- mcopt

  if (np == 1) {
    per.prim <- vapply(mcmc, function(z) z$best_row[1, ncol(z$best_row) - 1], numeric(1))
    ll.prim  <- llmax
  }

  # Subtract the best-fit for next signal
  Npar   <- ncol(mcopt) - 2L
  par.opt <- mcopt[which.max(mcopt[, 'loglike']), 1:Npar]
  if (np > 0) {
    kep.opt <- par.opt[1:(np * Nkeppar)]
    res <- loglikelihood(par.opt, predict = TRUE)$res
    if (out$Nrv > 0) {
      for (ins in out$ins.rv) out[[ins]]$RV[, 2] <- as.numeric(res[[ins]]$res1)
    }
  }

  # BF criterion
  ll <- max(mcopt[, ncol(mcopt)])
  if (np > 0) {
    if (is.na(ll0)) ll0 <- ll - 1e2
    lnbf3 <- ll - ll0 - 1.5 * log(length(unlist(res)))
    # decide to break or continue based on lnbf3 if desired
    Nsig <- Nsig + 1
  }
  ll0 <- ll
}

# Reset data
if (out$Nrv > 0) for (ins in out$ins.rv) out[[ins]]$RV[,2] <- out[[ins]]$data[,2]
if ('astrometry' %in% names(out)) out$astrometry <- data.astrometry0

mcopt <- mcmc.opt[[paste0('sig', Nsig)]]
Npar  <- ncol(mcopt) - 2L
par.opt <- mcopt[which.max(mcopt[, 'logpost']), 1:Npar]
Popt <- if (Nsig > 0) extract.par(par.opt, out1 = out, bases = bases)$P else NA

if (!save.memory) out[['mcmc.all']] <- mcmc.all
out[['mcmc.opt']] <- mcmc.opt
out$tauT <- RV.kepler(par.opt, out1 = out)$tauT

stopCluster(cl)
