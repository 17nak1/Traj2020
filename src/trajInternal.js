

let trajectory_internal = function (object, params, times, t0, as_data_frame = false, _getnativesymbolinfo = true, verbose) {

  times = [1,6,3,5];

  let ep = 'in \"trajectory\": ';

  verbose = verbose ? true : false;
  // as.data.frame <- as.logical(as.data.frame)

  if (times === undefined || !Array.isArray(times))
    times = object.times;
  // else
  //   times <- as.numeric(times)

  if (!Array.isArray(times) || times.length == 0)
    throw new Error(ep + '\"times\" is empty, there is no work to do');

  let isAscending = times.slice(1).map((e,i) => e > times[i]).every(x => x);
  if (!isAscending)
    throw new Error(ep + '\"times\" must be an increasing sequence of times');

  if (t0 === undefined)
    t0 = object.t0;
  // else
  //   t0 <- as.numeric(t0)

  if (t0 > times[0])  
    throw new Error(ep + '\"times\" the zero-time \"t0\"  must occur no later than the first observation');

  ntimes = times.length;

  if (params === undefined) params = object.params;
  // if (is.list(params)) params <- unlist(params)
  if (Object.keys(params).length === 0) 
    throw new Error(ep + '\"params\" must be supplied');
  // storage.mode(params) <- "double"

  // params <- as.matrix(params)
  // nrep <- ncol(params)
  // paramnames <- rownames(params)
  // if (is.null(paramnames))
  //   stop(ep,sQuote("params")," must have rownames",call.=FALSE)

  // x0 <- init.state(object,params=params,t0=t0)
  // nvar <- nrow(x0)
  // statenames <- rownames(x0)
  // dim(x0) <- c(nvar,nrep,1)
  // dimnames(x0) <- list(statenames,NULL,NULL)

  // type <- object@skeleton.type          # map or vectorfield?

  // pompLoad(object,verbose=verbose)
  // on.exit(pompUnload(object,verbose=verbose))

  // if (type=="map") {

  //   x <- tryCatch(
  //     .Call(iterate_map,object,times,t0,x0,params,.getnativesymbolinfo),
  //     error = function (e) {
  //       stop(ep,"in map iterator: ",
  //         conditionMessage(e),call.=FALSE)
  //     }
  //   )
  //   .getnativesymbolinfo <- FALSE

  // } else if (type=="vectorfield") {

  //   znames <- object@zeronames
  //   if (length(znames)>0) x0[znames,,] <- 0

  //   .Call(pomp_desolve_setup,object,x0,params,.getnativesymbolinfo)
  //   .getnativesymbolinfo <- FALSE

  //   X <- tryCatch(
  //     ode(
  //       y=x0,
  //       times=c(t0,times),
  //       method="lsoda",
  //       func="pomp_vf_eval",
  //       dllname="pomp",
  //       initfunc=NULL,
  //       parms=NULL,
  //       ...
  //     ),
  //     error = function (e) {
  //       stop(ep,"error in ODE integrator: ",conditionMessage(e),call.=FALSE)
  //     }
  //   )

  //   .Call(pomp_desolve_takedown)

  //   if (attr(X,"istate")[1L]!=2)
  //     warning(ep,"abnormal exit from ODE integrator, istate = ",attr(X,'istate')[1L],
  //       call.=FALSE) # nocov

  //   if (verbose) {
  //     deSolve::diagnostics(X)
  //   }

  //   x <- array(data=t(X[-1L,-1L]),dim=c(nvar,nrep,ntimes),
  //     dimnames=list(statenames,NULL,NULL))

  //   for (z in znames)
  //     for (r in seq_len(ncol(x)))
  //       x[z,r,-1] <- diff(x[z,r,])

  // } else {

  //   stop(ep,"deterministic skeleton has not been properly specified",call.=FALSE)

  // }

  // dimnames(x) <- setNames(dimnames(x),c("variable","rep","time"))

  // if (as.data.frame) {
  //   x <- lapply(
  //     seq_len(ncol(x)),
  //     function (k) {
  //       nm <- rownames(x)
  //       y <- x[,k,,drop=FALSE]
  //       dim(y) <- dim(y)[c(1L,3L)]
  //       y <- as.data.frame(t(y))
  //       names(y) <- nm
  //       y$time <- times
  //       y$traj <- as.integer(k)
  //       y
  //     }
  //   )
  //   x <- do.call(rbind,x)
  //   x$traj <- factor(x$traj)
  // }

  // x
}

trajectory_internal();

exports.trajectory_internal = trajectory_internal;