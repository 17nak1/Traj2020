const { initState } = require("./initState.js");
const { iterateMap } = require("./trajectory.js");

exports.trajectory = function (object, params, times, t0, asDataFrame, args){

  let ep = 'in \"trajectory\": ';

  if (times === undefined || !Array.isArray(times))
    times = object.times;

  if (!Array.isArray(times) || times.length == 0)
    throw new Error(ep + '\"times\" is empty, there is no work to do');

  let isAscending = times.slice(1).map((e,i) => e > times[i]).every(x => x);
  if (!isAscending)
    throw new Error(ep + '\"times\" must be an increasing sequence of times');

  if (t0 === undefined)
    t0 = object.t0;

  if (t0 > times[0])  
    throw new Error(ep + '\"times\" the zero-time \"t0\"  must occur no later than the first observation');

  let ntimes = times.length;

  if (params === undefined) params = object.params;
  
  if (Object.keys(params).length === 0) 
    throw new Error(ep + '\"params\" must be supplied');
  

  params = [params]//as.matrix(params)
  let nrep = params.length;
  let paramnames = Object.keys(params[0]);
  if (paramnames.length <= 0)
    throw new Error(ep + '\"params\" must have rownames');

  let x0 = initState(object, params);
  let nvar = x0[0].length;
  statenames = Object.keys(x0[0]);
  // dim(x0) <- c(nvar,nrep,1)
  // dimnames(x0) <- list(statenames,NULL,NULL)

  let type = object.skeletonDetail.type;          // map or vectorfield?

  let x;
  if (type === "map") {

    // try {
      x = iterateMap(object, times, t0, x0, params);
    // } catch (error) {
    //   throw new Error(` ${ep} in map iterator: ${error}`)
    // }
      

  } else if (type=="vectorfield") {
    throw new Error("vectorfield is not translated.")
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
    // }

  //   x <- array(data=t(X[-1L,-1L]),dim=c(nvar,nrep,ntimes),
  //     dimnames=list(statenames,NULL,NULL))

  //   for (z in znames)
  //     for (r in seq_len(ncol(x)))
  //       x[z,r,-1] <- diff(x[z,r,])

  } else {

    throw new Erro(`${ep} deterministic skeleton has not been properly specified`);
    return;

  }

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

  return x
}