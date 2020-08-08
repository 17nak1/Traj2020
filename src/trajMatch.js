
const { coef, partrans } = require("./helpers")
const tmofInternal = function (object, params, est, transform) {

  const ep = "In traj.match.objfun";
  if (!est) est = [];
  // est <- as.character(est)

  if (!params) params = coef(object);
  // if ((!is.numeric(params))||(is.null(names(params))))
  //   stop(ep,sQuote("params")," must be a named numeric vector",call.=FALSE)
  if (transform)
    params = partrans(object,params,dir="toEstimationScale");
  let parEstIdx = match(est,names(params))//????? it returns index, but maybe in the key of object order Changed
  if (any(is.na(par.est.idx)))
    stop(ep,"parameter(s): ",
         paste(sapply(est[is.na(par.est.idx)],sQuote),collapse=","),
         " not found in ",sQuote("params"),call.=FALSE)

  function (par) {
    pompLoad(object)
    params[par.est.idx] <- par
    if (transform)
      tparams <- partrans(object,params,dir="fromEstimationScale")
    d <- dmeasure(
      object,
      y=object@data,
      x=trajectory(
        object,
        params=if (transform) tparams else params,
        ...
      ),
      times=time(object),
      params=if (transform) tparams else params,
      log=TRUE
    )
    pompUnload(object)
    -sum(d)
  }
}


const trajMatch = function (object, start, est = [],
            method = c("Nelder-Mead","subplex","SANN","BFGS", "sannbox","nloptr"),transform = false)
{

  if (!start) start = coef(object);

  let m = minim.internal(
    objfun=traj.match.objfun(
      object=object,
      params=start,
      est=est,
      transform=transform
    ),
    start=start,
    est=est,
    object=object,
    method=method,
    transform=transform
  )

  // fill params slot appropriately
  coef(object) = m.params

  // fill states slot appropriately
  x = trajectory(object)
  // object.states = array(data=x,dim=dim(x)[c(1L,3L)])
  // rownames(object@states) <- rownames(x)

  return {
    // "traj.matched.pomp",
    object,
    transform: transform,
    est: est,
    value: -m.value,
    evals: m.evals,
    convergence: m.convergence,
    msg: m.msg
  }
}


setMethod(
  "traj.match",
  signature=signature(object="traj.matched.pomp"),
  function (object, est, transform, ...)
  {
    if (missing(est)) est <- object@est
    if (missing(transform)) transform <- object@transform

    f <- selectMethod("traj.match","pomp")

    f(
      object=object,
      est=est,
      transform=transform,
      ...
    )
  }
)
