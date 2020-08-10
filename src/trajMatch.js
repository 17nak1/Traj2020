
const { coef, partrans } = require("./helpers");
const { trajectory } = require("trajectory");
const snippet = require("exampleJs/modelSnippet.js");
const { minimInternal } = require("minimInternal.js");

exports.trajMatch = function (object, start, est = [],
   method = c("Nelder-Mead","subplex","SANN","BFGS", "sannbox","nloptr"),transform = false) {
  
  if (!start) start = coef(object);

  let m = minimInternal(
     objfun=trajMatchObjfun(
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

const trajMatchObjfun  = function (object, params, est, transform = FALSE, args) {
  tmofInternal(
    object=object,
    params=params,
    est=est,
    transform=transform,
    args
  )
}

const tmofInternal = function (object, params, est, transform, args) {

  ep = "in traj.match.objfun : "
  
  if (!est) est = [];

  if (Object.keys(params).length === 0) params = coef(object);
  
  // if ((!is.numeric(params))||(is.null(names(params))))
  //   throw new Error(ep+ "'params' must be a named numeric vector");
  if (transform) {
    params = partrans(object,params,dir="toEstimationScale");
  }
    
  // it does match(est,names(params))
  let parEstIdx = []; 
  for (let i = 0; i < est.length; i++) {
    for (let j = 0; j < Object.keys(params).length; j++) {
      if(est[i] === Object.keys(params)[j]) parEstIdx.push(j);
    }
  }

  if (parEstIdx.some(a => {return a === NaN}))
    throw new Error(ep + "est does not match with parameters")

  return function (par) {
    let d;
    params[parEstIdx] = par
    if (transform) tparams = partrans(object,params,dir="fromEstimationScale")
    d = snippet.dmeasure(
      object,
      y=object.data,
      x=trajectory(
        object,
        params = transform? tparams : params,
        args
      ),
      times = object.times,
      params = transform? tparams : params,
      log = true
    )
    return -d.reduce((a, b) => a + b, 0);
  }
}
