const { partrans } = require("./helpers");
const subplex = require("../subplex/run.js");
const { trajMatchObjfun } = require("./trajMatchObjfun.js");

exports.minimInternal = function(objfun, start, est, object, method, transform,
  lower = null, upper = null, lb = lower, ub = upper, args)
{

  let ep = "In minimInternal:";
  if (Object.keys(start).length < 1)
    throw new Error(ep +"start must be supplied")
  let guess = {};
  if (transform) {
    start = partrans(object,[start],dir="toEstimationScale")[0]
    if (Object.keys(start).some(a => {a === null}))//||(!all(est%in%names(start)))
      throw new Error("'est' must refer to parameters named in partrans(object,start,dir=toEstimationScale")
  } else {
    if (Object.keys(start).some(a => {a === null}))// ||(!all(est%in%names(start)))
      throw new Error("'est' must refer to parameters named in start.")
  }
  if (est.length > 0) guess = start.slice(...est);

  let val;
  if (est.length === 0) {
    val = trajMatchObjfun(object=object, params=start, est=guess, transform=transform)
    conv = NA
    evals = [1,0];
    msg = "no optimization performed"

  } else {

    let opts = args

    if (method == 'subplex') {
      opt = subplex(par=guess,fn=objfun,control=opts)
    } else if (method=="sannbox") {
      throw new Error ("Method 'sannbox' is not translated")
    } else if (method=="nloptr") {
      throw new Error ("Method 'nloptr' is not translated")
    } else {
      // opt <- optim(par=guess,fn=objfun,method=method,control=opts)
    }

    // msg <- as.character(opt$message)

    // if (method == "nloptr") {

    //   val <- opt$objective
    //   start[est] <- unname(opt$solution)
    //   conv <- opt$status
    //   evals <- opt$iterations

    // } else {

    //   val <- opt$value
    //   start[est] <- unname(opt$par)
    //   conv <- opt$convergence
    //   evals <- opt$counts

    // }
  }

  if (transform) start = partrans(object, [start], dir = "fromEstimationScale")[0];

  return {
    params : start,
    est : est,
    transform : transform,
    value : val,
    convergence : parseInt(conv),
    evals : parseInt(evals),
    msg : msg
  }
}
