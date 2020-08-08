const { partrans } = require("./helpers")
const subplex = require("../subplex/run.js")

exports.minimInternal = function(objfun, start, est, object, method, transform,
  lower = NULL, upper = NULL, lb = lower, ub = upper)
{

  let ep = "In minimInternal:";
  if (Object.keys(start).length < 1)
    throw new Error(ep +"start must be supplied")
  let guess = {};
  if (transform) {
    start = partrans(object,start,dir="toEstimationScale")
    // if (is.null(names(start))||(!all(est%in%names(start))))
    //   throw new Error("est")," must refer to parameters named in ",
    //     sQuote("partrans(object,start,dir=\"toEstimationScale\")"),call.=FALSE)
    if (est.length > 0) guess = start.slice(...est)
  } else {
    // if (is.null(names(start))||(!all(est%in%names(start))))
    //   throw new Error("est")," must refer to parameters named in ",
    //     sQuote("start"),call.=FALSE)
    if (est.length > 0) guess = start.slice(...est);
  }

  if (est.length === 0) {

    val <- objfun(guess)
    conv = NA
    evals = [1,0];
    msg = "no optimization performed"

  } else {

    opts <- list(...)

    if (method == 'subplex') {
      opt = subplex(par=guess,fn=objfun,control=opts)
    } else if (method=="sannbox") {
      throw new Error ("Method 'sannbox' is not translated")
    } else if (method=="nloptr") {
      throw new Error ("Method 'nloptr' is not translated")
    } else {
      // opt <- optim(par=guess,fn=objfun,method=method,control=opts)
    }

    msg <- as.character(opt$message)

    if (method == "nloptr") {

      val <- opt$objective
      start[est] <- unname(opt$solution)
      conv <- opt$status
      evals <- opt$iterations

    } else {

      val <- opt$value
      start[est] <- unname(opt$par)
      conv <- opt$convergence
      evals <- opt$counts

    }
  }

  if (transform)
    start = partrans(object,start,dir="fromEstimationScale");

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
