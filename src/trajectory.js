
exports.iterateMap = function (object, times, t0, x0, params) {

  let deltat = object.skeletonDetail.deltaT;
  let t = t0;

  let nvars = Object.keys(x0[0]).length;
  let nreps = x0.length;

  let npars = Object.keys(params[0]).length;

  if (nreps !== params.length)
    throw new Error("in 'trajectory': dimension mismatch between 'x0' and 'params'"); // # nocov

  
  let ntimes = times.length;

  let Snames = Object.keys(x0[0]);
  let Pnames = Object.keys(params[0]);
  let Cnames = object.covarnames;

  let zeronames = object.zeronames;
  let nzeros = zeronames.length;
  
  let X = [];
  for(let i = 0; i < ntimes; i++){
    a = [{...(x0[0])}];
    X.push(a);
  }
  // set up the computations
  let ff = object.skeleton;
  iterate_map_native(X, times, params, deltat, t, x0, ff, object);

  return X;
}

const iterate_map_native = function(X, times, params, deltat, t0, x0, func, object){  
  let zeronames = object.zeronames; 
  if (deltat <= 0)
    throw new Error("In euler.js: 'delta.t' should be a positive number");
  let nvars = Object.keys(X[0][0]).length;
  let nreps = X[0].length;
  let npars = Object.keys(params[0]).length;
  let ntimes = times.length;
  let t = t0;
  let xt = X;
  let args = object.globals;
  for (let step = 0; step < ntimes; step++) {
    
    // set accumulator variables to zero
    for (j = 0; j < nreps; j++) {
      for (i = 0; i < zeronames.length; i++) {
         x0[j][zeronames[i]] = 0;
      } 
    }
    dt = deltat;
    nstep = numMapSteps(t, times[step], dt);      
    for (let k = 0; k < nstep; k++) { // loop over Euler steps
      let  interpolatorObj = object.interpolator(t);
      for (let j = 0 ; j < nreps; j++) { // loop over replicates
         x0[j] = func(x0[j], params[0], t, dt, interpolatorObj,args);
      }
      t = (t * 10 + dt * 10)/10;
    }
    xt[step][0] = {...x0[0]};
  }
  
  return xt;

}

const numMapSteps = function (t1, t2, dt) {
  let DOUBLE_EPS = 10e-8
  let tol = Math.sqrt(DOUBLE_EPS)
  let nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 - tol))
  return (nstep > 0) ? nstep : 0
}