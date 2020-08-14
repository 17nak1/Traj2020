
// exports.trajectory = function(){
//   let v = {S:1, E:0, I:0, R:0, H:0};
//   return new Array(548).fill(0).map(a => v)
// }


exports.iterateMap = function (object, times, t0, x0, params) {
  // let fn;
  // let X;
  // let pompfun;
  // let zeronames;
  // let *zidx = 0;

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

  // // set up the covariate table
  // // covariate_table = make_covariate_table(object,&ncovars);

  // extract user-defined function
  //pompfun = trajMatchObjfun;
  // PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // // extract 'userdata' as pairlist
  // PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // get names and indices of accumulator variables
  let zeronames = object.zeronames;
  let nzeros = zeronames.length;
  // if (nzeros > 0) {
  //   zidx = INTEGER(PROTECT(matchnames(Snames,zeronames,"state variables"))); nprotect++;
  // }

  // // create array to store results
  // {
  //   let dim = [nvars, nreps, ntimes];
  //   PROTECT(X = makearray(3,dim)); nprotect++;
  //   setrownames(X,Snames,3);
  // }
  let X = [];
  for(let i = 0; i < ntimes; i++){
    a = [{...(x0[0])}];
    X.push(a);
  }
  // // set up the computations
  
  // // int *sidx, *pidx, *cidx;
  let ff = object.skeleton;
  // *((void **) (&ff)) = R_ExternalPtrAddr(fn);
  // // construct state, parameter, covariate indices
  // sidx = INTEGER(PROTECT(name_index(Snames,pompfun,"statenames","state variables"))); nprotect++;
  // pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames","parameters"))); nprotect++;
  // cidx = INTEGER(PROTECT(name_index(Cnames,pompfun,"covarnames","covariates"))); nprotect++;

  iterate_map_native(X, times, params, deltat, t, x0, ff, object);


  // UNPROTECT(nprotect);
  return X;
}

const iterate_map_native = function(X, times, params, deltat, t0, x0, func, object){  
  let zeronames = object.zeronames; 
  if (deltat <= 0)
    throw new Error("In euler.js: 'delta.t' should be a positive number");
  let nvars = Object.keys(X[0][0]).length;
  let nreps = X.length;
  let npars = Object.keys(params[0]).length;
  let ntimes = times.length;
  let t = t0;
  let xt = X;
  let args = object.globals;
  for (let step = 0; step < ntimes; step++) {
    
    // set accumulator variables to zero
    for (j = 0; j < nreps; j++) {
      for (i = 0; i < zeronames.length; i++) {
         xt[j][0][zeronames[i]] = 0;
      } 
    }
    dt = deltat;
    nstep = numMapSteps(t, times[step], dt);      
    for (let k = 0; k < nstep; k++) { // loop over Euler steps
      let  interpolatorObj = object.interpolator(t);
      for (let j = 0 ; j < nreps; j++) { // loop over replicates
         xt[j][0] = func(xt[j][0], params[0], t, dt, interpolatorObj,args);
      }
      t += dt;

      if ((k == nstep-2)) { // penultimate step
        dt = times[step] - t;
        t = times[step] - dt;
      }  
    }
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