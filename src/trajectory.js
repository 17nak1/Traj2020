
// exports.trajectory = function(){return {S:1, E:0, I:0, R:0, H:0}}
// -*- C++ -*-

// #include <Rdefines.h>
// #include <R_ext/Constants.h>

// #include "pomp_internal.h"

const iterate_map_native =function (X, time, p, deltat, t, x,
   ntimes, nvars, npars, ncovars, nzeros, nreps, sidx, pidx, cidx, zidx, covar_table, ff, args)
{
  let covars = null;
  int nsteps;
  let *Xs, *xs, *ps;
  int h, i, j, k;
  set_pomp_userdata(args);
  if (ncovars > 0) covars = (double *) Calloc(ncovars, double);
  for (k = 0; k < ntimes; k++, time++, X += nvars*nreps) {
    R_CheckUserInterrupt();
    for (i = 0; i < nzeros; i++)
      for (j = 0, xs = &x[zidx[i]]; j < nreps; j++, xs += nvars)
	*xs = 0.0;
    nsteps = num_map_steps(t,*time,deltat);
    for (h = 0; h < nsteps; h++) {
      table_lookup(covar_table,t,covars);
      for (j = 0, Xs = X, xs = x, ps = p; j < nreps; j++, Xs += nvars, xs += nvars, ps += npars) {
	(*ff)(Xs,xs,ps,sidx,pidx,cidx,ncovars,covars,t);
      }
      memcpy(x,X,nvars*nreps*sizeof(double));
      t += deltat;
    }
    if (nsteps == 0) memcpy(X,x,nvars*nreps*sizeof(double));
  }
  if (ncovars > 0) Free(covars);
  unset_pomp_userdata();
}


exports.iterate_map = function (object, times, t0, x0, params)
{
  let nprotect = 0;
  let fn;
  let X;
  let Snames, Pnames, Cnames;
  let pompfun;
  let zeronames;
  // let *zidx = 0;
  let nvars, npars, nreps, ntimes, ncovars, nzeros;

  let deltat = object.skeletonDetail.deltaT;
  let t = object.t0;

  nvars = x0[0].length; nreps = x0.length;

  npars = params[0].length;

  if (nreps !== params.length)
    throw new Error("in 'trajectory': dimension mismatch between 'x0' and 'params'"); // # nocov

  
  let ntimes = times.length;

  // PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x0))); nprotect++;
  // PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  // PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(GET_SLOT(object,install("covar"))))); nprotect++;

  // set up the covariate table
  // covariate_table = make_covariate_table(object,&ncovars);

  // extract user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("skeleton"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // get names and indices of accumulator variables
  let zeronames = object.zeronames;
  let nzeros = zeronames.length;
  if (nzeros > 0) {
    zidx = INTEGER(PROTECT(matchnames(Snames,zeronames,"state variables"))); nprotect++;
  }

  // create array to store results
  {
    let dim = [nvars, nreps, ntimes];
    PROTECT(X = makearray(3,dim)); nprotect++;
    setrownames(X,Snames,3);
  }
  // let X = new Array();
  // set up the computations
  
  // int *sidx, *pidx, *cidx;
  // pomp_skeleton *ff;
  *((void **) (&ff)) = R_ExternalPtrAddr(fn);
  // construct state, parameter, covariate indices
  sidx = INTEGER(PROTECT(name_index(Snames,pompfun,"statenames","state variables"))); nprotect++;
  pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames","parameters"))); nprotect++;
  cidx = INTEGER(PROTECT(name_index(Cnames,pompfun,"covarnames","covariates"))); nprotect++;

  iterate_map_native(X, times, params, deltat, t, x0, ntimes, nvars, npars, ncovars, nzeros, nreps,
    sidx, pidx, cidx, zidx, &covariate_table, ff, args);


  UNPROTECT(nprotect);
  return X;
}

