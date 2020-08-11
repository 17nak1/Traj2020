#library(pomp)
start_time <- Sys.time()
######################################################  Model Snippet
seasonal.sir.ode <- Csnippet("
                             double beta0;
                             double va, tt;
                             double seas;


                             // term-time seasonality
                             tt = (t-floor(t))*365.25;
                             if ((tt>=7&&tt<=100) || (tt>=115&&tt<=199) || (tt>=252&&tt<=300) || (tt>=308&&tt<=356))
                             seas = 1.0+amplitude*0.2411/0.7589;
                             else
                             seas = 1.0-amplitude;

                             // transmission rate
                             beta0 = R0*(gamma+mu)*(sigma+mu)/sigma;
                             // expected force of infection
                             double Beta = beta0*seas/pop;

                             va = 0;


                             DS =  birthrate*(1-va) -  Beta*S*I  - mu*S;
                             DE =  Beta*S*I - (sigma+mu)*E;
                             DI =  sigma*E - (gamma+mu)*I;
                             DR =  gamma*I  - mu*R + birthrate*va;
                             DH =  gamma*I;
                             ") #Seasonally forced snippet



# The above uses true incidence for reporting

initz <- Csnippet("
                  double m = pop/(S_0+E_0+I_0+R_0);
                  S = nearbyint(m*S_0);
                  E = nearbyint(m*E_0);
                  I = nearbyint(m*I_0);
                  R = nearbyint(m*R_0);
                  H = 0;
                  //printf(\"%f %f %f %f %f \\n\",S,E,R,I,H);

                  ")
# Sampling from the normal approximation of the binomial distribution
#' Measurement process, calling on both rmeasure and dmeasure on H
#' The Data is  binomial distributed
dmeas <- Csnippet("
                  double m = rho*H;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;

                  if (R_FINITE(cases)) {
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }

                  if (give_log) lik = log(lik);
                  } else {
                  lik = (give_log) ? 0 : 1;
                  }

                  ")


rmeas <- Csnippet("
                  double m = rho*H;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")




toEst <- Csnippet("
                  Tmu = log(mu);
                  Tsigma = log(sigma);
                  Tpsi = log(psi);
                  Tgamma = log(gamma);
                  TR0 = log(R0);
                  Trho = logit(rho);
                  Tamplitude = logit(amplitude);
                  ")

fromEst <- Csnippet("
                    Tmu = exp(mu);
                    Tsigma = exp(sigma);
                    Tpsi = exp(psi);
                    Tgamma = exp(gamma);
                    TR0 = exp(R0);
                    Trho = expit(rho);
                    Tamplitude = expit(amplitude);
                    ")

params_mod <- c("R0","amplitude","gamma","mu","sigma","rho","psi")

params_ic <- c("S_0","E_0","R_0","I_0")

statenames <- c("S","E","I","R","H")

zeronames <- c("H")
#########################################            data

London_BiData <- read.csv(file.path("London_BiData.csv"))
London_covar <- read.csv(file.path("London_covar.csv"))

#########################################       make pomp
pomp(
  data = London_BiData,
  times="time",
  t0=1940,#with(get(paste0(name,"_BiData")),2*time[1]-time[2]),
  skeleton=map(seasonal.sir.ode, delta.t = 1/365.25),
  rmeasure=rmeas,
  covar=London_covar, #covart,
  tcovar="time",
  dmeasure=dmeas,
  zeronames=zeronames,
  initializer=initz,
  toEstimationScale=toEst,
  fromEstimationScale=fromEst,
  statenames=statenames,
  paramnames=c(params_mod,params_ic)
) -> m1


params_fixed <- "R0"
params_fixed <- c("gamma","sigma",params_fixed)
params_mod_fit <- c("R0",  "mu", "rho" ,"psi")
params_ic_fit <- character(0)

ParamSetFile <- paste0("ParamSet_DeterministicSEIR_run28.csv")
full_set <- read.csv(file=ParamSetFile,header=T)
current_params <- unlist(full_set[1,])

coef(m1) <- c(current_params)

traj.match(m1,
           start=current_params,
           transform=T,
           est=c()) -> sets_traj
coef(sets_traj)
logLik(sets_traj)
end_time <- Sys.time()
end_time - start_time
