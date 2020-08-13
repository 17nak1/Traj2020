
snippet = {}
//let mathLib = require('./mathLib')

snippet.endTime = "2020-07-16";
snippet.T0 = 75;
snippet.T1 = 139;
snippet.modeltype = { nstageE: 3, nstageI: 3, nstageP: 3, nstageH: 3, nstageC: 3, nstageV: 3 };
snippet.pop = 10e6;

snippet.dObs = function(states, params, covar, t, dt, args){
  let HOSP = params.H1;////////////////////////////////////////// FIXME: change to array
  let CARE = params.C1;//////////////////////////////////////////
  let VENT = params.V1;//////////////////////////////////////////
  let lik_reports, lik_deaths, lik_hospital, lik_ICU, lik_ventilator;
  let reported_cases, reported_deaths;
  let TOT_H, TOT_C, TOT_V;
  let i;
  if (isFinite(reports)) {
    reported_cases = states.casesH + states.casesIQ;
    lik_reports = 0;///////////////////////////////////////////// TODO: dpois(reports, reported_cases + 1e-6, 1);
  } else {
    lik_reports = 0;
  }
  if (isFinite(deaths)) {
    reported_deaths = states.deathsCV + states.deathsIIQ;
    lik_deaths = 0;//////////////////////////////////////////////dpois(deaths, reported_deaths + 1e-6, 1);
  } else {
    lik_deaths = 0;
  }
  for(i = 0, TOT_H = 0; i < args.nstageH; i++) {
    TOT_H += HOSP[i];
  }
  for(i = 0, TOT_C = 0; i < args.nstageC; i++) {
    TOT_C += CARE[i];
  }
  for(i = 0, TOT_V = 0; i < args.nstageV; i++) {
    TOT_V += VENT[i];
  }
  if (isFinite(hospital) && isFinite(ICU)) {
    lik_hospital = 0;////////////////////////////////////dpois(hospital - ICU, TOT_H + 1e-6, 1);
  } else {
    lik_hospital = 0;
  }
  if (isFinite(ICU) && isFinite(ventilator)) {
    lik_ICU = 0;/////////////////////////////////////dpois(ICU - ventilator, TOT_C + 1e-6, 1);
  } else {
    lik_ICU = 0;
  }
  if (isFinite(ventilator)) {
    lik_ventilator = 0;//////////////////////////////////dpois(ventilator, TOT_V + 1e-6, 1);
  } else {
    lik_ventilator = 0;
  }
  lik = lik_reports + lik_deaths + lik_hospital + lik_ICU + lik_ventilator;
  if (args.give_log == 0 ) {
    lik = exp(lik);
  }
}

snippet.rObs = function(states, params, covar, t, dt, args){
  let HOSP = params.H1;
  let CARE = params.C1;
  let VENT = params.V1;
  let TOT_H, TOT_C, TOT_V;
  let reported_cases = casesH + casesIQ;
  reports = 0;///////////// TODO: rpois(reported_cases  + 1e-6);
  
  let reported_deaths = deathsCV + deathsIIQ;
  deaths = 0;///////////////rpois(reported_deaths  + 1e-6);
  
  let i;
  for(i = 0, TOT_H = 0; i < args.nstageH; i++) {
    TOT_H += HOSP[i];
  }
  for(i = 0, TOT_C = 0; i < args.nstageC; i++) {
    TOT_C += CARE[i];
  }                 
  for(i = 0, TOT_V = 0; i < args.nstageV; i++) {
    TOT_V += VENT[i];
  }
  ventilator = 0;////////////////////////////////rpois(TOT_V + 1e-6);
  ICU = 0;/////////////////////////////////rpois(TOT_C + 1e-6) + ventilator;
  hospital = 0;///////////////////////////////rpois(TOT_H + 1e-6) + ICU;
}

snippet.rSim = function(states, params, covar, t, dt, args){
  let SUSC = params.S;
  let EXPD = params.E1;
  let PRE = params.P1;
  let INFD = params.I1;
  let HOSP = params.H1;
  let CARE = params.C1;
  let VENT = params.V1;
  let DEAD = params.M;
  let RCVD = params.R;
  
  let EXPDQ =params.EQ1;
  let PREQ = params.PQ1;
  let INFDQ =params.IQ1;
  
  
  // Different transmission rates
  
  let dQdt;
  let TOT_PRE, TOT_INFD;
  let i;
  
  for(i = 0, TOT_PRE = 0; i < args.nstageP; i++) {
    TOT_PRE += PRE[i];
  }
  for(i = 0, TOT_INFD = 0; i < args.nstageI; i++) {
    TOT_INFD += INFD[i];
  }
  
  let PD, lambdaI, lambdaP, lambdaPQ, lambda, lambdaQ;
  if ( isFinite(tests) ) {
    PD = rho * (tests / (tests + TF));
  } else {
    PD = rho * (15.0 / (15.0 + TF)) ;
  }
  lambdaI = betaI * TOT_INFD;
  lambdaP = betaI * theta * TOT_PRE * (1-PD);
  lambdaPQ = betaI * theta * TOT_PRE * PD;


  if (t < T0) {
    dQdt  = 0;///////////////////////////////////// TODO: rgammawn(beta_sd, dt)/dt;
    lambda = ( (lambdaI + lambdaP + iota) / pop ) * dQdt;
    lambdaQ = ( (lambdaPQ) / pop ) * dQdt;
  } else if (t < T0 + 7.0) {
    let x = (t-T0) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = 0;//////////////////////////rgammawn(((1-ss) + ss*dB0)*beta_sd, dt)/dt;
    lambda = ( ( ((1-ss) + ss*dI0) * lambdaI + ((1-ss) + ss*dP0) * lambdaP + ((1-ss) + ss*dT0)*iota ) / pop ) * dQdt;
    lambdaQ = ( ( ((1-ss) + ss*dP0) * lambdaPQ  ) / pop ) * dQdt;
  } else if (t < T1) {
    dQdt  = 0;///////////////////////////////rgammawn(dB0*beta_sd, dt)/dt;
    lambda = ( ( dI0 * lambdaI + dP0 * lambdaP + dT0 * iota ) / pop ) * dQdt;
    lambdaQ = ( (  dP0 * lambdaPQ  ) / pop ) * dQdt;
  } else if (t < T1 +7.0) {
    let x = (t-T1) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = 0;//////////////////////////rgammawn(((1-ss)*dB0 + ss*dB1)*beta_sd, dt)/dt;
    lambda = ( ( ((1-ss)*dI0 + ss*dI1) * lambdaI + ((1-ss)*dP0 + ss*dP1) * lambdaP + ((1-ss)*dT0 + ss*dT1)*iota ) / pop ) * dQdt;
    lambdaQ = ( ( ((1-ss)*dP0 + ss*dP1) * lambdaPQ  ) / pop ) * dQdt;
  } else {
    dQdt  = 0;/////////////////////////////rgammawn(dB1*beta_sd, dt)/dt;
    lambda = ( ( dI1 * lambdaI + dP1 * lambdaP + dT1 * iota ) / pop ) * dQdt;
    lambdaQ = ( (  dP1 * lambdaPQ  ) / pop ) * dQdt;
  }

  // From class S
  let transS = Array(2);
  let rateS = Array(2);
  rateS[0] = lambda ;
  rateS[1] = lambdaQ;
  // TODO: reulermultinom(2,SUSC[0], &rateS[0], dt, &transS[0]);
  
  
  // From class EQ
  let transEQ = array(args.nstageE);
  let rateEQ = nstageE * sigma;
  for (i = 0; i < args.nstageE; i++) {
    0;////////// reulermultinom(1, EXPDQ[i], &rateEQ, dt, &transEQ[i]);
  }
  
  // From class PQ
  let transPQ = array(args.nstageP + 1);
  let ratePQ = nstageP * kappa;
  for (i = 0; i < args.nstageP-1; i++) {
    0;//////////////reulermultinom(1, PREQ[i], &ratePQ, dt, &transPQ[i]);
  }
  let ratePQIQH = array(2);
  ratePQIQH[0] = (1-qP) * args.nstageP * kappa;
  ratePQIQH[1] = qP * args.nstageP * kappa;
  0;////////////////reulermultinom(2, PREQ[args.nstageP-1], &ratePQIQH[0], dt, &transPQ[nstageP-1]);
  
  // From class IQ
  let transIQ = array(args.nstageI + 1);
  let rateIQ = args.nstageI * gammaI;
  for (i = 0; i < args.nstageI-1; i++) {
    0;/////////////////reulermultinom(1, INFDQ[i], &rateIQ, dt, &transIQ[i]);
  }
  
  let rateIQRD = array(2);
  rateIQRD[0] = (1-mI) * args.nstageI * gammaI;
  rateIQRD[1] = mI * args.nstageI * gammaI;
  0;////////////////////reulermultinom(2, INFDQ[nstageI-1], &rateIQRD[0], dt, &transIQ[nstageI-1]);
  
  // From class E
  let transE = array(nstageE);
  let rateE = args.nstageE * sigma;
  for (i = 0; i < args.nstageE; i++) {
    0;/////////////////reulermultinom(1, EXPD[i], &rateE, dt, &transE[i]);
  }
  
  // From class P
  let transP = array(nstageP + 2);
  let rateP = args.nstageP * kappa;
  for (i = 0; i < nstageP-1; i++) {
    0;///////////////reulermultinom(1, PRE[i], &rateP, dt, &transP[i]);
  }
  
  let ratePIHIQ = array(3);
  ratePIHIQ[0] = (1.0-PD) * (1.0-qP) * args.nstageP * kappa;
  ratePIHIQ[1] = qP * args.nstageP * kappa;
  ratePIHIQ[2] = PD * (1.0-qP) * args.nstageP * kappa;
  0;///////////////reulermultinom(3, PRE[nstageP-1], &ratePIHIQ[0], dt, &transP[nstageP-1]);
  
  // From class I
  let transI = array(args.nstageI + 1);
  let rateI = args.nstageI * gammaI;
  for (i = 0; i < args.nstageI-1; i++) {
    0;/////////////////reulermultinom(1, INFD[i], &rateI, dt, &transI[i]);
  }
  
  let rateIRD = array(2);
  rateIRD[0] = (1-mI) * args.nstageI * gammaI;
  rateIRD[1] = mI * args.nstageI * gammaI;
  0;/////////////////reulermultinom(2, INFD[nstageI-1], &rateIRD[0], dt, &transI[nstageI-1]);
  
  
  // From class H
  let transH = array(args.nstageH + 1);
  let rateH = args.nstageH * gammaH;
  for (i = 0; i < args.nstageH-1; i++) {
    0;/////////////reulermultinom(1, HOSP[i], &rateH, dt, &transH[i]);
  }
  let rateHRC = array(2);
  rateHRC[0] = (1-qH) * args.nstageH * gammaH;
  rateHRC[1] = qH * args.nstageH * gammaH;
  0;///////////////////////reulermultinom(2, HOSP[nstageH-1], &rateHRC[0], dt, &transH[nstageH-1]);
  
  
  // From class C
  let transC = array(args.nstageC + 2);
  let rateC = args.nstageC * gammaC;
  for (i = 0; i < args.nstageC-1; i++) {
    0;///////////////////////reulermultinom(1, CARE[i], &rateC, dt, &transC[i]);
  }
  let rateCRVM = array(3);
  rateCRVM[0] = (1-mC) * (1-qC) * args.nstageC * gammaC;
  rateCRVM[1] = qC * args.nstageC * gammaC;
  rateCRVM[2] = mC * (1-qC) * args.nstageC * gammaC;
  0;/////////////reulermultinom(3, CARE[nstageC-1], &rateCRVM[0], dt, &transC[nstageC-1]);
  
  // From class V
  let transV = array(args.nstageV + 1);
  let rateV = args.nstageV * gammaV;
  for (i = 0; i < args.nstageV-1; i++) {
    0;///////////////////reulermultinom(1, VENT[i], &rateV, dt, &transV[i]);
  }
  let rateVRD = array(2);
  rateVRD[0] = (1-mV) * args.nstageV * gammaV;
  rateVRD[1] = mV * args.nstageV * gammaV;
  0;/////////////////reulermultinom(2, VENT[nstageV-1], &rateVRD[0], dt, &transV[nstageV-1]);
  
  
  
  // Balance the equations
  SUSC[0] -= transS[0];
  EXPD[0] += transS[0];
  for (i = 0; i < args.nstageE; i++) EXPD[i] -= transE[i];
  for (i = 1; i < args.nstageE; i++) EXPD[i] += transE[i-1];
  PRE[0] += transE[args.nstageE-1];
  for (i = 0; i < args.nstageP; i++) PRE[i] -= transP[i];
  for (i = 1; i < args.nstageP; i++) PRE[i] += transP[i-1];
  INFD[0] += transP[args.nstageP-1];
  HOSP[0] += transP[args.nstageP];
  PRE[args.nstageP-1] -= transP[args.nstageP];
  for (i = 0; i < args.nstageI; i++) INFD[i] -= transI[i];
  for (i = 1; i < args.nstageI; i++) INFD[i] += transI[i-1];
  for (i = 0; i < args.nstageH; i++) HOSP[i] -= transH[i];
  for (i = 1; i < args.nstageH; i++) HOSP[i] += transH[i-1];
  CARE[0] += transH[args.nstageH];
  HOSP[args.nstageH-1] -= transH[args.nstageH];
  INFD[args.nstageI-1] -= transI[args.nstageI];
  for (i = 0; i < args.nstageC; i++) CARE[i] -= transC[i];
  for (i = 1; i < args.nstageC; i++) CARE[i] += transC[i-1];
  VENT[0] += transC[args.nstageC];
  CARE[args.nstageC-1] -= transC[args.nstageC] + transC[args.nstageC+1];
  for (i = 0; i < args.nstageV; i++) VENT[i] -= transV[i];
  for (i = 1; i < args.nstageV; i++) VENT[i] += transV[i-1];
  VENT[args.nstageV-1] -= transV[args.nstageV];               
  RCVD[0] += transI[args.nstageI-1] + transH[args.nstageH-1] + transC[args.nstageC-1] + transV[args.nstageV-1];
  DEAD[0] += transI[args.nstageI] + transC[args.nstageC+1] + transV[args.nstageV];
  
  SUSC[0] -= transS[1];
  EXPDQ[0] += transS[1];
  for (i = 0; i < args.nstageE; i++) EXPDQ[i] -= transEQ[i];
  for (i = 1; i < args.nstageE; i++) EXPDQ[i] += transEQ[i-1];
  PREQ[0] += transEQ[args.nstageE-1];
  for (i = 0; i < args.nstageP; i++) PREQ[i] -= transPQ[i];
  for (i = 1; i < args.nstageP; i++) PREQ[i] += transPQ[i-1];
  PREQ[args.nstageP-1] -= transPQ[args.nstageP];  
  PRE[args.nstageP-1] -= transP[args.nstageP+1];
  INFDQ[0] += transPQ[args.nstageP-1] + transP[args.nstageP+1];
  HOSP[0] += transPQ[args.nstageP];
  for (i = 0; i < args.nstageI; i++) INFDQ[i] -= transIQ[i];
  for (i = 1; i < args.nstageI; i++) INFDQ[i] += transIQ[i-1];
  INFDQ[args.nstageI-1] -= transIQ[args.nstageI];
  RCVD[0] += transIQ[args.nstageI-1];
  DEAD[0] += transIQ[args.nstageI];
  
  casesI += transP[args.nstageP-1];
  casesIQ += transPQ[args.nstageP-1] + transP[args.nstageP+1];
  casesH += transPQ[args.nstageP] + transP[args.nstageP];
  deathsIIQ += transI[args.nstageI] + transIQ[args.nstageI];
  deathsCV += transC[args.nstageC+1] + transV[args.nstageV];
}

snippet.skel = function(states, params, covar, t, dt, args){
  dt = 0.1;
  let  SUSC = params.S;
  let  EXPD = params.E1;
  let  PRE = params.P1;
  let  INFD = params.I1;
  let  HOSP = params.H1;
  let  CARE = params.C1;
  let  VENT = params.V1;
  let  DEAD = params.M;
  let  RCVD = params.R;
  
  let  DSUSC = params.DS;
  let  DEXPD = params.DE1;
  let  DPRE = params.DP1;
  let  DINFD = params.DI1;
  let  DHOSP = params.DH1;
  let  DCARE = params.DC1;
  let  DVENT = params.DV1;
  let  DDEAD = params.DM;
  let  DRCVD = params.DR;
  
  
  let  EXPDQ = params.EQ1;
  let  PREQ = params.PQ1;
  let  INFDQ = params.IQ1;
  let  DEXPDQ = params.DEQ1;
  let  DPREQ = params.DPQ1;
  let  DINFDQ = params.DIQ1;
  
  
  // Different transmission rates
  
  let  dQdt;
  let  TOT_INFD, TOT_PRE;
  let i;
  
  for(i = 0, TOT_PRE = 0; i < argsnstageP; i++) {
    TOT_PRE += PRE[i];
  }
  for(i = 0, TOT_INFD = 0; i < argsnstageI; i++) {
    TOT_INFD += INFD[i];
  }
  
  
  let  PD, lambdaI, lambdaP, lambdaPQ, lambda, lambdaQ;
  if ( isFinite(tests) ) {
    PD = rho * (tests / (tests + TF));
  } else {
    PD = rho * (15.0 / (15.0 + TF)) ;
  }
  lambdaI = betaI * TOT_INFD;
  lambdaP = betaI * theta * TOT_PRE * (1-PD);
  lambdaPQ = betaI * theta * TOT_PRE * PD;
  
  if (t < T0) {
    dQdt  = 0;/////////////// TODO: rgammawn(beta_sd, dt)/dt;
    lambda = ( (lambdaI + lambdaP + iota) / pop ) * dQdt;
    lambdaQ = ( (lambdaPQ) / pop ) * dQdt;
  } else if (t < T0 + 7.0) {
    let  x = (t-T0) / 7.0;
    let  ss = 3*x*x - 2*x*x*x;
    dQdt  = 0;/////////////// rgammawn(((1-ss) + ss*dB0)*beta_sd, dt)/dt;
    lambda = ( ( ((1-ss) + ss*dI0) * lambdaI + ((1-ss) + ss*dP0) * lambdaP + ((1-ss) + ss*dT0)*iota ) / pop ) * dQdt;
    lambdaQ = ( ( ((1-ss) + ss*dP0) * lambdaPQ  ) / pop ) * dQdt;
  } else if (t < T1) {
    dQdt  = 0;/////////////// rgammawn(dB0*beta_sd, dt)/dt;
    lambda = ( ( dI0 * lambdaI + dP0 * lambdaP + dT0 * iota ) / pop ) * dQdt;
    lambdaQ = ( (  dP0 * lambdaPQ  ) / pop ) * dQdt;
  } else if (t < T1 +7.0) {
    let  x = (t-T1) / 7.0;
    let  ss = 3*x*x - 2*x*x*x;
    dQdt  = 0;/////////////// rgammawn(((1-ss)*dB0 + ss*dB1)*beta_sd, dt)/dt;
    lambda = ( ( ((1-ss)*dI0 + ss*dI1) * lambdaI + ((1-ss)*dP0 + ss*dP1) * lambdaP + ((1-ss)*dT0 + ss*dT1)*iota ) / pop ) * dQdt;
    lambdaQ = ( ( ((1-ss)*dP0 + ss*dP1) * lambdaPQ  ) / pop ) * dQdt;
  } else {
    dQdt  = 0;/////////////// rgammawn(dB1*beta_sd, dt)/dt;
    lambda = ( ( dI1 * lambdaI + dP1 * lambdaP + dT1 * iota ) / pop ) * dQdt;
    lambdaQ = ( (  dP1 * lambdaPQ  ) / pop ) * dQdt;
  }
  
  // From class S
  let  transS = array(2);
  let  rateS = array(2);
  rateS[0] = lambda ;
  rateS[1] = lambdaQ;
  let  rateS_tot = rateS[0] + rateS[1];
  let  transS_tot = (1.0 - exp(- rateS_tot * dt)) * SUSC[0];
  transS[0] = rateS[0] / rateS_tot * transS_tot;
  transS[1] = rateS[1] / rateS_tot * transS_tot;
  
  // From class EQ
  let  transEQ = array(args.nstageE);
  let  rateEQ = args.nstageE * sigma;
  for (i = 0; i < args.nstageE; i++) {
    transEQ[i] =  (1.0 - exp(- rateEQ * dt))* EXPDQ[i];
  }
  
  // From class PQ
  let  transPQ = array(args.nstageP + 1);
  let  ratePQ = args.nstageP * kappa;
  for (i = 0; i < args.nstageP-1; i++) {
    transPQ[i] =  (1.0 - exp(- ratePQ * dt))* PREQ[i];
  }
  let  transPQIQH;
  transPQIQH = (1.0 - exp(- ratePQ * dt))* PREQ[args.nstageP-1];
  transPQ[args.nstageP-1] = (1-qP) * transPQIQH;
  transPQ[args.nstageP] = qP * transPQIQH;
  
  // From class IQ
  let  transIQ = array(args.nstageI + 1);
  let  rateIQ = args.nstageI * gammaI;
  for (i = 0; i < args.nstageI-1; i++) {
    transIQ[i] =  (1.0 - exp(- rateIQ * dt))* INFDQ[i];
  }
  let  transIQRD;
  transIQRD = (1.0 - exp(- rateIQ * dt))* INFDQ[args.nstageI-1];
  transIQ[args.nstageI-1] = (1-mI) * transIQRD;
  transIQ[args.nstageI] = mI * transIQRD;
  
  
  // From class E
  let  transE = array(args.nstageE);
  let  rateE = args.nstageE * sigma;
  for (i = 0; i < args.nstageE; i++) {
    transE[i] =  (1.0 - exp(- rateE * dt))* EXPD[i];
  }
  
  // From class P
  let  transP = array(args.nstageP + 2);
  let  rateP = args.nstageP * kappa;
  for (i = 0; i < args.nstageP-1; i++) {
    transP[i] =  (1.0 - exp(- rateP * dt))* PRE[i];
  }
  let  transPIHIQ;
  transPIHIQ = (1.0 - exp(- rateP * dt))* PRE[args.nstageP-1];
  transP[args.nstageP-1] = (1-PD) * (1-qP) * transPIHIQ;
  transP[args.nstageP] = qP * transPIHIQ;
  transP[args.nstageP+1] = PD * (1-qP) * transPIHIQ;
  
  // From class I
  let  transI = array(args.nstageI + 1);
  let  rateI = args.nstageI * gammaI;
  for (i = 0; i < args.nstageI-1; i++) {
    transI[i] =  (1.0 - exp(- rateI * dt))* INFD[i];
  }
  
  let  transIRD;
  transIRD = (1.0 - exp(- rateI * dt))* INFD[args.nstageI-1];
  transI[args.nstageI-1] = (1-mI) * transIRD;
  transI[args.nstageI] = mI * transIRD;
  
  
  // From class H
  let  transH = array(args.nstageH + 1);
  let  rateH = args.nstageH * gammaH;
  for (i = 0; i < args.nstageH-1; i++) {
    transH[i] =  (1.0 - exp(- rateH * dt))* HOSP[i];
  }
  
  let  transHRC;
  transHRC = (1.0 - exp(- rateH * dt))* HOSP[args.nstageH-1];
  transH[args.nstageH-1] = (1-qH) * transHRC;
  transH[args.nstageH] = qH * transHRC;
  
  
  // From class C
  let  transC = array(args.nstageC + 2);
  let  rateC = args.nstageC * gammaC;
  for (i = 0; i < args.nstageC-1; i++) {
    transC[i] =  (1.0 - exp(- rateC * dt))* CARE[i];
  }
  let  transCRVM;
  transCRVM = (1.0 - exp(- rateC * dt))* CARE[args.nstageC-1];
  transC[args.nstageC-1] =  (1-mC) * (1-qC) * transCRVM;
  transC[args.nstageC] = qC * transCRVM;
  transC[args.nstageC+1] = mC * (1-qC) * transCRVM;
  
  // From class V
  let  transV = array(args.nstageV + 1);
  let  rateV = args.nstageV * gammaV;
  for (i = 0; i < args.nstageV-1; i++) {
    transV[i] =  (1.0 - exp(- rateV * dt))* VENT[i];
  }
  let  transVRD;
  transVRD = (1.0 - exp(- rateV * dt))* VENT[args.nstageV-1];
  transV[args.nstageV-1] = (1-mV) * transVRD;
  transV[args.nstageV] = mV * transVRD;
  
  // Balance the equations
  DSUSC[0] = SUSC[0];
  for (i = 0; i < args.nstageE; i++) DEXPD[i] = EXPD[i];
  for (i = 0; i < args.nstageP; i++) DPRE[i] = PRE[i];
  for (i = 0; i < args.nstageI; i++) DINFD[i] = INFD[i];
  for (i = 0; i < args.nstageH; i++) DHOSP[i] = HOSP[i];
  for (i = 0; i < args.nstageC; i++) DCARE[i] = CARE[i];
  for (i = 0; i < args.nstageV; i++) DVENT[i] = VENT[i];
  DDEAD[0] = DEAD[0];
  DRCVD[0] = RCVD[0];
  DcasesI = casesI;
  DcasesIQ = casesIQ;
  DcasesH = casesH;
  DdeathsIIQ = deathsIIQ;
  DdeathsCV = deathsCV;
  
  for (i = 0; i < args.nstageE; i++) DEXPDQ[i] = EXPDQ[i];
  for (i = 0; i < args.nstageP; i++) DPREQ[i] = PREQ[i];
  for (i = 0; i < args.nstageI; i++) DINFDQ[i] = INFDQ[i];
  
  
  
  // Balance the equations
  DSUSC[0] -= transS[0];
  DEXPD[0] += transS[0];
  for (i = 0; i < args.nstageE; i++) DEXPD[i] -= transE[i];
  for (i = 1; i < args.nstageE; i++) DEXPD[i] += transE[i-1];
  DPRE[0] += transE[args.nstageE-1];
  for (i = 0; i < args.nstageP; i++) DPRE[i] -= transP[i];
  for (i = 1; i < args.nstageP; i++) DPRE[i] += transP[i-1];
  DINFD[0] += transP[args.nstageP-1];
  DHOSP[0] += transP[args.nstageP];
  DPRE[args.nstageP-1] -= transP[args.nstageP];
  for (i = 0; i < args.nstageI; i++) DINFD[i] -= transI[i];
  for (i = 1; i < args.nstageI; i++) DINFD[i] += transI[i-1];
  for (i = 0; i < args.nstageH; i++) DHOSP[i] -= transH[i];
  for (i = 1; i < args.nstageH; i++) DHOSP[i] += transH[i-1];
  DCARE[0] += transH[args.nstageH];
  DHOSP[args.nstageH-1] -= transH[args.nstageH];
  DINFD[args.nstageI-1] -= transI[args.nstageI];
  for (i = 0; i < args.nstageC; i++) DCARE[i] -= transC[i];
  for (i = 1; i < args.nstageC; i++) DCARE[i] += transC[i-1];
  DVENT[0] += transC[args.nstageC];
  DCARE[args.nstageC-1] -= transC[args.nstageC] + transC[args.nstageC+1];
  for (i = 0; i < args.nstageV; i++) DVENT[i] -= transV[i];
  for (i = 1; i < args.nstageV; i++) DVENT[i] += transV[i-1];
  DVENT[args.nstageV-1] -= transV[args.nstageV];               
  DRCVD[0] += transI[args.nstageI-1] + transH[args.nstageH-1] + transC[args.nstageC-1] + transV[args.nstageV-1];
  DDEAD[0] += transI[args.nstageI] + transC[args.nstageC+1] + transV[args.nstageV];
  
  DSUSC[0] -= transS[1];
  DEXPDQ[0] += transS[1];
  for (i = 0; i < args.nstageE; i++) DEXPDQ[i] -= transEQ[i];
  for (i = 1; i < args.nstageE; i++) DEXPDQ[i] += transEQ[i-1];
  DPREQ[0] += transEQ[args.nstageE-1];
  for (i = 0; i < args.nstageP; i++) DPREQ[i] -= transPQ[i];
  for (i = 1; i < args.nstageP; i++) DPREQ[i] += transPQ[i-1];
  DPREQ[args.nstageP-1] -= transPQ[args.nstageP];  
  DPRE[args.nstageP-1] -= transP[args.nstageP+1];
  DINFDQ[0] += transPQ[args.nstageP-1] + transP[args.nstageP+1];
  DHOSP[0] += transPQ[args.nstageP];
  for (i = 0; i < args.nstageI; i++) DINFDQ[i] -= transIQ[i];
  for (i = 1; i < args.nstageI; i++) DINFDQ[i] += transIQ[i-1];
  DINFDQ[args.nstageI-1] -= transIQ[args.nstageI];
  DRCVD[0] += transIQ[args.nstageI-1];
  DDEAD[0] += transIQ[args.nstageI];
  
  DcasesI += transP[args.nstageP-1];
  DcasesIQ += transPQ[args.nstageP-1] + transP[args.nstageP+1];
  DcasesH += transPQ[args.nstageP] + transP[args.nstageP];
  DdeathsIIQ += transI[args.nstageI] + transIQ[args.nstageI];
  DdeathsCV += transC[args.nstageC+1] + transV[args.nstageV];
}

snippet.rInit = function(states, params, covar, t, dt, args){
  let SUSC = params.S;
  let EXPD = params.E1;
  let PRE = params.P1;
  let INFD = params.I1;
  let HOSP = params.H1;
  let CARE = params.C1;
  let VENT = params.V1;
  let DEAD = params.M;
  let RCVD = params.R;
  
  let EXPDQ = params.EQ1;
  let PREQ = params.PQ1;
  let INFDQ = params.IQ1;
  
  let i;
  
  let fS, fEQ, fPQ, fIQ, fE, fP, fI, fH, fC, fR, fM, fV; 
  fS = S0;
  fEQ = EQ0/args.nstageE;
  fPQ = PQ0/args.nstageP;
  fIQ = IQ0/args.nstageI;
  fE = E0/args.nstageE;
  fP = P0/args.nstageP;
  fI = I0/args.nstageI;
  fH = H0/args.nstageH;
  fC = C0/args.nstageC;
  fV = V0/args.nstageV;
  fM = M0;
  fR = 1 - fS - args.nstageE*fE - args.nstageP*fP - args.nstageI*fI - args.nstageH*fH - args.nstageC*fC - args.nstageV*fV - args.nstageE*fEQ - args.nstageP*fPQ - args.nstageI*fIQ - fM;
  
  SUSC[0] = 0;//////TODO: nearbyint(pop*fS);
  for (i = 0; i < args.nstageE; i++) EXPDQ[i] = 0;//////// nearbyint(pop*fEQ);
  for (i = 0; i < args.nstageP; i++) PREQ[i] = 0;//////// nearbyint(pop*fPQ);
  for (i = 0; i < args.nstageI; i++) INFDQ[i] = 0;//////// nearbyint(pop*fIQ);
  for (i = 0; i < args.nstageE; i++) EXPD[i] = 0;//////// nearbyint(pop*fE);
  for (i = 0; i < args.nstageP; i++) PRE[i] = 0;//////// nearbyint(pop*fP);
  for (i = 0; i < args.nstageI; i++) INFD[i] = 0;//////// nearbyint(pop*fI);
  for (i = 0; i < args.nstageH; i++) HOSP[i] = 0;//////// nearbyint(pop*fH);
  for (i = 0; i < args.nstageC; i++) CARE[i] = 0;//////// nearbyint(pop*fC);
  for (i = 0; i < args.nstageV; i++) VENT[i] = 0;//////// nearbyint(pop*fV);
  DEAD[0] = 0;//////// nearbyint(pop*fM);
  RCVD[0] = 0;//////// nearbyint(pop*fR);
  
  casesI = 0.0;
  casesIQ = 0.0;
  casesH = 0.0;
  deathsIIQ = 0.0;
  deathsCV = 0.0;
}




snippet.toEstimationScale = function(params) {
  let esimParams = Object.assign({}, params);
  esimParams.mu = Math.log(params.mu);
  esimParams.psi = Math.log(params.psi);
  esimParams.sigma = Math.log(params.sigma);
  esimParams.gamma = Math.log(params.gamma);
  esimParams.R0 = Math.log(params.R0);
  esimParams.rho = mathLib.logit(params.rho);
  esimParams.amplitude = mathLib.logit(params.amplitude);
  
  return esimParams;
}

snippet.fromEstimationScale = function(params) {
  let esimParams = Object.assign({}, params);
  esimParams.mu = Math.exp(params.mu);
  esimParams.psi = Math.exp(params.psi);
  esimParams.sigma = Math.exp(params.sigma);
  esimParams.gamma = Math.exp(params.gamma);
  esimParams.R0 = Math.exp(params.R0);
  esimParams.rho = mathLib.expit(params.rho);
  esimParams.amplitude = mathLib.expit(params.amplitude);
  
  return esimParams
}






snippet.statenames = function (args = { nstageE: 3, nstageI: 3, nstageP: 3, nstageH: 3, nstageC: 3, nstageV: 3 }) {
  return [ "S",
    ...[...Array(args.nstageE)].map((_, i) => 'EQ' + (i + 1)),
    ...[...Array(args.nstageP)].map((_, i) => 'PQ' + (i + 1)),
    ...[...Array(args.nstageI)].map((_, i) => 'IQ' + (i + 1)),
    ...[...Array(args.nstageE)].map((_, i) => 'E' + (i + 1)),
    ...[...Array(args.nstageP)].map((_, i) => 'P' + (i + 1)),
    ...[...Array(args.nstageI)].map((_, i) => 'I' + (i + 1)),
    ...[...Array(args.nstageH)].map((_, i) => 'H' + (i + 1)),
    ...[...Array(args.nstageC)].map((_, i) => 'C' + (i + 1)),
    ...[...Array(args.nstageV)].map((_, i) => 'V' + (i + 1)),
    // paste0("EQ",seq_len(nstageE)),
    // paste0("PQ",seq_len(nstageP)),
    // paste0("IQ",seq_len(nstageI)),
    // paste0("E",seq_len(nstageE)),
    // paste0("P",seq_len(nstageP)),
    // paste0("I",seq_len(nstageI)),
    // paste0("H",seq_len(nstageH)),
    // paste0("C",seq_len(nstageC)),
    // paste0("V",seq_len(nstageV)),
    "M", "R"
  ];                
}

snippet.icnames = function() { return ["S0","EQ0", "PQ0", "IQ0", "E0", "P0", "I0", "H0", "C0", "V0", "M0"]; };

snippet.zeronames = ["casesI", "casesIQ", "casesH", "deathsIIQ", "deathsCV"];

snippet.params_log = ["betaI", "iota","beta_sd", "sigma", "kappa", "gammaI", "gammaH", "gammaC", "gammaV","TF"];

snippet.params_logit = ["rho", "theta","dI0", "dP0", "dT0", "dB0","dI1", "dP1", "dT1", "dB1","qP", "qH", "qC", "mI", "mC", "mV"];

snippet.params_mod = [ ...snippet.params_log, ...snippet.params_logit];

snippet.params_ic = snippet.icnames();

snippet.dObs();
module.exports = snippet