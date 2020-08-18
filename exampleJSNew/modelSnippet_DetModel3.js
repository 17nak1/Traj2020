let snippet = {}
let mathLib = require('./mathLib');

let nstageE = 3;
let nstageP = 3;
let nstageI = 3;
let nstageH = 3;
let nstageC = 3;
let nstageV = 3;
snippet.endTime = "2020-07-16"
snippet.T0 = 75
snippet.T1 = 139

snippet.skeleton = function (states, params, t, dt, covar, args) {
  if(t==0){
    console.log(params.betaI, params.iota);
  }
  let d ={};
  dt = 0.1;
  let SUSC = [states.S];
  let DEAD = [states.M];
  let RCVD = [states.R];

  let EXPD = [], PRE = [], INFD = [], HOSP = [], CARE = [], VENT = [];
  for(let i = 0; i < args.nstageE; i++) EXPD.push(states[`E${i + 1}`]);
  for(let i = 0; i < args.nstageP; i++) PRE.push(states[`P${i + 1}`]);
  for(let i = 0; i < args.nstageI; i++) INFD.push(states[`I${i + 1}`]);
  for(let i = 0; i < args.nstageH; i++) HOSP.push(states[`H${i + 1}`]);
  for(let i = 0; i < args.nstageC; i++) CARE.push(states[`C${i + 1}`]);
  for(let i = 0; i < args.nstageV; i++) VENT.push(states[`V${i + 1}`]);
  
  let EXPDQ = [], PREQ = [], INFDQ = [];
  for(let i = 0; i < args.nstageE; i++) EXPDQ.push(states[`EQ${i + 1}`]);
  for(let i = 0; i < args.nstageP; i++) PREQ.push(states[`PQ${i + 1}`]);
  for(let i = 0; i < args.nstageI; i++) INFDQ.push(states[`IQ${i + 1}`]);

  // let SUSC = &S;
  // let EXPD = &E1;
  // let PRE = &P1;
  // let INFD = &I1;
  // let HOSP = &H1;
  // let CARE = &C1;
  // let VENT = &V1;
  // let DEAD = &M;
  // let RCVD = &R;
  
  let DSUSC = new Array(1) //&DS;
  let DEXPD = new Array(args.nstageE) //&DE1;
  let DPRE = new Array(args.nstageP) //&DP1;
  let DINFD = new Array(args.nstageI) //&DI1;
  let DHOSP = new Array(args.nstageH) //&DH1;
  let DCARE = new Array(args.nstageC) //&DC1;
  let DVENT = new Array(args.nstageV) //&DV1;
  let DDEAD = new Array(1) //&DM;
  let DRCVD = new Array(1) //&DR;
  
  
  // let EXPDQ = &EQ1;
  // let PREQ = &PQ1;
  // let INFDQ = &IQ1;
  let DEXPDQ = new Array(args.nstageE)//&DEQ1;
  let DPREQ = new Array(args.nstageP)//&DPQ1;
  let DINFDQ =new Array(args.nstageI)// &DIQ1;
  
  
  // Different transmission rates
  
  let dQdt;
  
  let TOT_PRE = 0;
  for(let i = 0; i < args.nstageP; i++) {
    TOT_PRE += PRE[i];
  }

  let TOT_INFD = 0;
  for(i = 0; i < args.nstageI; i++) {
  TOT_INFD += INFD[i];
  }
  
  let PD, lambdaI, lambdaP, lambdaPQ, lambda, lambdaQ;
  if ( isFinite(covar.tests) ) {
    PD = params.rho * (covar.tests / (covar.tests + params.TF));
  } else {
    PD = params.rho * (15.0 / (15.0 + params.TF)) ;
  }
  lambdaI = params.betaI * TOT_INFD;
  lambdaP = params.betaI * params.theta * TOT_PRE * (1-PD);
  lambdaPQ = params.betaI * params.theta * TOT_PRE * PD;
  
  if (t < args.T0) {
    dQdt  = mathLib.rgammawn(params.beta_sd, dt)/dt;
    lambda = ( (lambdaI + lambdaP + params.iota) / args.pop ) * dQdt;
    lambdaQ = ( (lambdaPQ) / args.pop ) * dQdt;
  } else if (t < args.T0 + 7.0) {
    let x = (t-args.T0) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = mathLib.rgammawn(((1-ss) + ss*params.dB0)*params.beta_sd, dt)/dt;
    lambda = ( ( ((1-ss) + ss*params.dI0) * lambdaI + ((1-ss) + ss*params.dP0) * lambdaP + ((1-ss) + ss*params.dT0)*params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( ( ((1-ss) + ss*params.dP0) * lambdaPQ  ) / args.pop ) * dQdt;
  } else if (t < args.T1) {
    dQdt  = mathLib.rgammawn(params.dB0*params.beta_sd, dt)/dt;
    lambda = ( ( params.dI0 * lambdaI + params.dP0 * lambdaP + params.dT0 * params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( (  params.dP0 * lambdaPQ  ) / args.pop ) * dQdt;
  } else if (t < args.T1 +7.0) {
    let x = (t-args.T1) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = mathLib.rgammawn(((1-ss)*params.dB0 + ss*params.dB1)*params.beta_sd, dt)/dt;
    lambda = ( ( ((1-ss)*params.dI0 + ss*params.dI1) * lambdaI + ((1-ss)*params.dP0 + ss*params.dP1) * lambdaP + ((1-ss)*params.dT0 + ss*params.dT1)*params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( ( ((1-ss)*params.dP0 + ss*params.dP1) * lambdaPQ  ) / args.pop ) * dQdt;
  } else {
    dQdt  = mathLib.rgammawn(params.dB1*params.beta_sd, dt)/dt;
    lambda = ( ( params.dI1 * lambdaI + params.dP1 * lambdaP + params.dT1 * params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( (  params.dP1 * lambdaPQ  ) / args.pop ) * dQdt;
  }
  
  // From class S
  let transS = new Array(2);
  let rateS = new Array(2);
  rateS[0] = lambda ;
  rateS[1] = lambdaQ;
  let rateS_tot = rateS[0] + rateS[1];
  let transS_tot = (1 - Math.exp(- rateS_tot * dt)) * SUSC[0];
  transS[0] = rateS[0] / rateS_tot * transS_tot;
  transS[1] = rateS[1] / rateS_tot * transS_tot;
  
  // From class EQ
  let transEQ = new Array(args.nstageE);
  let rateEQ = args.nstageE * params.sigma;
  for (let i = 0; i < args.nstageE; i++) {
    transEQ[i] =  (1 - Math.exp(- rateEQ * dt))* EXPDQ[i];
  }
  
  // From class PQ
  let transPQ = new Array(args.nstageP + 1);
  let ratePQ = args.nstageP * params.kappa;
  for (let i = 0; i < args.nstageP-1; i++) {
    transPQ[i] =  (1 - Math.exp(- ratePQ * dt))* PREQ[i];
  }
  let transPQIQH;
  transPQIQH = (1 - Math.exp(- ratePQ * dt))* PREQ[args.nstageP-1];
  transPQ[args.nstageP-1] = (1-params.qP) * transPQIQH;
  transPQ[args.nstageP] = params.qP * transPQIQH;
  
  // From class IQ
  let transIQ = new Array(args.nstageI + 1);
  let rateIQ = args.nstageI * params.gammaI;
  for (let i = 0; i < args.nstageI-1; i++) {
    transIQ[i] =  (1 - Math.exp(- rateIQ * dt))* INFDQ[i];
  }
  let transIQRD;
  transIQRD = (1 - Math.exp(- rateIQ * dt))* INFDQ[args.nstageI-1];
  transIQ[args.nstageI-1] = (1 - params.mI) * transIQRD;
  transIQ[args.nstageI] = params.mI * transIQRD;
  
  
  // From class E
  let transE = new Array(args.nstageE);
  let rateE = args.nstageE * params.sigma;
  for (let i = 0; i < args.nstageE; i++) {
    transE[i] =  (1 - Math.exp(- rateE * dt))* EXPD[i];
  }
  
  // From class P
  let transP = new Array(args.nstageP + 2);
  let rateP = args.nstageP * params.kappa;
  for (let i = 0; i < args.nstageP-1; i++) {
    transP[i] =  (1 - Math.exp(- rateP * dt))* PRE[i];
  }
  let transPIHIQ;
  transPIHIQ = (1.0 - Math.exp(- rateP * dt))* PRE[args.nstageP-1];
  transP[args.nstageP-1] = (1 - PD) * (1 - params.qP) * transPIHIQ;
  transP[args.nstageP] = params.qP * transPIHIQ;
  transP[args.nstageP+1] = PD * (1 - params.qP) * transPIHIQ;
  
  // From class I
  let transI = new Array(args.nstageI + 1);
  let rateI = args.nstageI * params.gammaI;
  for (let i = 0; i < args.nstageI-1; i++) {
    transI[i] =  (1 - Math.exp(- rateI * dt))* INFD[i];
  }
  
  let transIRD;
  transIRD = (1 - Math.exp(- rateI * dt))* INFD[args.nstageI-1];
  transI[args.nstageI-1] = (1 - params.mI) * transIRD;
  transI[args.nstageI] = params.mI * transIRD;
  
  
  // From class H
  let transH = new Array(args.nstageH + 1);
  let rateH = args.nstageH * params.gammaH;
  for (let i = 0; i < args.nstageH-1; i++) {
    transH[i] =  (1.0 - Math.exp(- rateH * dt))* HOSP[i];
  }
  
  let transHRC;
  transHRC = (1.0 - Math.exp(- rateH * dt))* HOSP[args.nstageH-1];
  transH[args.nstageH-1] = (1 - params.qH) * transHRC;
  transH[args.nstageH] = params.qH * transHRC;
  
  
  // From class C
  let transC= new Array(args.nstageC + 2);
  let rateC = args.nstageC * params.gammaC;
  for (let i = 0; i < args.nstageC-1; i++) {
    transC[i] =  (1.0 - Math.exp(- rateC * dt))* CARE[i];
  }
  let transCRVM;
  transCRVM = (1.0 - Math.exp(- rateC * dt))* CARE[args.nstageC-1];
  transC[args.nstageC-1] =  (1 - params.mC) * (1 - params.qC) * transCRVM;
  transC[args.nstageC] = params.qC * transCRVM;
  transC[args.nstageC+1] = params.mC * (1 - params.qC) * transCRVM;
  
  // From class V
  let transV = new Array(args.nstageV + 1);
  let rateV = args.nstageV * params.gammaV;
  for (let i = 0; i < args.nstageV-1; i++) {
    transV[i] =  (1 - Math.exp(- rateV * dt))* VENT[i];
  }
  let transVRD;
  transVRD = (1.0 - Math.exp(- rateV * dt))* VENT[args.nstageV-1];
  transV[args.nstageV-1] = (1 - params.mV) * transVRD;
  transV[args.nstageV] = params.mV * transVRD;
  
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
  d.casesI = states.casesI;
  d.casesIQ = states.casesIQ;
  d.casesH = states.casesH;
  d.deathsIIQ = states.deathsIIQ;
  d.deathsCV = states.deathsCV;
  
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
  
  d.casesI += transP[args.nstageP-1];
  d.casesIQ += transPQ[args.nstageP-1] + transP[args.nstageP+1];
  d.casesH += transPQ[args.nstageP] + transP[args.nstageP];
  d.deathsIIQ += transI[args.nstageI] + transIQ[args.nstageI];
  d.deathsCV += transC[args.nstageC+1] + transV[args.nstageV];

  // TODO: result should these parameters
  // let DSUSC = new Array(1) //&DS;
  d['S'] = DSUSC[0];
  // let DEXPD = new Array(args.nstageE) //&DE1;
  d['E1'] = DEXPD[0];
  d['E2'] = DEXPD[1];
  d['E3'] = DEXPD[2];
  // let DPRE = new Array(args.nstageP) //&DP1;
  d['P1'] = DPRE[0];
  d['P2'] = DPRE[1];
  d['P3'] = DPRE[2];
  // let DINFD = new Array(args.nstageI) //&DI1;
  d['I1'] = DINFD[0];
  d['I2'] = DINFD[1];
  d['I3'] = DINFD[2];
  // let DHOSP = new Array(args.nstageH) //&DH1;
  d['H1'] = DHOSP[0];
  d['H2'] = DHOSP[1];
  d['H3'] = DHOSP[2];
  // let DCARE = new Array(args.nstageC) //&DC1;
  d['C1'] = DCARE[0];
  d['C2'] = DCARE[1];
  d['C3'] = DCARE[2];
  // let DVENT = new Array(args.nstageV) //&DV1;
  d['V1'] = DVENT[0];
  d['V2'] = DVENT[1];
  d['V3'] = DVENT[2];
  // let DDEAD = new Array(1) //&DM;
  d['M'] = DDEAD[0];
  // let DRCVD = new Array(1) //&DR;
  d['R'] = DRCVD[0];
  // let DEXPDQ = new Array(args.nstageE)//&DEQ1;
  d['EQ1'] = DEXPDQ[0];
  d['EQ2'] = DEXPDQ[1];
  d['EQ3'] = DEXPDQ[2];
  // let DPREQ = new Array(args.nstageP)//&DPQ1;
  d['PQ1'] = DPREQ[0];
  d['PQ2'] = DPREQ[1];
  d['PQ3'] = DPREQ[2];
  // let DINFDQ =new Array(args.nstageI)// &DIQ1;
  d['IQ1'] = DINFDQ[0];
  d['IQ2'] = DINFDQ[1];
  d['IQ3'] = DINFDQ[2];

  return d;
}
// Note: for translating "double *SUSC = &S" we directly fill S and no need to define SUSC
snippet.initializer = function(params, covar, args) {
  let initObj = {};
  
  let fS = params.S0;
  let fEQ = params.EQ0 / args.nstageE;
  let fPQ = params.PQ0 / args.nstageP;
  let fIQ = params.IQ0 / args.nstageI;
  let fE = params.E0 / args.nstageE;
  let fP = params.P0 / args.nstageP;
  let fI = params.I0 / args.nstageI;
  let fH = params.H0 / args.nstageH;
  let fC = params.C0 / args.nstageC;
  let fV = params.V0 / args.nstageV;
  let fM = params.M0;
  let fR = 1 - fS - nstageE*fE - nstageP*fP - nstageI*fI - nstageH*fH - nstageC*fC - nstageV*fV - nstageE*fEQ - nstageP*fPQ - nstageI*fIQ - fM;
  
  initObj["S"] = Math.round(args.pop*fS);
  for (let i = 0; i < nstageE; i++) initObj[`EQ${i + 1}`] = Math.round(args.pop*fEQ);
  for (let i = 0; i < nstageP; i++) initObj[`PQ${i + 1}`] = Math.round(args.pop*fPQ);
  for (let i = 0; i < nstageI; i++) initObj[`IQ${i + 1}`] = Math.round(args.pop*fIQ);
  for (let i = 0; i < nstageE; i++) initObj[`E${i + 1}`] = Math.round(args.pop*fE);
  for (let i = 0; i < nstageP; i++) initObj[`P${i + 1}`] = Math.round(args.pop*fP);
  for (let i = 0; i < nstageI; i++) initObj[`I${i + 1}`] = Math.round(args.pop*fI);
  for (let i = 0; i < nstageH; i++) initObj[`H${i + 1}`] = Math.round(args.pop*fH);
  for (let i = 0; i < nstageC; i++) initObj[`C${i + 1}`] = Math.round(args.pop*fC);
  for (let i = 0; i < nstageV; i++) initObj[`V${i + 1}`] = Math.round(args.pop*fV);
  initObj["M"] = Math.round(args.pop*fM);
  initObj["R"] = Math.round(args.pop*fR);
  
  initObj["casesI"] = 0;
  initObj["casesIQ"] = 0;
  initObj["casesH"] = 0;
  initObj["deathsIIQ"] = 0;
  initObj["deathsCV"] = 0;
  
  return initObj;
}

snippet.dmeasure = function (data ,states, params, giveLog = 1) {
  let HOSP = [], CARE = [], VENT = [];
  for(let i = 0; i < nstageH; i++) HOSP.push(states[`H${i + 1}`]);
  for(let i = 0; i < nstageC; i++) CARE.push(states[`C${i + 1}`]);
  for(let i = 0; i < nstageV; i++) VENT.push(states[`V${i + 1}`]);
  
  let lik_reports, lik_deaths, lik_hospital, lik_ICU, lik_ventilator;
  
  if (isFinite(data.reports)) {
    let reported_cases = states.casesH + states.casesIQ;
    lik_reports = mathLib.dpois(data.reports, reported_cases + 1e-6, 1);
  } else {
    lik_reports = 0;
  }
  
  if (isFinite(data.deaths)) {
    let reported_deaths = states.deathsCV + states.deathsIIQ;
    lik_deaths = mathLib.dpois(data.deaths, reported_deaths + 1e-6, 1); 
  } else {
    lik_deaths = 0;
  }
  
  let TOT_H = 0;
  for(let i = 0; i < nstageH; i++) {
    TOT_H += HOSP[i];
  }

  let TOT_C = 0;
  for(let i = 0; i < nstageC; i++) {
    TOT_C += CARE[i];
  }

  let TOT_V = 0;
  for(let i = 0; i < nstageV; i++) {
    TOT_V += VENT[i];
  }

  if (isFinite(data.hospital) & isFinite(data.ICU)) {
    lik_hospital = mathLib.dpois(data.hospital - data.ICU, TOT_H + 1e-6, 1);
  } else {
    lik_hospital = 0;
  }
  if (isFinite(data.ICU) & isFinite(data.ventilator)) {
    lik_ICU = mathLib.dpois(data.ICU - data.ventilator, TOT_C + 1e-6, 1);
  } else {
    lik_ICU = 0;
  }
  if (isFinite(data.ventilator)) {
    lik_ventilator = mathLib.dpois(data.ventilator, TOT_V + 1e-6, 1);
  } else {
    lik_ventilator = 0;
  }
  
  let lik = lik_reports + lik_deaths + lik_hospital + lik_ICU + lik_ventilator;
  
  if (giveLog == 0 ) {
    lik = Math.exp(lik);
  }
  
  return lik;
}


let params_log = ["betaI", "iota","beta_sd", "sigma", "kappa", "gammaI", "gammaH", "gammaC", "gammaV"];
let params_logit = ["rho", "theta", "dI0", "dP0", "dT0", "dB0","dI1", "dP1", "dT1", "dB1","qP", "qH", "qC", "mI", "mC", "mV"];
let statenamesFn = function() {
  let sn = ["S"]
  for (let i = 0; i < nstageE; i++) sn.push("EQ"+i);
  for (let i = 0; i < nstageP; i++) sn.push("PQ"+i);
  for (let i = 0; i < nstageI; i++) sn.push("IQ"+i);
  for (let i = 0; i < nstageE; i++) sn.push("E"+i);
  for (let i = 0; i < nstageP; i++) sn.push("P"+i);
  for (let i = 0; i < nstageI; i++) sn.push("I"+i);
  for (let i = 0; i < nstageH; i++) sn.push("H"+i);
  for (let i = 0; i < nstageC; i++) sn.push("C"+i);
  for (let i = 0; i < nstageV; i++) sn.push("V"+i);
  sn.push("M", "R");
  return sn;
}    

snippet.paramsMod = [...params_log, ...params_logit];
snippet.paramsIc = ["S0","EQ0", "PQ0", "IQ0", "E0", "P0", "I0", "H0", "C0", "V0", "M0"];
snippet.zeronames = ["casesI", "casesIQ", "casesH", "deathsIIQ", "deathsCV"];
snippet.statenames = [...statenamesFn(), ...snippet.zeronames];

snippet.toEstimationScale = function(params) {
  let estiParams = Object.assign({}, params);
  for (let i = 0; i < params_log.length; i++) {
    estiParams[params_log[i]] = Math.log(params[params_log[i]]);
  }

  for (let i = 0; i < params_logit.length; i++) {
    estiParams[params_logit[i]] = mathLib.qlogis(params[params_logit[i]]);
  }
  
  return estiParams;
}

snippet.fromEstimationScale = function(params) {
  let estiParams = Object.assign({}, params);
  for (let i = 0; i < params_log.length; i++) {
    estiParams[params_log[i]] = Math.exp(params[params_log[i]]);
  }

  for (let i = 0; i < params_logit.length; i++) {
    estiParams[params_logit[i]] = mathLib.plogis(params[params_logit[i]]);
  }
  
  return estiParams;
}



module.exports = snippet

