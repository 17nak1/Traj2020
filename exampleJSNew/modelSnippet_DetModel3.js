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

snippet.skeleton = function (states, params, t, covar) {
  let d ={};
  let dt = 0.1;
  let SUSC = [states.S];
  let DEAD = [states.M];
  let RCVD = [states.R];

  let EXPD = [], PRE = [], INFD = [], HOSP = [], CARE = [], VENT = [];
  for(let i = 0; i < nstageE; i++) EXPD.push(states[`E${i + 1}`]);
  for(let i = 0; i < nstageP; i++) PRE.push(states[`P${i + 1}`]);
  for(let i = 0; i < nstageI; i++) INFD.push(states[`I${i + 1}`]);
  for(let i = 0; i < nstageH; i++) HOSP.push(states[`H${i + 1}`]);
  for(let i = 0; i < nstageC; i++) CARE.push(states[`C${i + 1}`]);
  for(let i = 0; i < nstageV; i++) VENT.push(states[`V${i + 1}`]);
  
  let EXPDQ = [], PREQ = [], INFDQ = [];
  for(let i = 0; i < nstageE; i++) EXPDQ.push(states[`EQ${i + 1}`]);
  for(let i = 0; i < nstageP; i++) PREQ.push(states[`PQ${i + 1}`]);
  for(let i = 0; i < nstageI; i++) INFDQ.push(states[`IQ${i + 1}`]);

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
  let DEXPD = new Array(nstageE) //&DE1;
  let DPRE = new Array(nstageP) //&DP1;
  let DINFD = new Array(nstageI) //&DI1;
  let DHOSP = new Array(nstageH) //&DH1;
  let DCARE = new Array(nstageC) //&DC1;
  let DVENT = new Array(nstageV) //&DV1;
  let DDEAD = new Array(1) //&DM;
  let DRCVD = new Array(1) //&DR;
  
  
  // let EXPDQ = &EQ1;
  // let PREQ = &PQ1;
  // let INFDQ = &IQ1;
  let DEXPDQ = new Array(nstageE)//&DEQ1;
  let DPREQ = new Array(nstageP)//&DPQ1;
  let DINFDQ =new Array(nstageI)// &DIQ1;
  
  
  // Different transmission rates
  
  let dQdt;
  
  let TOT_PRE = 0;
  for(let i = 0; i < nstageP; i++) {
    TOT_PRE += PRE[i];
  }

  let TOT_INFD = 0;
  for(i = 0; i < nstageI; i++) {
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
  
  if (t < T0) {
    dQdt  = mathLib.rgammawn(params.beta_sd, dt)/dt;
    lambda = ( (lambdaI + lambdaP + params.iota) / pop ) * dQdt;
    lambdaQ = ( (lambdaPQ) / pop ) * dQdt;
  } else if (t < T0 + 7.0) {
    let x = (t-T0) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = mathLib.rgammawn(((1-ss) + ss*dB0)*params.beta_sd, dt)/dt;
    lambda = ( ( ((1-ss) + ss*dI0) * lambdaI + ((1-ss) + ss*dP0) * lambdaP + ((1-ss) + ss*dT0)*params.iota ) / pop ) * dQdt;
    lambdaQ = ( ( ((1-ss) + ss*dP0) * lambdaPQ  ) / pop ) * dQdt;
  } else if (t < T1) {
    dQdt  = mathLib.rgammawn(dB0*params.beta_sd, dt)/dt;
    lambda = ( ( dI0 * lambdaI + dP0 * lambdaP + dT0 * params.iota ) / pop ) * dQdt;
    lambdaQ = ( (  dP0 * lambdaPQ  ) / pop ) * dQdt;
  } else if (t < T1 +7.0) {
    let x = (t-T1) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = mathLib.rgammawn(((1-ss)*dB0 + ss*dB1)*params.beta_sd, dt)/dt;
    lambda = ( ( ((1-ss)*dI0 + ss*dI1) * lambdaI + ((1-ss)*dP0 + ss*dP1) * lambdaP + ((1-ss)*dT0 + ss*dT1)*params.iota ) / pop ) * dQdt;
    lambdaQ = ( ( ((1-ss)*dP0 + ss*dP1) * lambdaPQ  ) / pop ) * dQdt;
  } else {
    dQdt  = mathLib.rgammawn(dB1*params.beta_sd, dt)/dt;
    lambda = ( ( dI1 * lambdaI + dP1 * lambdaP + dT1 * params.iota ) / pop ) * dQdt;
    lambdaQ = ( (  dP1 * lambdaPQ  ) / pop ) * dQdt;
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
  let transEQ = new Array(nstageE);
  let rateEQ = nstageE * sigma;
  for (let i = 0; i < nstageE; i++) {
    transEQ[i] =  (1 - Math.exp(- rateEQ * dt))* EXPDQ[i];
  }
  
  // From class PQ
  let transPQ = new Array(nstageP + 1);
  let ratePQ = nstageP * kappa;
  for (let i = 0; i < nstageP-1; i++) {
    transPQ[i] =  (1 - Math.exp(- ratePQ * dt))* PREQ[i];
  }
  let transPQIQH;
  transPQIQH = (1 - Math.exp(- ratePQ * dt))* PREQ[nstageP-1];
  transPQ[nstageP-1] = (1-params.qP) * transPQIQH;
  transPQ[nstageP] = params.qP * transPQIQH;
  
  // From class IQ
  let transIQ = new Array(nstageI + 1);
  let rateIQ = nstageI * gammaI;
  for (let i = 0; i < nstageI-1; i++) {
    transIQ[i] =  (1 - Math.exp(- rateIQ * dt))* INFDQ[i];
  }
  let transIQRD;
  transIQRD = (1 - Math.exp(- rateIQ * dt))* INFDQ[nstageI-1];
  transIQ[nstageI-1] = (1 - params.mI) * transIQRD;
  transIQ[nstageI] = params.mI * transIQRD;
  
  
  // From class E
  let transE = new Array(nstageE);
  let rateE = nstageE * sigma;
  for (let i = 0; i < nstageE; i++) {
    transE[i] =  (1 - Math.exp(- rateE * dt))* Math.EXPD[i];
  }
  
  // From class P
  let transP = new Array(nstageP + 2);
  let rateP = nstageP * kappa;
  for (let i = 0; i < nstageP-1; i++) {
    transP[i] =  (1 - Math.exp(- rateP * dt))* PRE[i];
  }
  let transPIHIQ;
  transPIHIQ = (1.0 - Math.exp(- rateP * dt))* PRE[nstageP-1];
  transP[nstageP-1] = (1 - PD) * (1 - params.qP) * transPIHIQ;
  transP[nstageP] = params.qP * transPIHIQ;
  transP[nstageP+1] = PD * (1 - params.qP) * transPIHIQ;
  
  // From class I
  let transI = new Array(stageI + 1);
  let rateI = nstageI * gammaI;
  for (let i = 0; i < nstageI-1; i++) {
    transI[i] =  (1 - Math.exp(- rateI * dt))* INFD[i];
  }
  
  let transIRD;
  transIRD = (1 - Math.exp(- rateI * dt))* INFD[nstageI-1];
  transI[nstageI-1] = (1 - Math.mI) * transIRD;
  transI[nstageI] = Math.mI * transIRD;
  
  
  // From class H
  let transH = new Array(nstageH + 1);
  let rateH = nstageH * gammaH;
  for (let i = 0; i < nstageH-1; i++) {
    transH[i] =  (1.0 - Math.exp(- rateH * dt))* HOSP[i];
  }
  
  let transHRC;
  transHRC = (1.0 - Math.exp(- rateH * dt))* HOSP[nstageH-1];
  transH[nstageH-1] = (1 - params.qH) * transHRC;
  transH[nstageH] = params.qH * transHRC;
  
  
  // From class C
  let transC= new Array(nstageC + 2);
  let rateC = nstageC * gammaC;
  for (let i = 0; i < nstageC-1; i++) {
    transC[i] =  (1.0 - Math.exp(- rateC * dt))* CARE[i];
  }
  let transCRVM;
  transCRVM = (1.0 - Math.exp(- rateC * dt))* CARE[nstageC-1];
  transC[nstageC-1] =  (1 - params.mC) * (1 - params.qC) * transCRVM;
  transC[nstageC] = params.qC * transCRVM;
  transC[nstageC+1] = params.mC * (1 - params.qC) * transCRVM;
  
  // From class V
  let transV = new Array(nstageV + 1);
  let rateV = nstageV * gammaV;
  for (let i = 0; i < nstageV-1; i++) {
    transV[i] =  (1 - Math.exp(- rateV * dt))* VENT[i];
  }
  let transVRD;
  transVRD = (1.0 - Math.exp(- rateV * dt))* VENT[nstageV-1];
  transV[nstageV-1] = (1 - params.mV) * transVRD;
  transV[nstageV] = params.mV * transVRD;
  
  // Balance the equations
  DSUSC[0] = SUSC[0];
  for (i = 0; i < nstageE; i++) DEXPD[i] = EXPD[i];
  for (i = 0; i < nstageP; i++) DPRE[i] = PRE[i];
  for (i = 0; i < nstageI; i++) DINFD[i] = INFD[i];
  for (i = 0; i < nstageH; i++) DHOSP[i] = HOSP[i];
  for (i = 0; i < nstageC; i++) DCARE[i] = CARE[i];
  for (i = 0; i < nstageV; i++) DVENT[i] = VENT[i];
  DDEAD[0] = DEAD[0];
  DRCVD[0] = RCVD[0];
  DcasesI = casesI;
  DcasesIQ = casesIQ;
  DcasesH = casesH;
  DdeathsIIQ = deathsIIQ;
  DdeathsCV = deathsCV;
  
  for (i = 0; i < nstageE; i++) DEXPDQ[i] = EXPDQ[i];
  for (i = 0; i < nstageP; i++) DPREQ[i] = PREQ[i];
  for (i = 0; i < nstageI; i++) DINFDQ[i] = INFDQ[i];
  
  
  
  // Balance the equations
  DSUSC[0] -= transS[0];
  DEXPD[0] += transS[0];
  for (i = 0; i < nstageE; i++) DEXPD[i] -= transE[i];
  for (i = 1; i < nstageE; i++) DEXPD[i] += transE[i-1];
  DPRE[0] += transE[nstageE-1];
  for (i = 0; i < nstageP; i++) DPRE[i] -= transP[i];
  for (i = 1; i < nstageP; i++) DPRE[i] += transP[i-1];
  DINFD[0] += transP[nstageP-1];
  DHOSP[0] += transP[nstageP];
  DPRE[nstageP-1] -= transP[nstageP];
  for (i = 0; i < nstageI; i++) DINFD[i] -= transI[i];
  for (i = 1; i < nstageI; i++) DINFD[i] += transI[i-1];
  for (i = 0; i < nstageH; i++) DHOSP[i] -= transH[i];
  for (i = 1; i < nstageH; i++) DHOSP[i] += transH[i-1];
  DCARE[0] += transH[nstageH];
  DHOSP[nstageH-1] -= transH[nstageH];
  DINFD[nstageI-1] -= transI[nstageI];
  for (i = 0; i < nstageC; i++) DCARE[i] -= transC[i];
  for (i = 1; i < nstageC; i++) DCARE[i] += transC[i-1];
  DVENT[0] += transC[nstageC];
  DCARE[nstageC-1] -= transC[nstageC] + transC[nstageC+1];
  for (i = 0; i < nstageV; i++) DVENT[i] -= transV[i];
  for (i = 1; i < nstageV; i++) DVENT[i] += transV[i-1];
  DVENT[nstageV-1] -= transV[nstageV];               
  DRCVD[0] += transI[nstageI-1] + transH[nstageH-1] + transC[nstageC-1] + transV[nstageV-1];
  DDEAD[0] += transI[nstageI] + transC[nstageC+1] + transV[nstageV];
  
  DSUSC[0] -= transS[1];
  DEXPDQ[0] += transS[1];
  for (i = 0; i < nstageE; i++) DEXPDQ[i] -= transEQ[i];
  for (i = 1; i < nstageE; i++) DEXPDQ[i] += transEQ[i-1];
  DPREQ[0] += transEQ[nstageE-1];
  for (i = 0; i < nstageP; i++) DPREQ[i] -= transPQ[i];
  for (i = 1; i < nstageP; i++) DPREQ[i] += transPQ[i-1];
  DPREQ[nstageP-1] -= transPQ[nstageP];  
  DPRE[nstageP-1] -= transP[nstageP+1];
  DINFDQ[0] += transPQ[nstageP-1] + transP[nstageP+1];
  DHOSP[0] += transPQ[nstageP];
  for (i = 0; i < nstageI; i++) DINFDQ[i] -= transIQ[i];
  for (i = 1; i < nstageI; i++) DINFDQ[i] += transIQ[i-1];
  DINFDQ[nstageI-1] -= transIQ[nstageI];
  DRCVD[0] += transIQ[nstageI-1];
  DDEAD[0] += transIQ[nstageI];
  
  DcasesI += transP[nstageP-1];
  DcasesIQ += transPQ[nstageP-1] + transP[nstageP+1];
  DcasesH += transPQ[nstageP] + transP[nstageP];
  DdeathsIIQ += transI[nstageI] + transIQ[nstageI];
  DdeathsCV += transC[nstageC+1] + transV[nstageV];
}
// Note: for translating "double *SUSC = &S" we directly fill S and no need to define SUSC
snippet.initializer = function(args, covar) {
  let initObj = {};
  
  let fS = args.S0;
  let fEQ = args.EQ0 / nstageE;
  let fPQ = args.PQ0 / nstageP;
  let fIQ = args.IQ0 / nstageI;
  let fE = args.E0 / nstageE;
  let fP = args.P0 / nstageP;
  let fI = args.I0 / nstageI;
  let fH = args.H0 / nstageH;
  let fC = args.C0 / nstageC;
  let fV = args.V0 / nstageV;
  let fM = args.M0;
  let fR = 1 - fS - nstageE*fE - nstageP*fP - nstageI*fI - nstageH*fH - nstageC*fC - nstageV*fV - nstageE*fEQ - nstageP*fPQ - nstageI*fIQ - fM;
  
  initObj["S"] = Math.round(pop*fS);
  for (let i = 0; i < nstageE; i++) initObj[`EQ${i + 1}`] = Math.round(pop*fEQ);
  for (let i = 0; i < nstageP; i++) initObj[`PQ${i + 1}`] = Math.round(pop*fPQ);
  for (let i = 0; i < nstageI; i++) initObj[`IQ${i + 1}`] = Math.round(pop*fIQ);
  for (let i = 0; i < nstageE; i++) initObj[`E${i + 1}`] = Math.round(pop*fE);
  for (let i = 0; i < nstageP; i++) initObj[`P${i + 1}`] = Math.round(pop*fP);
  for (let i = 0; i < nstageI; i++) initObj[`I${i + 1}`] = Math.round(pop*fI);
  for (let i = 0; i < nstageH; i++) initObj[`H${i + 1}`] = Math.round(pop*fH);
  for (let i = 0; i < nstageC; i++) initObj[`C${i + 1}`] = Math.round(pop*fC);
  for (let i = 0; i < nstageV; i++) initObj[`V${i + 1}`] = Math.round(pop*fV);
  initObj["M"] = Math.round(pop*fM);
  initObj["R"] = Math.round(pop*fR);
  
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
    let reported_deaths = states.deathsV + states.deathsI;
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

