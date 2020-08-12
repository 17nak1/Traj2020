
snippet = {}
let mathLib = require('./mathLib')

snippet.skeleton = function (states, params, t, dt, covar) {
  let seas, dy = {};
  
  let beta0 = params.R0 * (params.gamma + params.mu) * (params.sigma + params.mu) / params.sigma;
  let va = 0
  
  let tt = (t - Math.floor(t)) * 365.25
  if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
    seas = 1 + params.amplitude * 0.2411 / 0.7589
  } else {
    seas = 1 - params.amplitude
  }
  let Beta = beta0 * seas / covar.pop
  dy.S = covar.birthrate * (1 - va) - Beta * states.S * states.I - params.mu * states.S
  dy.E = Beta * states.S * states.I - (params.sigma + params.mu) * states.E
  dy.I = params.gamma * states.I - params.mu * states.R + covar.birthrate * va
  dy.R = params.sigma * states.E - (params.gamma + params.mu) * states.I
  dy.H = params.gamma * states.I
  return dy;
}

snippet.initializer = function(states, covar) {
  
  let m = covar.pop / (states.S_0 + states.E_0 + states.R_0 + states.I_0);
  let S = Math.round(m * states.S_0);
  let E = Math.round(m * states.E_0);
  let I = Math.round(m * states.I_0);
  let R = Math.round(m * states.R_0);
  let H = 0;
  return {S: S, E: E, I: I, R: R, H: H};
}

snippet.dmeasure = function (data ,hiddenState, params, giveLog = 1) {
  let lik
  let rho = params.rho;
  let psi = params.psi;
  let H = hiddenState.H;
  let cases = data.cases;

  let tol = 1.0e-18
  let mn = rho * H;
  let v = mn * (1.0 - rho + psi * psi * mn);
  
  let modelCases = Number(cases);
  if(!isNaN(modelCases)){
    if (modelCases > 0.0) {
      lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
    } else {
      lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol;
    }
    if (giveLog) lik = Math.log(lik);
  } else {
    lik = (giveLog) ? 0 : 1;
  }
  return lik
}


snippet.paramsMod = ["R0","amplitude","gamma","mu","sigma","rho","psi"];
snippet.paramsIc = ["S_0", "E_0", "I_0", "R_0"];
snippet.zeronames = ["H"];
snippet.statenames = ["S","E","I","R","H"];

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

module.exports = snippet