let rootDir ='.'
const fs = require('fs');
const { trajMatch } = require('../src/trajMatch.js');

const create_dataset = require('./CreateDataset.js');
const create_covars = require('./CreateCovars.js');
const snippet = require('./modelSnippet_DetModel3.js');

let rootDirData ='./private_data';

run = 1;
job = 2;


total_cores = 1;
no_points = 1;

no_cores = total_cores;
job = job - (Math.ceil(job/no_cores)-1)*no_cores;

params_fixed = ["beta_sd","dB0", "dB1","sigma","kappa"];

// read all rows and chaeck the border time and convert the selected colnames
let temp, file;
file = fs.readFileSync(rootDirData+'/ParamSet_run1.csv').toString();
let lines = file.split(/\r\n|\n/);

let paramsetData = [];
let paramsetHeader = lines[0].replace(/['"]+/g, '').split(',');
for (let i = 1; i < lines.length; i++) {
  let temp = lines[i].split(',');
  if(temp.length > 1) {
    let tempParamset =	{};
      for(let j = 0; j < temp.length; j++){
        tempParamset[paramsetHeader[j]] = Number(temp[j]);
      }
      paramsetData.push(tempParamset);
    
  }
}

select_set = [...snippet.paramsMod, ...snippet.paramsIc];

// Generate covars, data and pomp object
data = create_dataset(snippet.endTime)
covars = create_covars(snippet.endTime)

globals = { nstageE: 3, nstageP: 3, nstageI: 3, nstageH: 3, nstageC: 3, nstageV: 3, pop: 10e6, T0: 75, T1: 139 };

let dataHeader = data.shift();
let dataCases = data.map((x,i,arr) => { 
  a = {};
  a[dataHeader[1]] = x[1];
  a[dataHeader[2]] = x[2];
  a[dataHeader[3]] = x[3];
  a[dataHeader[4]] = x[4];
  a[dataHeader[5]] = x[5];
  return a;});
let dataCasesTimes = data.map((x,i,arr) => x[0]);
dataHeader.shift();

let covarHeader = covars.shift();
let dataCovar = covars.map((x,i,arr) => { 
  a = {};
  a[covarHeader[1]] = x[1];
  return a;});
let dataCovarTimes = covars.map((x,i,arr) => x[0]);
covarHeader.shift();

const pompData = {
  data :  dataCases,
  times:  dataCasesTimes,
  t0: 0,
  skeletonDetail:  { type:"map", deltaT: 0.1 },
  covar: dataCovar,
  tcovar: dataCovarTimes,
  zeronames: snippet.zeronames,
  statenames: snippet.statenames,
  paramnames: [...snippet.paramsMod, ...snippet.paramsIc],
  covarnames: covarHeader,
  obsnames: dataHeader,
  globals: globals,
};

let tm = trajMatch(paramsetData[0],{object: pompData, est: [],transform: true, method: "subplex"})

console.log('finished.');