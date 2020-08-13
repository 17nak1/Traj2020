let rootDir ='.'
const fs = require('fs');
const { trajMatch } = require('../src/trajMatch.js');

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
      for(let j = 1; j < temp.length; j++){
        tempParamset[paramsetHeader[j]] = Number(temp[j]);
      }
      paramsetData.push(tempParamset);
    
  }
}

select_set = [...snippet.params_mod, ...snippet.params_ic];

console.log('finished.');