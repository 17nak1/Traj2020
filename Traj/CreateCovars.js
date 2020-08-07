const fs = require('fs');
let rootDirData ='./private_data';

let create_covars = function(endTime="2020-04-27",predTime=null){
	let t0 = new Date("2019-12-31");
	let tf = new Date(endTime);

	// read all rows and chaeck the border time and convert the selected colnames
	let temp, file;
	file = fs.readFileSync(rootDirData+'/covidtesting.csv').toString();
	let lines = file.split(/\r\n|\n/);
	let selected_colnames = ["Reported Date",
		"Total patients approved for testing as of Reporting Date",
		"Under Investigation"];
	let dataCovar_Index = [];	
	let dataCovar = [["Date","TotalTests","Pending"]];
	let dataCovar_temp = lines[0].replace(/['"]+/g, '').split(',');
	for(let i = 0; i < selected_colnames.length; i++){
		dataCovar_Index.push(dataCovar_temp.indexOf(selected_colnames[i]));
	}
  let startdate = new Date("2020-02-03");
	for (let i = 1; i < lines.length; i++) {
		let temp = lines[i].split(',');
		if(temp.length > 1) {
			tempDate = new Date(temp[dataCovar_Index[0]]);
			if( tempDate > startdate && tempDate <= tf){
				tempDayCount = Math.ceil(Math.abs(tempDate - t0) / (1000 * 3600 * 24));
			  let tempcovar =	[tempDayCount];
				for(let j = 1; j < selected_colnames.length; j++){
					tempcovar.push(Number(temp[dataCovar_Index[j]]));
				}
				dataCovar.push(tempcovar);
			}
		}
	}

	time = new Array(dataCovar[dataCovar.length-1][dataCovar_Index[0]]).fill(0)
	var timeDiff = Math.abs(date2.getTime() - date1.getTime());
	var diffDays = Math.ceil(timeDiff / (1000 * 3600 * 24)); 
	alert(diffDays);	
}

create_covars();
console.log("finish");

module.exports = create_covars;