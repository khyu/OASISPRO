/* pld = Persistent locoregional disease */
/* nieod = No imaging evidence of disease */
/* pdm = Persistent distant metastases */
function chart28(data) {
	var NA = 0;
	var unknown = 0;
	var pld = 0;
	var nieod = 0;
	var pdm = 0;
	console.log(data);
	for (var x in data) {
		if (data[x] == '[Not Available]') {
			NA++;
		}
		if (data[x] == '[Unknown]') {
			unknown++;
		}
		if (data[x] == 'Persistent locoregional disease') {
			pld++;
		}
		if (data[x] == 'No imaging evidence of disease') {
			nieod++;
		}
		if (data[x] == 'Persistent distant metastases') {
			pdm++;
		}
	}
	
	$('#container').highcharts({
		title: {
			text: 'Clinical Status Within 3 Months of Surgery'
		},
		tooltip: {
			pointFormat: '<b>{point.percentage:.1f}%</b>'
		},
		plotOptions: {
			pie: {
				allowPointSelect: true,
				cursor: 'pointer',
				dataLabels: {
					enabled: true
				}
			}
		},
		series: [{
			type: 'pie',
			name: '',
			data: [
				['[Not Available]', NA],
				['[Unknown]', unknown],
				['Persistent locoregional disease', pld],
				['No imaging evidence of disease', nieod],
				['Persistent distant metastases', pdm]
			]
		}]
	});
}