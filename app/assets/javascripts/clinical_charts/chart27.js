function chart27(data) {
	var NA = 0;
	var unknown = 0;
	var yes = 0;
	var no = 0;
	console.log(data);
	for (var x in data) {
		if (data[x] == '[Not Available]') {
			NA++;
		}
		if (data[x] == '[Unknown]') {
			unknown++;
		}
		if (data[x] == 'YES') {
			yes++;
		}
		if (data[x] == 'NO') {
			no++;
		}
	}
	
	$('#container').highcharts({
		title: {
			text: 'Pharmaceutical tx Adjuvant'
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
				['Yes', yes],
				['No', no]
			]
		}]
	});
}