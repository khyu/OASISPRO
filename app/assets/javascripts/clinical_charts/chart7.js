function chart7(data) {
	var distribution = [0,.5,1,2,3,4,5,6,7,8,9,10,20];
	var generated_bar = generate_continuous_bar(distribution, data, function(data) { return data/365; })
	
	$('#container').highcharts({
		chart: {
			type: 'column'
		},
		title: {
			text: 'Last Contact Years To'
		},
		xAxis: {
			categories: generated_bar['labels']
		},
		yAxis: {
			title: {
				text: 'Number of patients'
			}
		},
		plotOptions: {
			column: {
				pointPadding: 0.2,
				borderWidth: 0
			}
		},
		tooltip: {
			pointFormat: '{point.y}'
		},
		series: [{
			showInLegend: false,
			data: generated_bar['groups']
		}]
	});	
}