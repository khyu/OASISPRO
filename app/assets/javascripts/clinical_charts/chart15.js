function chart15(data) {
	var distribution = [0,10,1];
	var generated_bar = generate_continuous_bar(distribution, data);
	
	$('#container').highcharts({
		chart: {
			type: 'column'
		},
		title: {
			text: 'Tumor Size Width'
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