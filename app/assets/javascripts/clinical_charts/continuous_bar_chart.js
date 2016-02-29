function chart_c(data, title, y_axis) {
	var distribution = generate_distribution(data);
	var generated_bar = generate_continuous_bar(distribution, data);
	
	$('#container').highcharts({
		chart: {
			type: 'column'
		},
		title: {
			text: title
		},
		xAxis: {
			categories: generated_bar['labels']
		},
		yAxis: {
			title: {
				text: y_axis
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