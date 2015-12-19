function chart1(data) {
	// [start, end, increment]
	var distribution = [20,90,10];
	var generated_bar = generate_continuous_bar(distribution, data, function(data) { return data/-365; })
	
	$('#container').highcharts({
		chart: {
			type: 'column'
		},
		title: {
			text: 'Age'
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
				pointPadding: 0,
				borderWidth: 0,
				groupPadding: 0
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