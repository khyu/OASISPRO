function chart9(data) {
	$('#container').highcharts({
		chart: {
			type: 'column'
		},
		title: {
			text: 'History Thyroid Disease'
		},
		xAxis: {
			categories: get_bar_X_axis_categories(data)
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
			data: generate_bar_dataPoints(get_bar_X_axis_categories(data), data)
		}]
	});	
}