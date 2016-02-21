function discreteBarChart(data, title, y_axis) {
	$('#container').highcharts({
		chart: {
			type: 'column'
		},
		title: {
			text: title
		},
		xAxis: {
			categories: get_bar_X_axis_categories(data)
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
			data: generate_bar_dataPoints(get_bar_X_axis_categories(data), data)
		}]
	});	
}