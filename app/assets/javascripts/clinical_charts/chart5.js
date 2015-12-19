function chart5(data) {
	$('#container').highcharts({
		title: {
			text: 'History Neoadjuvant Treatment'
		},
		tooltip: {
			pointFormat: '<b>{point.y} ({point.percentage:.1f}%)</b>'
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
			data: generate_pie_categories(data)
		}]
	});
}

/*
function chart5(data) {
	$('#container').highcharts({
		chart: {
			type: 'column'
		},
		title: {
			text: 'History Neoadjuvant Treatment'
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
*/