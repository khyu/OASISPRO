function chart4(data) {
	$('#container').highcharts({
		title: {
			text: 'Ethnicity'
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