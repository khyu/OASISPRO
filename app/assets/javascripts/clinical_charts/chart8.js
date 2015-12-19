/*function chart8(data) {
	$('#container').highcharts({
		title: {
			text: 'Death Days To'
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
}*/

function chart8(data) {
	var distribution = [0,400,800,1200,1600,2000];
	var generated_bar = generate_continuous_bar(distribution, data)
	
	$('#container').highcharts({
		chart: {
			type: 'column'
		},
		title: {
			text: 'Death Days To'
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