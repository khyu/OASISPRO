/*function chart7(data) {
	$('#container').highcharts({
		title: {
			text: 'Last Contact Years To'
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