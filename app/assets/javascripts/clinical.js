function get_data(tumor_type, clinical_variable) {
	$.get('/clinical/chart_data', {tumor_type: tumor_type, clinical_variable: clinical_variable}, function (data) {
		window[chart_type](data);
	});
}

function generate_pie_categories(data) {
	var histogram = {};
	
	for (var x in data) {
		if (!histogram[data[x]]) {
			histogram[data[x]] = 0;
		}
		histogram[data[x]]++;
	}
	
	var categories = [];
	for (var h in histogram) {
		categories.push([h, histogram[h]]);
	}
	
	return categories;
}

/* Returns an array of all of the data values in the data */
function get_bar_X_axis_categories(data) {
	var categories = [];
	for (var x in data) {
		if (categories.indexOf(data[x]) == -1) {
			categories.push(data[x])
		}
	}
	categories.sort();
	return categories;
}

function generate_bar_dataPoints(categories, data) {
	dataPoints = [];
	indexMap = {};
	currIndex = 0;
	for (var x in categories) {
		indexMap[categories[x]] = currIndex;
		dataPoints[currIndex] = 0;
		currIndex++;
	}
	for (var y in data) {
		dataPoints[indexMap[data[y]]]++;
	}
	console.log(dataPoints);
	return dataPoints;
}

function generate_continuous_bar(distribution, data, transform) {
	var labels = ['N/A'];
	var groups = [0];
	
	var bins = [];
	if (distribution.length == 3) {
		for (var x = distribution[0]; x <= distribution[1]; x += distribution[2]) {
			bins.push(x);
		}
	}
	else {
		bins = distribution;
	}
	
	for (var x = 1; x < bins.length; x++) {
		groups.push(0);
		labels.push('[' +bins[x-1] + ', ' + bins[x]+')');
	}
	
	for (var x in data) {
		if (parseFloat(data[x]) != data[x]) {
			groups[0]++;
		}
		else {
			if (typeof transform == 'function') {
				data[x] = transform(data[x]);
			}
			
			for (var y = 1; y < bins.length; y++) {
				if (data[x]*1 < bins[y]) {
					groups[y]++;
					break;
				}
			}
		}
	}
	
	if (groups[0] == 0) {
		labels = labels.slice(1);
		groups = groups.slice(1);
	}
	
	return {'labels': labels, 'groups': groups}	
}
