function BatchParamsManager() {
	var change_partition = function() {
		$('.random-partition-params').toggle($("input[name=partition]:checked").val() == "random");
		$('.batch-partition-params').toggle($("input[name=partition]:checked").val() == "batch");
		$('.kfold-params').toggle($("input[name=partition]:checked").val() == "kfold");
	};
	var specify_seed = function() {
		$('.random-partition-params input[name=random_seed]').toggle($("input[name=seed_type]:checked").val() == "specify_seed");
	}

	$("input[name=partition]").click(change_partition);
		change_partition();

	$("input[name=seed_type]").click(specify_seed);
		specify_seed();

	$('#tumor_type').on('change', function() {
		var tbl = document.getElementById('batch-options-table');
		for (var i = tbl.rows.length - 1; i >= 0; i--) {
			$(tbl.rows[i]).hide();
		};
		$.get("/analysis/batch_ids", { tumor_type: $('#tumor_type').val() }, function(data) {
			// return data as a json array of batch-ids
			for (var i = data['ids'].length - 1; i >= 0; i--) {
				$('#batch-' + data['ids'][i].toString()).show();
			};
		});
	});
}