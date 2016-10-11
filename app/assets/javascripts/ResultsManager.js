function ResultsManager(tumor_type, data_source, prediction_target, partition, var1, var2, var3, var4, session_id) {
	this.tumor_type = tumor_type;
	this.data_source = data_source;
	this.prediction_target = prediction_target;
	this.partition = partition;
	if (var1=="-1") {
		var1 = "N/A";
	}
	if (var2=="-1") {
		var2 = "N/A";
	}
	if (var3=="-1") {
		var3 = "N/A";
	}
	if (var4=="-1") {
		var4 = "N/A";
	}
	this.var1 = var1;
	this.var2 = var2;
	this.var3 = var3;
	this.var4 = var4;
	this.session_id = session_id;
}


ResultsManager.prototype.run = function(type) {
	var self = this;

	$("#progressbar").progressbar({ value: 0 });
	$("#progressbar-status").text("Initializing...");

	var interval = setInterval(function() {
		$.get( "/results/progress", {session_id: self.session_id}, function(data) {
			$("#progressbar").progressbar({ value: data.percent });
			if (data.done) {
				clearInterval(interval);
				location.href = "/analysis/" + type + "?done=1&tumor_type=" + self.tumor_type + "&data_source=" + self.data_source + "&prediction_target=" + self.prediction_target + "&partition=" + self.partition + "&var1=" + self.var1 + "&var2=" + self.var2 + "&var3=" + self.var3 + "&var4=" + self.var4 + "&session_id=" + self.session_id + "#results";
			}

			if (data.status) {
				$("#progressbar-status").text(data.status);
			}

			if (data.error) {
				$("#progressbar-status").html('<div class="text-error">ERROR: <br><br>' + data.error);
				clearInterval(interval);
		  	}
		});
	}, 2000);
};

ResultsManager.prototype.run_stringdb = function() {
	var feature_weights = [];
	
	$('input[name="feature_weights"]:checked').each(function() {
		feature_weights.push($(this).val());
	});
	$('#stringdb_input').val(feature_weights.join("\n"));
	$('#stringdb_form').submit();
};
