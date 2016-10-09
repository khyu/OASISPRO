function ResultsManager(session_id) {
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
				location.href = "/analysis/" + type + "?done=1&session_id=" + self.session_id + "#results";
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
