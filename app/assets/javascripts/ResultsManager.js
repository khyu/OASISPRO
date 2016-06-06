function ResultsManager(session_id) {
	this.session_id = session_id;
}


ResultsManager.prototype.run = function(type) {
	//$("#results-modal").modal();
	var self = this;

	$("#progressbar").progressbar({ value: 0 });
	$("#progressbar-status").text("Initializing...");

	var interval = setInterval(function() {
		$.get( "/results/progress", {session_id: self.session_id}, function(data) {
			//$("#results-modal-content").html(data.milestones.join("<br>"));
			console.log(data);
			$("#progressbar").progressbar({ value: data.percent });
			if (data.done) {
				clearInterval(interval);
				if (type == 'stage') {
					location.href = "/analysis/stage?done=1&session_id=" + self.session_id + "#results";
				} else if (type == 'survival') {
					location.href = "/analysis/survival?done=1&session_id=" + self.session_id + "#results";
				}
				
			}

			$("#progressbar-status").text(data.status);
			if (data.error) {
				//$("#results-modal-content").html("ERROR: <br><br>" + data.error);
				$("#progress-status").html("ERROR: <br><br>" + data.error);
				clearInterval(interval);
		  	}
		});
	}, 2000);
};

ResultsManager.prototype.run_stringdb = function() {
	var feature_weights = [];
	
	$('input[name="feature_weights"]:checked').each(function() {
		console.log($(this));
		feature_weights.push($(this).val());
	});
	$('#stringdb_input').val(feature_weights.join("\n"));
	$('#stringdb_form').submit();
};
