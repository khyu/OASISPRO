function ResultsManager(session_id) {
	this.lines = 0;
	this.session_id = session_id;
}


ResultsManager.prototype.run = function() {
	$("#results-modal").modal();
	var self = this;

	var interval = setInterval(function() {
		$.get( "/results/progress", {session_id: self.session_id, lines: self.lines}, function(data) {
		  $("#results-modal-content").html(data.milestones.join("<br>"));
		  if (data.done) {
		  	location.href = "/analysis/stage?done=1&session_id=" + self.session_id + "#results";
		  }
		  if (data.error) {
		  	$("#results-modal-content").html("ERROR: <br><br>" + data.error);
		  	clearInterval(interval);
		  }
		});
	}, 2000);
}