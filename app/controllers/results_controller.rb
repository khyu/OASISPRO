class ResultsController < ApplicationController
	def index
	end

	def progress
		return if params[:session_id] !~ /\A[0-9]+\z/

		milestones = []
		percent = 0
		done = false
		error_file = ""

		if File.exists?("public/sessions/#{params[:session_id]}/error.txt")
			error_file = File.open("public/sessions/#{params[:session_id]}/error.txt", "r").read
		end

		error_file = "" if !(error_file.include?('Execution halted') || error_file.include?('Rscript: command not found'))
	
		if File.exists?("public/sessions/#{params[:session_id]}/milestones.txt")
			milestones = File.open("public/sessions/#{params[:session_id]}/milestones.txt", "r").read.split("\n")

			percent = milestones.last.split(',')[1].to_i

			if milestones.last == 'Completed!,100'
				done = true
			end


		end		

		render json: {percent: percent, done: done, error: error_file.gsub("\n", "<br>")}
	end
end
