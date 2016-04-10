class ResultsController < ApplicationController
	def index
	end

	def progress
		return if params[:session_id] !~ /\A[0-9]+\z/

		milestones = []
		done = false
		error_file = ""

		t = Time.now.to_i
		while milestones.length <= params[:lines].to_i && Time.now.to_i - t < 10
			if File.exists?("public/sessions/#{params[:session_id]}/error.txt")
				error_file = File.open("public/sessions/#{params[:session_id]}/error.txt", "r").read
			end

			if error_file.include?('Execution halted') || error_file.include?('Rscript: command not found')
				break
			else
				error_file = ""
			end

			if File.exists?("public/sessions/#{params[:session_id]}/milestones.txt")
				milestones = File.open("public/sessions/#{params[:session_id]}/milestones.txt", "r").read.split("\n")

				if milestones.last == 'Completed!'
					done = true
					break
				elsif milestones.length > params[:lines].to_i
					break
				end
			end

			sleep 1
		end

		render json: {milestones: milestones, done: done, error: error_file.gsub("\n", "<br>")}
	end
end
