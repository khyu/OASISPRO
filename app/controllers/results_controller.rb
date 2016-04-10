class ResultsController < ApplicationController
	def index
	end

	def progress
		return if params[:session_id] !~ /\A[0-9]+\z/

		milestones = []
		done = false

		puts "public/sessions/#{params[:session_id]}/milestones.txt"
		t = Time.now.to_i
		while milestones.length <= params[:lines].to_i && Time.now.to_i - t < 10
			if File.exists?("public/sessions/#{params[:session_id]}/milestones.txt")
				milestones = File.open("public/sessions/#{params[:session_id]}/milestones.txt", "rb").read.split("\n")

				if milestones.last == 'Completed!'
					done = true
					break
				elsif milestones.length > params[:lines].to_i
					break
				end
			end

			sleep 1
			break
		end

		render json: {milestones: milestones, done: done}
	end
end
