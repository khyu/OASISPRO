class AdminController < ApplicationController
	def index
	end

	def run_downloader
		@running = false
		@password_incorrect = false
		if params[:password] == "oasispro"
			@running = true
			system("bash public/administrator/updateDatabase.sh")
		else
			@password_incorrect = true
		end

		render :index
	end
end