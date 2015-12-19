class ClinicalController < ApplicationController
	def index
	end
	
	def chart_data
		data = Array.new
		
		f = File.open("public/Clinical-Cut.txt", "r")
		f.each_line do |line|
		  line = line.strip.split("\t")
		  data << line[params[:chart].to_i - 1]
		end
		f.close
		
		render json: data[3..-1]
	end
end
