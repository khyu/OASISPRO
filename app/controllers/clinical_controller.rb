class ClinicalController < ApplicationController
	require 'site_constants.rb' 

	def index
	end

	def get_clinical_variables
		data = []
		f = File.open("../data/" + params[:tumor_type] + "-clinical-cut.txt", "r")
		x = 0
		f.each_line do |line|
			if x > 0
				data = line.strip.replace("_", " ").capitolize.split("\t")
				break
			end
			x += 1
		end

		render json: data
	end

	def chart_data
		data = Array.new
		
		f = File.open("../data/" + params[:tumor_type] + "-clinical-cut.txt", "r")
		f.each_line do |line|
		  line = line.strip.split("\t")
		  data << line[params[:chart].to_i - 1]
		end
		f.close
		
		render json: data[3..-1]
	end
end
