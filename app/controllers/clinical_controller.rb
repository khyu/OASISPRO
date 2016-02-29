class ClinicalController < ApplicationController
	require 'site_constants.rb' 

	def index
	end

	def get_clinical_variables
		f = File.open("../data/nationwidechildrens.org_clinical_patient_" + params[:tumor_type] + ".txt", "r")
		data = get_variable_names_row(f).gsub("_", " ").split("\t")

		f.close

		render json: data.sort_by!{ |var| var.upcase }
	end

	def chart_data
		data = Array.new
		
		clinical_variable = params[:clinical_variable].gsub(" ", "_")
		f = File.open("../data/nationwidechildrens.org_clinical_patient_" + params[:tumor_type] + ".txt", "r")
		clinical_variable_index = get_col_index(clinical_variable, f)
		f.close

		g = File.open("../data/nationwidechildrens.org_clinical_patient_" + params[:tumor_type] + ".txt", "r")
		g.each_line do |line|
		  line = line.strip.split("\t")
		  data << line[clinical_variable_index]
		end
		g.close

		render json: {data: data[3..-1], chart_type: SiteConstants::CLINICAL_VARIABLE_TYPES[clinical_variable]}
	end

	private # -------------------------------------

	# returns the line in the file containing the varibale names
	def get_variable_names_row(file)
		x = 0
		file.each_line do |line|
			if x == 1
				return line
			end
			x += 1
		end
	end

	def get_col_index(target_var_name, file)
		var_names_row = get_variable_names_row(file)
		variables = var_names_row.split(' ')
		index = 0
		for var in variables
			if var.downcase == target_var_name.downcase
				return index
			end
			index += 1
		end
	end

end