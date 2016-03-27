class ClinicalController < ApplicationController
	require 'site_constants.rb' 

	def index
	end

	def get_clinical_variables
		f = File.open("../data/nationwidechildrens.org_clinical_patient_" + params[:tumor_type] + ".txt", "r")
		data = get_variable_names_row(f).gsub("_", " ").split("\t")
		num_records = get_num_records(f)
		f.close
		data = data.sort_by!{ |var| var.upcase }
		render json: {"tumor_type" => SiteConstants::TUMOR_TYPES[params[:tumor_type].to_sym], "num_records" => num_records, "vars" => data}
	end

	def chart_data
		data = get_data(params[:clinical_variable], params[:tumor_type])

		render json: {data: data[3..-1], chart_type: SiteConstants::CLINICAL_VARIABLE_TYPES[params[:clinical_variable].gsub(" ", "_")]}
	end

	def table_data
		data = get_data(params[:clinical_variable], params[:tumor_type])
		render :partial => "data_table", :locals => { data: data_table_data(data)}
	end

	private # -------------------------------------

	def get_data(clinical_variable, tumor_type)
		data = Array.new
		clinical_variable = clinical_variable.gsub(" ", "_")
		f = File.open("../data/nationwidechildrens.org_clinical_patient_" + tumor_type + ".txt", "r")
		clinical_variable_index = get_col_index(clinical_variable, f)
		f.close

		g = File.open("../data/nationwidechildrens.org_clinical_patient_" + tumor_type + ".txt", "r")
		g.each_line do |line|
		  line = line.strip.split("\t")
		  data << line[clinical_variable_index]
		end
		g.close
		return data
	end

	def data_table_data(data)
		elems = {}
		for line in data
			if !elems[line]
				elems[line] = 0
			end
			elems[line] = elems[line] + 1
		end
		data_array = Array.new()
		for elem in elems
			data_array << [elem, elems[elem]]
		end
		return data_array
	end

	def get_num_records(file)
		x = 0
		file.each_line do |line|
			x += 1
		end
		return x - 1
	end

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