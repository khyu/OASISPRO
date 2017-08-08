class AnalysisController < ApplicationController
	require 'site_constants.rb'
	require 'mail'

	def stage
		if params[:done]
			@tumor_type = "#{params[:tumor_type]}"
			@data_source = "#{params[:data_source]}"
			@prediction_target = "#{params[:prediction_target]}"
			@partition = "#{params[:partition]}"
			@var1 = "#{params[:var1]}"
			@var2 = "#{params[:var2]}"
			@var3 = "#{params[:var3]}"
			@var4 = "#{params[:var4]}"
			@sessionID = "#{params[:session_id]}"

			# number of samples in the training set
			@nSamplesTraining = Array.new
			f = File.open("public/sessions/#{params[:session_id]}/nSamplesTraining.txt", "r")
			f.each_line do |line|
				@nSamplesTraining << line
			end
			f.close

			# number of samples in the test set
			@nSamplesTest = Array.new
			f = File.open("public/sessions/#{params[:session_id]}/nSamplesTest.txt", "r")
			f.each_line do |line|
				@nSamplesTest << line
			end
			f.close

			# group 1
			@pred_group1 = Array.new
			f = File.open("public/sessions/#{params[:session_id]}/outcomeLabel1.txt", "r")
			f.each_line do |line|
				@pred_group1 << line
			end
			f.close

			# group 2
			@pred_group2 = Array.new
			f = File.open("public/sessions/#{params[:session_id]}/outcomeLabel2.txt", "r")
			f.each_line do |line|
				@pred_group2 << line
			end
			f.close

			# area under curve
			@aucs = get_aucs(params[:session_id])

			# feature weights
			@feature_weights = get_feature_weights(params[:session_id])
		end

		if params[:generate]
			return if params[:session_id] !~ /\A[0-9]+\z/

			validators = {
				tumor_type: SiteConstants::TUMOR_TYPES.keys.map { |x| x.to_s },
				data_source: SiteConstants::DATA_TYPES.keys.map { |x| x.to_s },
				prediction_target: SiteConstants::PREDICTION_TARGTS.keys.map { |x| x.to_s },
				partition: SiteConstants::PARTITION_TYPES,
				random_seed: '*',
				training_percentage: 'FLOAT',
				feature_selection_method: SiteConstants::FEATURE_SELECTION_METHOD.keys.map { |x| x.to_s },
				num_top_features: 'INTEGER',
				session_id: '*',
				email: '*'
			}

			if params[:partition] == 'kfold'
				params[:random_seed] = params[:num_folds]
			end

			if params[:partition] == 'batch'
				params[:random_seed] = "batch.txt"
				Dir.mkdir("public/sessions") unless File.exists?("public/sessions")
				Dir.mkdir("public/sessions/#{params[:session_id]}") unless File.exists?("public/sessions/#{params[:session_id]}")
				File.write("public/sessions/#{params[:session_id]}/" + params[:random_seed], (params[:medical_centers] || []).join("\n") + "\n")
			else
				if params[:random_seed].blank?
					params[:random_seed] = -1
				end
				Dir.mkdir("public/sessions") unless File.exists?("public/sessions")
				Dir.mkdir("public/sessions/#{params[:session_id]}") unless File.exists?("public/sessions/#{params[:session_id]}")
			end

			File.write("public/sessions/#{params[:session_id]}/outcomeLabel1.txt", params[:prediction_target_group1].join("\n") + "\n")
			File.write("public/sessions/#{params[:session_id]}/outcomeLabel2.txt", params[:prediction_target_group2].join("\n") + "\n")

			
			File.write("public/sessions/#{params[:session_id]}/mlparameters.txt", (params[:ml_methods1] || "off") + "\n" + (params[:ml_methods2] || "off") + "\n" + (params[:ml_methods3] || "off") + "\n" + (params[:ml_methods4] || "off") + "\n" + (params[:ml_methods5] || "off") + "\n" + (params[:ml_methods6] || "off") + "\n" + (params[:ml_methods7] || "off") + "\n" + (params[:ml_methods8] || "off") + "\n" + (params[:ml_methods9] || "off") + "\n" + (params[:ml_methods10] || "off") + "\n")
			File.write("public/sessions/#{params[:session_id]}/mlparameters.txt", params[:ml1laplace] + "\n" + params[:ml2treedep] + "\n" + params[:ml2cp] + "\n" + params[:ml3treedep] + "\n" + params[:ml4bs] + "\n" + params[:ml5ntrees] + "\n" + params[:ml6ntrees] + "\n" + params[:ml7costmin] + "\n" + params[:ml7costmax] + "\n" + params[:ml8costmin] + "\n" + params[:ml8costmax] + "\n" + params[:ml9costmin] + "\n" + params[:ml9costmax] + "\n" + params[:ml10costmin] + "\n" + params[:ml10costmax] + "\n", mode: 'a')


			if params[:training_percentage].present?
				params[:training_percentage] = Float(params[:training_percentage])/100
			else
				params[:training_percentage] = -1
			end

			if ['LOOCV', 'batch', 'kfold'].include?(params[:partition])
				params[:training_percentage] = -1
			end

			result = build_command(validators, params)
			@valid_command = result[:valid_command]
		
			if @valid_command
				@command = "Rscript --max-ppsize=500000 public/RCodes/binaryClassification.R#{result[:command]} 2>public/sessions/#{params[:session_id]}/error.txt &"
				system(@command)
			end
		end
	end

	def survival
		if params[:done]
			@tumor_type = "#{params[:tumor_type]}"
			@data_source = "#{params[:data_source]}"
			@partition = "#{params[:partition]}"
			@var1 = "#{params[:var1]}"
			@var2 = "#{params[:var2]}"
			@var3 = "#{params[:var3]}"
			@var4 = "#{params[:var4]}"
			@email = "#{params[:email]}"
			@sessionID = "#{params[:session_id]}"
			

			# number of samples
			@nSamples = Array.new
			f = File.open("public/sessions/#{params[:session_id]}/nSamples.txt", "r")
			f.each_line do |line|
				@nSamples << line
			end
			f.close

			# p value
			@pTest = Array.new
			f = File.open("public/sessions/#{params[:session_id]}/pTest.txt", "r")
			x = 0
			f.each_line do |line|
				@pTest << [line]
			end
			f.close
			@feature_weights = get_feature_weights(params[:session_id])
		end

		if params[:generate]
			params[:session_id] = rand.to_s.sub("0.", "") 
			validators = {
				tumor_type: SiteConstants::TUMOR_TYPES.keys.map { |x| x.to_s },
				data_source: SiteConstants::DATA_TYPES.keys.map { |x| x.to_s },
				partition: SiteConstants::PARTITION_TYPES,
				random_seed: '*',
				training_percentage: 'FLOAT',
				alpha_lower_bound: 'FLOAT',
				alpha_upper_bound: 'FLOAT',
				lambda_lower_bound: 'FLOAT',
				lambda_upper_bound: 'FLOAT',
				clinical_variable_file: '*',
				session_id: '*',
				email: '*'
			}

			if params[:partition] == 'batch'
				params[:random_seed] = "batch.txt"
				Dir.mkdir("public/sessions") unless File.exists?("public/sessions")
				Dir.mkdir("public/sessions/#{params[:session_id]}") unless File.exists?("public/sessions/#{params[:session_id]}")
				File.write("public/sessions/#{params[:session_id]}/" + params[:random_seed], (params[:medical_centers] || []).join("\n") + "\n")
			else
				if params[:partition] == 'kfold'
					params[:random_seed] = params[:num_folds]
				end
				if params[:random_seed].blank?
					params[:random_seed] = -1
				end
				Dir.mkdir("public/sessions") unless File.exists?("public/sessions")
				Dir.mkdir("public/sessions/#{params[:session_id]}") unless File.exists?("public/sessions/#{params[:session_id]}")
			end

			if params[:training_percentage].present?
				params[:training_percentage] = Float(params[:training_percentage])/100
			else
				params[:training_percentage] = -1
			end

			# Default for alpha/lambda bound params
			validators.keys[5..8].each do |param|
				params[param] = -1 if params[param].blank?
			end

			if params[:clinical_variables]
				params[:clinical_variable_file] = "variables.txt"
				Dir.mkdir("public/sessions") unless File.exists?("public/sessions")
				Dir.mkdir("public/sessions/#{params[:session_id]}") unless File.exists?("public/sessions/#{params[:session_id]}")
				File.write("public/sessions/#{params[:session_id]}/" + params[:clinical_variable_file], params[:clinical_variables].join("\n") + "\n")
			else
				params[:clinical_variable_file] = -1
			end

			result = build_command(validators, params)
			@valid_command = result[:valid_command]

			@command = "Rscript --max-ppsize=500000 public/RCodes/elasticNetCox.R#{result[:command]} 2>public/sessions/#{params[:session_id]}/error.txt &"
			system(@command)
			
			@sessionID = "#{params[:session_id]}"
		end
	end

	def get_features
		source = params[:source]
		method = params[:method]

		if source && method
			lines = Array.new
			f = File.open("public/feature_selection/"+source+"-"+method+".txt", "r")
			f.each_line do |line|
			  lines << line
			end
			f.close

			render json: lines
		else
			render text: ''
		end
	end

	def batch_ids
		ids = Array.new
		filename = 'nationwidechildrens.org_clinical_patient_' + params[:tumor_type] + '.txt'
		path = '../data/' # will need to alter to source of data pulled by file downloader
		f = File.open(path + filename, "r")
		numLines = 0
		f.each_line do |line|
			if (numLines > 2 && line != '')
				fullID = line.strip.split(' ')[1]
				ids.insert(-1, fullID.split('-')[1])
			end
			numLines = numLines + 1
		end
		f.close
		render json: {ids: ids}
	end

	def get_prediction_targets
		path = '../data/nationwidechildrens.org_clinical_patient_' + params[:tumor_type] + '.txt'

		targets = []
		i = 0
		File.open(path, "r") do |file_handle|
			file_handle.each_line do |line|
				if i == 1
					targets = line.strip.split("\t")
				elsif i > 1
					break
				end
				i += 1
			end
		end

		path = '../data/not_prediction_target.txt'
		non_targets = []
		File.open(path, "r") do |file_handle|
			file_handle.each_line do |line|
				non_targets << line.strip
			end
		end

		data = []
		targets.each do |target|
			if !non_targets.include?(target)
				data << target
			end
		end

		data = data.sort_by!{ |var| var.upcase }
		# Always place "gender" first, as it is a good introductory example.
		data.delete("gender")
		data.unshift("gender")

		render json: data
	end

	def get_clinical_variables
		path = '../data/nationwidechildrens.org_clinical_patient_' + params[:tumor_type] + '.txt'

		clinical_variables = []
		i = 0
		File.open(path, "r") do |file_handle|
			file_handle.each_line do |line|
				if i == 1
					clinical_variables = line.strip.split("\t")
				elsif i > 1
					break
				end
				i += 1
			end
		end

		render json: clinical_variables
	end

	def get_unique_prediction_target_values
		path = '../data/nationwidechildrens.org_clinical_patient_' + params[:tumor_type] + '.txt'

		column_num = 0
		data = []
		i = 0
		File.open(path, "r") do |file_handle|
			file_handle.each_line do |line|
				if i == 1
					clinical_variables = line.strip.split("\t")
					column_num = clinical_variables.index(params[:prediction_target])
				elsif i > 2
					clinical_variables = line.strip.split("\t")
					data << clinical_variables[column_num]
				end
				i += 1
			end
		end
		data = data.sort_by!{ |var| var.upcase }
		render json: data.uniq
	end

	private #------------------------------------------------------------------------------#

	def build_command(validators, params)
		command = ""
		valid_command = true
		validators.each do |param, possible_values|
			arg = params[param]
			arg = Integer(arg) if possible_values == 'INTEGER'
			arg = Float(arg) if possible_values == 'FLOAT'
			arg = -1 if arg == -1.0

			command += " #{arg}"
		end

		{command: command, valid_command: valid_command}
	end

	def get_aucs(session_id)
		aucs = File.open("public/sessions/#{session_id}/AUCs.txt", "r").read.split("\n")
		aucs.each_with_index do |auc, index|
			tmp = auc.split(",")
			aucs[index] = {name: tmp[0], value: tmp[1]}
		end
		return aucs
	end

	def get_feature_weights(session_id)
		feature_weights = File.open("public/sessions/#{session_id}/featureWeights.txt", "r").read.split("\n")
		feature_weights.each_with_index do |feature_weight, index|
			tmp = feature_weight.split(",")
			feature_weights[index] = {name: tmp[0], value: tmp[0].gsub(/-.-.\z/, '').gsub(/_.*\z/, ''), weight: tmp[1]}
		end
		return feature_weights
	end
end
