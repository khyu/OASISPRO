class AnalysisController < ApplicationController
	require 'site_constants.rb'

	def stage
		if params[:generate]
			validators = {
				tumor_type: SiteConstants::TUMOR_TYPES.keys.map { |x| x.to_s },
				data_source: SiteConstants::DATA_TYPES.keys.map { |x| x.to_s },
				prediction_target: SiteConstants::PREDICTION_TARGTS.keys.map { |x| x.to_s },
				partition: SiteConstants::PARTITION_TYPES,
				random_seed: 'INTEGER',
				training_percentage: 'FLOAT',
				feature_selection_method: SiteConstants::FEATURE_SELECTION_METHOD.keys.map { |x| x.to_s },
				num_top_features: 'INTEGER'
			}
			params[:random_seed] = -1 if params[:random_seed].blank?

			if params[:training_percentage].present?
				params[:training_percentage] = Float(params[:training_percentage])/100
			else
				params[:training_percentage] = -1
			end

			@valid_command = true
			@command = "Rscript binaryClassification.R"

			# Validate and build the command
			validators.each do |param, possible_values|
				arg = params[param]
				arg = Integer(arg) if possible_values == 'INTEGER'
				arg = Float(arg) if possible_values == 'FLOAT'
				arg = -1 if arg == -1.0

				if ((possible_values.class == Array && possible_values.include?(arg)) ||
					(possible_values.class != Array && arg.is_a?(Numeric)))
					@command += " #{arg}"
				else
					@valid_command = false
					break
				end
			end

			if @valid_command
				system(@command)
			end

			# TBD
			@area_under_curve = []
		end
	end

	def survival
		if params[:generate]
			validators = {
				tumor_type: SiteConstants::TUMOR_TYPES.keys.map { |x| x.to_s },
				data_source: SiteConstants::DATA_TYPES.keys.map { |x| x.to_s },
				alpha_lower_bound: 'FLOAT',
				alpha_upper_bound: 'FLOAT',
				lambda_lower_bound: 'FLOAT',
				lambda_upper_bound: 'FLOAT'
			}

			# Default for alpha/lambda bound params
			validators.keys[2..5].each do |param|
				params[param] = -1 if params[param].blank?
			end

			@valid_command = true
			@command = "Rscript elasticNetCox.R"

			validators.each do |param, possible_values|
				arg = params[param]
				arg = Integer(arg) if possible_values == 'INTEGER'
				arg = Float(arg) if possible_values == 'FLOAT'
				arg = -1 if arg == -1.0

				if ((possible_values.class == Array && possible_values.include?(arg)) ||
					(possible_values.class != Array && arg.is_a?(Numeric)))
					@command += " #{arg}"
				else
					@valid_command = false
					break
				end
			end
			
			if @valid_command
				system(@command)
			end
		end
	end

	def get_features
		source = params[:source]
		method = params[:method]

		lines = Array.new
		f = File.open("public/feature_selection/"+source+"-"+method+".txt", "r")
		f.each_line do |line|
		  lines << line
		end
		f.close

		render json: lines
	end
end
