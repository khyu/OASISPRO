class OmicsController < ApplicationController
	require 'site_constants.rb'

	def index
		if params[:generate]
			tumor_type = params[:tumor_type]
			data_source = params[:data_source]
			gene_name = params[:gene_name]
			clinical_variable = params[:clinical_variable]

			session_id = rand.to_s.sub("0.", "")
			@command = "Rscript public/RCodes/omicsVisualization.R #{tumor_type} #{data_source} #{gene_name} #{clinical_variable} #{session_id}"
			system(@command)
		end
	end

	def get_gene_names
		tumor_type = params[:tumor_type]
		data_source = params[:data_source]
		if SiteConstants::TUMOR_TYPES[tumor_type.to_sym] && SiteConstants::DATA_TYPES[data_source.to_sym]
			gene_names = File.open("../data/#{tumor_type}_#{data_source}_elemid.txt", "r").read
			gene_names = gene_names.split(',')
			render json: gene_names
		end
	end
end