class OmicsController < ApplicationController
	require 'site_constants.rb'

	def index
		if params[:generate]
			tumor_type = params[:tumor_type]
			data_source = params[:data_source]
			gene_name = params[:gene_name]
			clinical_variable = params[:clinical_variable]

			@session_id = rand.to_s.sub("0.", "")

			Dir.mkdir("public/sessions") unless File.exists?("public/sessions")
			Dir.mkdir("public/sessions/#{@session_id}") unless File.exists?("public/sessions/#{@session_id}")

			@command = "Rscript public/RCodes/omicsVisualization.R #{tumor_type} #{data_source} #{gene_name} #{clinical_variable} #{@session_id}"
			system(@command)
			@tumor_type = tumor_type
			@data_source = data_source
			@gene_name = gene_name
			@clinical_variable = clinical_variable
			@pvalue = File.open("public/sessions/#{@session_id}/pValue.txt", "r").read.strip
		end
	end

	def filter_word_in(word,array)
    	array.delete_if { |data| data.match(/[0-9][0-9][0-9][0-9]/) }
    	array.delete_if { |data| data.to_s == '' }
    	array.delete_if { |data| data.match(word) }
    	return array
	end

	def get_gene_names
		tumor_type = params[:tumor_type]
		data_source = params[:data_source]
		if SiteConstants::TUMOR_TYPES[tumor_type.to_sym] && SiteConstants::DATA_TYPES[data_source.to_sym]
			gene_names = File.open("../data/#{tumor_type}_#{data_source}_elemid.txt", "r").read
			gene_names = gene_names.split(',')
			#gene_names = filter_word_in("100130426", gene_names)
			gene_names = filter_word_in("NA", gene_names)
			render json: gene_names
		end
	end
end