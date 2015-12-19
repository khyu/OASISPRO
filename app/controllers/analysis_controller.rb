class AnalysisController < ApplicationController
	def stage
		if params[:generate]
			method = params[:feature_selection_method]
			k = Integer(params[:k]).to_s
			
			if params[:data_source] == 'rnaseq'
				file = 'runRnaseq.R'
			else
				file = 'runProteomics.R'
			end
			
			if ['ig', 'gain', 'su', 'chi'].include?(method)
				system 'Rscript public/dropbox/analysis/StageOverall/'+file+' '+method+' '+k
			end
			
			@area_under_curve = Array.new
			f = File.open("public/dropbox/analysis/StageOverall/outputAUC.csv", "r")
			
			aoc_labels = [
				'Recursive Partitioning Trees',
				'Conditional Inference Trees (CITs)', 
				'Random Forest with CITs', 
				'Bagging',
				'SVMs with Gaussian Kernel',
				'SVMs with Linear Kernel',
				'SVMs with Polynomial Kernel',
				'SVMs with Sigmoid Kernel',
				'Decision Trees',
				'Random Forest',
				'Naive Bayes',
				'Naive Bayes with Laplace Smoothing'
			]
			x = 0
			f.each_line do |line|
				@area_under_curve << [line, aoc_labels[x]]
				x += 1
			end
			f.close
		end
	end

	def survival
		if params[:generate]
			k = Integer(params[:k]).to_s
			
			if params[:data_source] == 'rnaseq'
				file = 'runRNAseqSurvival.R'
			else
				file = 'runProteomicsSurvival.R'
			end
			
			system 'Rscript public/dropbox/analysis/Survival/'+file+' '+k
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
