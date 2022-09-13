require 'bigdecimal'

class Performancer
	def initialize()
		@control = {}
	end

	def load_control(ref_array)
		ref_array.each do |node1, node2|
			if node2 != '-'
				query = @control[node1]
				if query.nil?
					@control[node1] = [node2]
				else
					query << node2
				end
			end
		end
		return control
	end

	# Pandey 2007, Association Analysis-based Transformations for Protein Interaction Networks: A Function Prediction Case Study
	def get_pred_rec(predictions, cut_number = 100, top_number = 10000)
		performance = [] #cut, pred, rec
		preds, limits = load_prediction(predictions)
		cuts = get_cuts(limits, cut_number)
		cuts.each do |cut|
			prec, rec = pred_rec(preds, cut, top_number)
			performance << [cut, prec, rec]
		end
		return performance
	end

	def load_prediction(pairs_array)
		pred = {}
		min = nil
		max = nil
		pairs_array.each do |key, label, score|
			query = pred[key]
			if !min.nil? && !max.nil?
				min = score if score < min
				max = score if score > max
			else
				min = score; max = score
			end
			if query.nil?
				pred[key] = [[label], [score]]
			else
				query.first << label
				query.last << score
			end
		end
		return pred, [min, max]
	end

	def pred_rec(preds, cut, top)
		predicted_labels = 0 #m
		true_labels = 0 #n
		common_labels = 0 # k
		@control.each do |key, c_labels|
			true_labels += c_labels.length #n
			pred_info = preds[key]
			if !pred_info.nil?
				labels, scores = pred_info
				reliable_labels = get_reliable_labels(labels, scores, cut, top)
				predicted_labels += reliable_labels.length #m
				common_labels += (c_labels & reliable_labels).length #k
			end
		end
		#puts "cut: #{cut} trueL: #{true_labels} predL: #{predicted_labels} commL: #{common_labels}"
		prec = common_labels.to_f/predicted_labels
		rec = common_labels.to_f/true_labels
		prec = 0.0 if prec.nan?
		rec = 0.0 if rec.nan?
		return prec, rec
	end


	def get_cuts(limits, n_cuts)
		cuts = []
		range = (limits.last - limits.first).abs.fdiv(n_cuts)
		range = BigDecimal(range, 10)
		cut = limits.first
		(n_cuts + 1).times do |n|
			cuts << (cut + n * range).to_f
		end
		return cuts
	end

	def get_reliable_labels(labels, scores, cut, top)
		reliable_labels = []
		scores.each_with_index do |score, i|
			reliable_labels << [labels[i], score] if score >= cut
		end
		reliable_labels = reliable_labels.sort!{|l1,l2| l2.last <=> l1.last}[0..top-1].map{|pred| pred.first}
		return reliable_labels
	end
end