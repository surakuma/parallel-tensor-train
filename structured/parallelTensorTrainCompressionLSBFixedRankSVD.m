function [G] = parallelTensorTrainCompressionLSBFixedRankSVD(A, rank)
    G = parallelTensorTrainCompression_h3_gFixedRankSVD(A, rank);
end
