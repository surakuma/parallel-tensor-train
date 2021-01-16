function [G] = parallelTensorTrainCompression_h1_ng(A, accur)
    ndims = size(A);
    d = length(ndims);
    accur_per_dim = accur / sqrt(d-1);
    G = helper_parallelTensorTrainCompression_h1_ng(1, A, 1, ndims, 1:length(ndims), accur_per_dim);   
end


