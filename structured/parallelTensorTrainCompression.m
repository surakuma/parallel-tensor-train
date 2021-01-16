function [G] = parallelTensorTrainCompression(A, accur)
    ndims = size(A);
    accur_per_decomposition = accur/sqrt(length(ndims)-1);
    G = helper_parallelTensorTrainCompression(1, A, 1, ndims, 1:length(ndims), accur_per_decomposition);   
end


