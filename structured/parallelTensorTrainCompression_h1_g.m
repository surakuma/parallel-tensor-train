function [G] = parallelTensorTrainCompression_h1_g(A, accur)
    ndims = size(A);
    G = helper_parallelTensorTrainCompression_h1_g(1, A, 1, ndims, 1:length(ndims), accur);   
end

