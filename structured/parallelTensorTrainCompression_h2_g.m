function [G] = parallelTensorTrainCompression_h2_g(A, accur)
    ndims = size(A);
    G = helper_parallelTensorTrainCompression_h2_g(1, A, 1, ndims, 1:length(ndims), accur);   
end


