function [G INVS] = parallelTensorTrainCompression_h3_g(A, accur)
    ndims = size(A);
    [G INVS] = helper_parallelTensorTrainCompression_h3_g(1, A, 1, ndims, 1:length(ndims), accur);   
    d_minus1 = length(INVS);
    %d_minus1 = length(ndims) -1;

    for k=1:d_minus1
    [dim1 dim2 dim3] = size(G{k});
    temp = reshape(G{k}, dim1*dim2, dim3) * diag(INVS{k});
    G{k} = reshape(temp, dim1, dim2, dim3);
    end


end


