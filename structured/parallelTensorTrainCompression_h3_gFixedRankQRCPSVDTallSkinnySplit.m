function [G] = parallelTensorTrainCompression_h3_gFixedRankQRCPSVDTallSkinnySplit(A, rank, first_part_frac)
    ndims = size(A);
    [G INVS] = helper_parallelTensorTrainCompression_h3_gFixedRankQRCPSVDTallSkinnySplit(1, A, 1, ndims, 1:length(ndims), rank, first_part_frac);   
    d_minus1 = length(INVS);
    %d_minus1 = length(ndims) -1;

    for k=1:d_minus1
    [dim1 dim2 dim3] = size(G{k});
    temp = reshape(G{k}, dim1*dim2, dim3) * diag(INVS{k});
    G{k} = reshape(temp, dim1, dim2, dim3);
    end


end


