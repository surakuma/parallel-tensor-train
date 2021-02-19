function [G INVS] = helper_parallelTensorTrainCompression_h3_gFixedRankSVD(r1, A, r2, dims, dims_identity, fixed_rank)
    
    if length(dims_identity) == 1
    G{dims_identity} = reshape(A, r1, numel(A)/(r1*r2), r2);
    INVS = {};
    return;
    end

    d = length(dims);
    d_half = idivide(d+1, int32(2));

    first_dim = r1 * prod(dims(1:d_half));

    C = reshape(A, first_dim, numel(A)/first_dim);
    tic
    [U S V] = svd(C, 'econ');
    toc
    r = findRank(S);
    r = min(r, fixed_rank);

    d1 = length(dims(1:d_half));
    d2 = length(dims(d_half+1:d));

    U=U(:,1:r);
    S=S(1:r,1:r);
    V=V(:,1:r);
    %INVS{dims_identity(d_half:d_half)} = 1./diag(S);
    %INVS = 1./diag(S);


    [G1 INVS1] = helper_parallelTensorTrainCompression_h3_gFixedRankSVD(r1, U*S, r, dims(1:d_half), dims_identity(1:d_half), fixed_rank);
    [G2 INVS2] = helper_parallelTensorTrainCompression_h3_gFixedRankSVD(r, S*V', r2, dims(d_half+1:d), dims_identity(d_half+1:d), fixed_rank);

    G = {G1{1:dims_identity(d_half)}, G2{dims_identity(d_half+1:d)}};
    INVS = {INVS1{1:length(INVS1)}, 1./(diag(S)'), INVS2{1:length(INVS2)}};
    
end


