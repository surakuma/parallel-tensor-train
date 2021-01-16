function [G] = helper_parallelTensorTrainCompression(r1, A, r2, dims, dims_identity, accur)
    
    if length(dims_identity) == 1
    %dims_identity
    G{dims_identity} = reshape(A, r1, numel(A)/(r1*r2), r2);
    return;
    end

    d = length(dims);
    d_half = idivide(d+1, int32(2));

    first_dim = r1 * prod(dims(1:d_half));

    C = reshape(A, first_dim, numel(A)/first_dim);
    [U S V] = svd(C, 'econ');

    r = findsuitableRankWithAccuracies(S, accur);

    U=U(:,1:r);
    S=S(1:r,1:r);
    V=V(:,1:r);

    G1 = helper_parallelTensorTrainCompression(r1, U, r, dims(1:d_half), dims_identity(1:d_half), accur);
    G2 = helper_parallelTensorTrainCompression(r, S*V', r2, dims(d_half+1:d), dims_identity(d_half+1:d), accur);
    %dims_identity(1:d_half)
    %dims_identity(d_half+1:d)
    G = {G1{1:dims_identity(d_half)}, G2{dims_identity(d_half+1:d)}};
    %dims_identity(1:d_half)
    %dims_identity(d_half+1:d)
    
end


