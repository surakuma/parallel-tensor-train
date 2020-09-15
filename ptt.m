n=4;
d=5;
ndims=ones(1,d)*n;
clear A
A=rand(ndims);

clear A
A=B;


%G = tensorTrainDecomposition(A);
%
%%approx_error = computeError(A, G);
%
allowed_error = 0.0001;
G = tensorTrainCompression(A, allowed_error);
%
approx_error_tt = computeError(A, G)
%
allowed_error = 0.0001;
G1 = parallelTensorTrainCompression_NPG1(A, allowed_error);
approx_error_ptt = computeError(A, G1)

function [r_trunc] = findsuitableRankWithAccuracies(S, accr)
    [rows columns] = size(S);
    %r = min(rows(S), columns(S));
    r = min(rows, columns);
    r_trunc = 1;
    accr_2 = accr * accr;

    trunc_error_2 = 0;

    for k=r:-1:1

    trunc_error_2 = trunc_error_2 + S(k,k)*S(k,k);

    if trunc_error_2 > accr_2
    r_trunc = k;
    break
    end

    end
end

function [r] = findsuitableRank(S)
    [rows columns] = size(S);
    %r = min(rows(S), columns(S));
    r = min(rows, columns);
end

function [G] = tensorTrainDecomposition(A)
    ndims = size(A);
    d = length(ndims);

    C=A;
    r=1;
    for k=1:d-1

    C= reshape(C, r*ndims(k), numel(C)/(r*ndims(k)));
    [U S V] = svd(C);
    %r2 = min(rows(S), columns(S));
    r2 = findsuitableRank(S);

    U=U(:,1:r2);
    S=S(1:r2,1:r2);
    V=V(:,1:r2);

    G{k} = reshape(U, r, ndims(k), r2);
    C = S*V';
    r=r2;
    end

    G{d} = reshape(C, r, ndims(d), 1);

end

function [G] = helper_parallelTensorTrainCompression_NPG1(r1, A, r2, dims, dims_identity, accur)
    
    if length(dims_identity) == 1
    %dims_identity
    G{dims_identity} = reshape(A, r1, numel(A)/(r1*r2), r2);
    return;
    end

    d = length(dims);
    d_half = idivide(d+1, int32(2));

    first_dim = r1 * prod(dims(1:d_half));

    C = reshape(A, first_dim, numel(A)/first_dim);
    [U S V] = svd(C);

    r = findsuitableRankWithAccuracies(S, accur);

    U=U(:,1:r);
    S=S(1:r,1:r);
    V=V(:,1:r);

    G1 = helper_parallelTensorTrainCompression_NPG1(r1, U, r, dims(1:d_half), dims_identity(1:d_half), accur);
    G2 = helper_parallelTensorTrainCompression_NPG1(r, S*V', r2, dims(d_half+1:d), dims_identity(d_half+1:d), accur);
    %dims_identity(1:d_half)
    %dims_identity(d_half+1:d)
    G = {G1{1:dims_identity(d_half)}, G2{dims_identity(d_half+1:d)}};
    %dims_identity(1:d_half)
    %dims_identity(d_half+1:d)
    
end

function [G] = parallelTensorTrainCompression_NPG1(A, accur)
    ndims = size(A);
    accur_per_decomposition = accur/sqrt(length(ndims)-1);
    G = helper_parallelTensorTrainCompression_NPG1(1, A, 1, ndims, 1:length(ndims), accur_per_decomposition);   
end

function [G] = tensorTrainCompression(A, accur)
    ndims = size(A);
    d = length(ndims);

    C=A;
    r=1;
    accur_per_decomp = accur/sqrt(d-1);
    for k=1:d-1

    C= reshape(C, r*ndims(k), numel(C)/(r*ndims(k)));
    [U S V] = svd(C);
    %r2 = min(rows(S), columns(S));
    %r2 = findsuitableRank(S);
    r2 = findsuitableRankWithAccuracies(S, accur_per_decomp);

    U=U(:,1:r2);
    S=S(1:r2,1:r2);
    V=V(:,1:r2);

    G{k} = reshape(U, r, ndims(k), r2);
    C = S*V';
    r=r2;
    end

    G{d} = reshape(C, r, ndims(d), 1);

end



function [X] = flatApproximationForTTDecomposition(G)

    %Not sure
    d = length(G);
    [dim1 dim2 dim3] = size(G{1});
    X = reshape(G{1}, dim2, numel(G{1})/dim2);

    for k=2:d
    [dim1 dim2 dim3] = size(G{k});
    temp = reshape(G{k}, dim1, dim2*dim3);
    X = X*temp;
    X = reshape(X, numel(X)/dim3, dim3);
    end
    X=X';
end

function [numel_cf] = number_of_elements(G)
    d = length(G);
    
    for k=1:d
    numel(G{k})
    end

    numel_cf = numel(G{1});
    for k=2:d
     numel_cf = numel_cf + numel(G{k});
    end
end

function [approx_error] = computeError(A, G)
    error_vector = reshape(A, 1, numel(A)) - flatApproximationForTTDecomposition(G);
    approx_error = norm(error_vector, "fro");
    compression_ratio_percentage = (1 - number_of_elements(G)/numel(A))* 100
end




