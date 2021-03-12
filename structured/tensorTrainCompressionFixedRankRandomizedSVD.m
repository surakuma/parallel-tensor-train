function [G] = tensorTrainCompressionFixedRankRandomizedSVD(A, rank)
    ndims = size(A);
    d = length(ndims);

    C=A;
    r=1;
    for k=1:d-1

    C= reshape(C, r*ndims(k), numel(C)/(r*ndims(k)));
    tic
    [U S V] = randomized_svd(C, rank);
    toc
    r2 = findRank(S);
    %%r2 = min(r2, rank);

    U=U(:,1:r2);
    S=S(1:r2,1:r2);
    V=V(:,1:r2);

    %%step=k
    %%rank = r2
    %%step_error = norm(C-U*S*V', "fro")
    %%accur_per_decomp

    G{k} = reshape(U, r, ndims(k), r2);
    C = S*V';
    r=r2;
    end

    G{d} = reshape(C, r, ndims(d), 1);

end



