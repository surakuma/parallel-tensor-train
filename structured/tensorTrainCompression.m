function [G] = tensorTrainCompression(A, accur)
    ndims = size(A);
    d = length(ndims);

    C=A;
    r=1;
    accur_per_decomp = accur/sqrt(d-1);
    for k=1:d-1

    C= reshape(C, r*ndims(k), numel(C)/(r*ndims(k)));
    [U S V] = svd(C, 'econ');
    %r2 = min(rows(S), columns(S));
    %r2 = findsuitableRank(S);
    r2 = findsuitableRankWithAccuracies(S, accur_per_decomp);
    %eco_error = norm(C-U*S*V', "fro")

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



