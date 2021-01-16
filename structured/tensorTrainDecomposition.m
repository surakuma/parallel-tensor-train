function [G] = tensorTrainDecomposition(A)
    ndims = size(A);
    d = length(ndims);

    C=A;
    r=1;
    for k=1:d-1

    C= reshape(C, r*ndims(k), numel(C)/(r*ndims(k)));
    [U S V] = svd(C, 'econ');
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


