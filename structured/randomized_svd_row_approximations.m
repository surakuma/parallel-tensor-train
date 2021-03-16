function [U S V] = randomized_svd_row_approximations(A, r)
    [rows columns] = size(A);

    p = 2;
    l = r+p;
    l = min(l, min(rows, columns));

    
    %% random Sketch matrix
    RS = randn(l, rows);

    Y = RS*A;

    [Q ~] = qr(Y',0);

    B = A*Q;

    [U S V] = svd(B, "ECO");
    
    V = Q*V;

    r = min(r, l);
    U=U(:,1:r);
    S=S(1:r,1:r);
    V=V(:,1:r);

end
