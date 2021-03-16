function [U S V] = randomized_svd(A, r)
    [rows columns] = size(A);

    p = 10;
    p = 2;
    l = r+p;
    l = min(l, min(rows, columns));

    
    %% random Sketch matrix
    RS = randn(columns, l);

    Y = A*RS;

    q=5;
    q=0;
    %% power scheme method to select Q
    %% for k=1:q:1
    %% Y = transpose(A)*Y;
    %% Y = A*Y;
    %% end

    %% subspace iteration method to select Q
    for k=1:q:1
    [Q ~] = qr(Y,0);
    [Q ~] = qr(transpose(A)*Q,0);
    Y = A*Q;
    end

    [Q ~] = qr(Y,0);

    B = transpose(Q) *A;

    [U S V] = svd(B, "ECO");
    
    U = Q*U;

    r = min(r, l);
    U=U(:,1:r);
    S=S(1:r,1:r);
    V=V(:,1:r);

end
