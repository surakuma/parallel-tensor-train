function [U S V] = randomized_svd(A, r)
    [rows columns] = size(A);

    p = 5;
    l = r+p;
    l = min(l, min(rows, columns));

    
    %% random Sketch matrix
    RS = randn(columns, l);

    Y = A*RS;

    q=2;
    for k=1:q:1
    Y = transpose(A)*Y;
    Y = A*Y;
    end

    [Q ~] = qr(Y,0);

    B = transpose(Q) *A;

    [U S V] = svd(B, "ECO");
    
    U = Q*U;

end
