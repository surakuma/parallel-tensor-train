function [U S V] = qrcp_before_svd(A, r)
    [rows columns] = size(A);
    r = min(r, min(rows, columns));
    [Q R P] = qr(A, 0);

    [U S V] = svd(R(1:r,:), "ECO");

    for k=r:-1:1
    if S(k,k) > 0
    final_r = k;
    break;
    end
    end

    %% inverse of permutation matrix
    pt(P) = 1:length(P);

    U = Q(:,1:r)*U(:,1:final_r);
    S = S(1:final_r, 1:final_r);
    V = V(pt,1:final_r);

end

    
