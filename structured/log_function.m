function B = log_function(dim)


    switch dim
    case 12
%    B = rand(3,5,2, 3, 6, 7, 8, 9);
%    return;

    for i=1:4
    for j=1:4
    for k=1:4
    for l=1:4
    for m=1:4
    for n=1:4
    for o=1:4
    for p=1:4
    for q=1:4
    for r=1:4
    for s=1:4
    for t=1:4
    B(i, j, k, l, m, n, o, p, q, r, s, t) = log (i + 2*j + 3*k + 4*l + 5*m + 6*n + 7*o + 8*p + 9*q + 10*r + 11*s + 12*t);
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    return;

    case 6
    for i=1:16
    for j=1:16
    for k=1:16
    for l=1:16
    for m=1:16
    for n=1:16
    B(i, j, k, l, m, n) = log (i + 2*j + 3*k + 4*l + 5*m + 6*n);
    end
    end
    end
    end
    end
    end
    return;

    end

end
