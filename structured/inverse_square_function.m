function B = inverse_square_function(dim)

    switch dim
    case 12
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
    B(i, j, k, l, m, n, o, p, q, r, s, t) = (i^2 + j^2 + k^2 + l^2 + m^2 + n^2 + o^2 + p^2 + q^2 + r^2 + s^2 + t^2)^(-1/2);
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

    case 6
    for i=1:16
    for j=1:16
    for k=1:16
    for l=1:16
    for m=1:16
    for n=1:16
    B(i, j, k, l, m, n) = (i^2 + j^2 + k^2 + l^2 + m^2 + n^2)^(-1/2);
    end
    end
    end
    end
    end
    end

    end

end
