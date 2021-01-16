function B = inverse_penta_function(dim)

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
    B(i, j, k, l, m, n, o, p, q, r, s, t) = (i^5 + j^5 + k^5 + l^5 + m^5 + n^5 + o^5 + p^5 + q^5 + r^5 + s^5 + t^5)^(-1/5);
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
    B(i, j, k, l, m, n) = (i^5 + j^5 + k^5 + l^5 + m^5 + n^5)^(-1/5);
    end
    end
    end
    end
    end
    end

    end

end
