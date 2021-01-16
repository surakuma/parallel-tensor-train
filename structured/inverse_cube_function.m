function B = inverse_cube_function(dim)

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
    B(i, j, k, l, m, n, o, p, q, r, s, t) = (i^3 + j^3 + k^3 + l^3 + m^3 + n^3 + o^3 + p^3 + q^3 + r^3 + s^3 + t^3)^(-1/3);
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
    B(i, j, k, l, m, n) = (i^3 + j^3 + k^3 + l^3 + m^3 + n^3)^(-1/3);
    end
    end
    end
    end
    end
    end

    end

end
