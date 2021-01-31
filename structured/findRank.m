function [rank] = findRank(S)
    [rows columns] = size(S);
    %r = min(rows(S), columns(S));
    r = min(rows, columns);
    rank = 0;

    for k=r:-1:1

    if S(k,k) > 0
    rank = k;
    break
    end

    end
end
