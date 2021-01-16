function [r] = findsuitableRank(S)
    [rows columns] = size(S);
    %r = min(rows(S), columns(S));
    r = min(rows, columns);
end


