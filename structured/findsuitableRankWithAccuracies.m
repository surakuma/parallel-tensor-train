function [r_trunc] = findsuitableRankWithAccuracies(S, accr)
    [rows columns] = size(S);
    %r = min(rows(S), columns(S));
    r = min(rows, columns);
    r_trunc = r;
    accr_2 = accr * accr;

    trunc_error_2 = 0;

    for k=r:-1:1

    trunc_error_2 = trunc_error_2 + S(k,k)*S(k,k);

    if trunc_error_2 > accr_2
    r_trunc = k;
    break
    end

    end
end


