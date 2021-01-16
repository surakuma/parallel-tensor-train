function [numel_cf] = number_of_elements(G)
    d = length(G);
    ne = [];
    
    for k=1:d
    ne = [ne numel(G{k})];
%    numel(G{k})
    end
    ne

    numel_cf = numel(G{1});
    for k=2:d
     numel_cf = numel_cf + numel(G{k});
    end
    total_elements = numel_cf
end


