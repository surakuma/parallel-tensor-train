function [X] = flatApproximationForTTDecomposition(G)

    %Not sure
    d = length(G);
    [dim1 dim2 dim3] = size(G{1});
    fprintf('Core tensor size = %d X %d X %d\n', dim1, dim2, dim3);
    X = reshape(G{1}, dim2, numel(G{1})/dim2);

    for k=2:d
    [dim1 dim2 dim3] = size(G{k});
    fprintf('Core tensor size = %d X %d X %d\n', dim1, dim2, dim3);
    temp = reshape(G{k}, dim1, dim2*dim3);
    X = X*temp;
    X = reshape(X, numel(X)/dim3, dim3);
    end
    X=X';
end


