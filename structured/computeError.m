function [approx_error] = computeError(A, G)
    error_vector = reshape(A, 1, numel(A)) - flatApproximationForTTDecomposition(G);
    approx_error = norm(error_vector, "fro");
    compression_ratio_percentage = (1 - number_of_elements(G)/numel(A))* 100
end


