function caller()
    str_tensor_choice = 'Enter the tensor choice:\n1) A tensor with 4^12 elements\n2) A tensor with 16^6 elements\n';
    fprintf(str_tensor_choice);
    n = input('');

    str_tensor_function = 'Enter the low rank function\n1) Log function\n2) Sin function\n3) Inverse square function\n4) Inverse cube function\n5) Inverse penta function\n';
    fprintf(str_tensor_function);
    f = input('');


    switch n
    case 1

    switch f
    case 1
    A = log_function(12);
    case 2
    A = sin_function(12);
    case 3
    A = inverse_square_function(12);
    case 4
    A = inverse_cube_function(12);
    case 5
    A = inverse_penta_function(12);
    otherwise
    disp('Invalid input, please input a valid option')
    return;
    end

    case 2

    switch f
    case 1
    A = log_function(6);
    case 2
    A = sin_function(6);
    case 3
    A = inverse_square_function(6);
    case 4
    A = inverse_cube_function(6);
    case 5
    A = inverse_penta_function(6);
    otherwise
    disp('Invalid input, please input a valid option')
    return;
    end


    otherwise
    disp('Invalid input, please input a valid option')
    return;
    end

    allowed_error = input('Enter the prescribed accuracy:');
    G = tensorTrainCompression(A, allowed_error);
    %G1ng = parallelTensorTrainCompression_h1_ng(A, allowed_error);
    G1 = parallelTensorTrainCompression_h1_g(A, allowed_error);
    %G2ng = parallelTensorTrainCompression_h2_ng(A, allowed_error);
    G2 = parallelTensorTrainCompression_h2_g(A, allowed_error);
    G3 = parallelTensorTrainCompression_h3_g(A, allowed_error);

    disp('Results for sequential tensor train approximation');
    approx_error_tt = computeError(A, G)
    %disp('Results for parallel tensor train approximation with LSR approach (no guarantee)');
    %approx_error_ptt_h1ng = computeError(A, G1ng)
    disp('Results for parallel tensor train approximation with LSR approach');
    approx_error_ptt_h1 = computeError(A, G1)
    %disp('Results for parallel tensor train approximation with SLSB approach (no guarantee)');
    %approx_error_ptt_h2ng = computeError(A, G2ng)
    disp('Results for parallel tensor train approximation with SLSB approach');
    approx_error_ptt_h2 = computeError(A, G2)
    disp('Results for parallel tensor train approximation with LSB approach');
    approx_error_ptt_h3 = computeError(A, G3)
    %size(A)

end 

