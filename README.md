# parallel-tensor-train
This repository discusses parallel and scalable tensor train decomposition techiniques for high dimensional tensors. Presently parallel algorithms are implemented in Matlab. Distributed memory implementations of the parallel algorithms in C++ are currently underway.

We proposed mainly 4 approaches to obtain tensor train approximations of high dimensional tensors. These are
1. STTA
2. LSR
3. SLSB
4. LSB


STTA is a sequential approach, while others are parallel. We refer the user to our [paper](https://hal.inria.fr/hal-03081555) for more details about these approaches.
 


## How to run code of this repository?

## Requirements (Matlab code)
* Matlab/Octave

### Enter into the structured folder (all matlab files are in this directory)
    >> cd structured

### Rank based examples

    >> A = rand(4,4,4,4);
    >> r=2

    >> Gsvd = tensorTrainCompressionFixedRankSVD(A, r);
    >> approx_error = computeError(A, Gsvd)
    Core tensor size = 1 X 4 X 2
    Core tensor size = 2 X 4 X 2
    Core tensor size = 2 X 4 X 2
    Core tensor size = 2 X 4 X 1
    ne =

        8   16   16    8

    total_elements =  48
    compression_ratio_percentage =  81.250
    approx_error =  4.2098

    >> Gqrcpsvd = tensorTrainCompressionFixedRankQRCPSVD(A, r);
    >> computeError(A, Gqrcpsvd)
    Core tensor size = 1 X 4 X 2
    Core tensor size = 2 X 4 X 2
    Core tensor size = 2 X 4 X 2
    Core tensor size = 2 X 4 X 1
    ne =

        8   16   16    8

    total_elements =  48
    compression_ratio_percentage =  81.250
    ans =  4.4020

    >> G3qrsvd = parallelTensorTrainCompressionLSBFixedRankQRCPSVD(A, r);
    >> approx_error = computeError(A, G3qrsvd)
    Core tensor size = 1 X 4 X 2
    Core tensor size = 2 X 4 X 2
    Core tensor size = 2 X 4 X 2
    Core tensor size = 2 X 4 X 1
    ne =

        8   16   16    8

    total_elements =  48
    compression_ratio_percentage =  81.250
    approx_error =  4.6781

### Long run  example(takes a long time):
     >> caller()
    Enter the tensor choice:
    1) A tensor with 4^12 elements
    2) A tensor with 16^6 elements
    1
    Enter the low rank function
    1) Log function
    2) Sin function
    3) Inverse square function
    4) Inverse cube function
    5) Inverse penta function
    1
    Enter the prescribed accuracy:0.001
    Results for sequential tensor train approximation
    Core tensor size = 1 X 4 X 3
    Core tensor size = 3 X 4 X 3
    Core tensor size = 3 X 4 X 3
    Core tensor size = 3 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 1
    ne =

       12   36   36   48   64   64   64   64   64   64   64   16

    total_elements =  596
    compression_ratio_percentage =  99.996
    approx_error_tt =    3.6922e-04
    
    Results for parallel tensor train approximation with LSR approach
    Core tensor size = 1 X 4 X 4
    Core tensor size = 4 X 4 X 12
    Core tensor size = 12 X 4 X 6
    Core tensor size = 6 X 4 X 15
    Core tensor size = 15 X 4 X 7
    Core tensor size = 7 X 4 X 4
    Core tensor size = 4 X 4 X 16
    Core tensor size = 16 X 4 X 22
    Core tensor size = 22 X 4 X 7
    Core tensor size = 7 X 4 X 15
    Core tensor size = 15 X 4 X 4
    Core tensor size = 4 X 4 X 1
    ne =

         16    192    288    360    420    112    256   1408    616    420    240     16

    total_elements =  4344
    compression_ratio_percentage =  99.974
    approx_error_ptt_h1 =    7.4966e-05
    
    Results for parallel tensor train approximation with SLSB approach
    Core tensor size = 1 X 4 X 4
    Core tensor size = 4 X 4 X 6
    Core tensor size = 6 X 4 X 5
    Core tensor size = 5 X 4 X 11
    Core tensor size = 11 X 4 X 7
    Core tensor size = 7 X 4 X 4
    Core tensor size = 4 X 4 X 9
    Core tensor size = 9 X 4 X 8
    Core tensor size = 8 X 4 X 6
    Core tensor size = 6 X 4 X 8
    Core tensor size = 8 X 4 X 4
    Core tensor size = 4 X 4 X 1
    ne =

        16    96   120   220   308   112   144   288   192   192   128    16

    total_elements =  1832
    compression_ratio_percentage =  99.989
    approx_error_ptt_h2 =    7.0418e-05
    
    Results for parallel tensor train approximation with LSB approach
    Core tensor size = 1 X 4 X 3
    Core tensor size = 3 X 4 X 3
    Core tensor size = 3 X 4 X 3
    Core tensor size = 3 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 4
    Core tensor size = 4 X 4 X 1
    ne =

       12   36   36   48   64   64   64   64   64   64   64   16

    total_elements =  596
    compression_ratio_percentage =  99.996
    approx_error_ptt_h3 =    3.6922e-04


   

/***** Implementations based on lapack and scalapack calls are in progress ******/

## Requirements (C++ sequential code)
* LAPACK
* BLAS

## Requirements (C++ parallel code)
* Scalapack
