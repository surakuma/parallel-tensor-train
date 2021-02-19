# parallel-tensor-train
This repository discusses parallel and scalable tensor train decomposition techiniques for high dimensional tensors. Presently parallel algorithms are implemented in Matlab. Distributed memory implementations of the parallel algorithms in C++ are currently underway. 


## How to run code of this repository?

## Requirements (Matlab code)
* Matlab/Octave
    
### Run:
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
   

## Requirements (C++ sequential code)
* LAPACK
* BLAS
