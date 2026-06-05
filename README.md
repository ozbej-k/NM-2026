# SOR iteration for sparse matrices

Author: Ožbej Kresal

In this exercise we implemented a sparse matrix datatype RedkaMatrika and the necessary interfaces to use the datatype to solve a $Ax = b$ system using SOR iteration. We then tested our sparse matrix and SOR implementation on two graphs, where we used the physical method for graph embedding to embed them in a plane.

## Sparse matrix implementation

We defined a new data type RedkaMatrika, which stores the non-zero elements of the matrix due to space requirements in two double arrays V and I, such that:

$$V[i][j] = a_{i,\,I[i][j]}$$

## SOR iteration to solve sparse linear systems
To solve sparse systems $Ax = b$ we implemented a SOR iteration function adapted to work with our sparse matrix datatype:

$$x,\; it = \mathrm{sor}(A,b,x_0,\omega,\texttt{tol}=10^{-10}),$$

where $x_0$ is the initial approximation, $\texttt{tol}$ is the condition for stopping the iteration, and $\omega$ is the parameter in the SOR iteration. The iteration then stops when

$$\|Ax^{(k)}-b\|_{\infty} < \delta,$$

where $\delta$ is given by the $\texttt{tol}$ argument.

# Run instructions
Start the interactive Julia loop in the folder where this readme file is located.

Some examples using our implemented functions can be found in 
[graph_embed.jl](doc/graph_embed.jl), which can be run with
```jl
# ] to enter package mode
activate .
# RETURN to exit package mode
include("doc/graph_embed.jl")
```
in the interactive Julia loop.

## Tests
To run the tests run the following command in the interactive loop:
```julia
# ] to enter package mode
activate .
test dn1
```

## Report PDF
To create the PDF report run the following command:
```julia
# ] to enter package mode
activate .
# RETURN to exit package mode
include("doc/make_report.jl")
```

