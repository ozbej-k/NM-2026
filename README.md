# Force of gravity

Author: Ožbej Kresal

In this exercise we implemented a Julia package which calculates the gravitational force between two parallel-placed unit homogeneous cubes at a distance of 1 (between the closest sides). Where the force between two bodies $T_1,T_2 \subset \mathbb{R}^3$ is equal to

$$\mathbf{F} = \int_{T_1}\int_{T_2}\frac{\mathbf{r}_1-\mathbf{r}_2}{\left\|\mathbf{r}_1-\mathbf{r}_2\right\|_2^3}\, d\mathbf{r}_1 \, d\mathbf{r}_2.$$

To do that we used the Gauss-Legendre quadrature, where we implemented both the calculation of nodes and weights using the Golub-Welsch algorithm.

# Run instructions
Start the interactive Julia loop in the folder where this readme file is located.

Some examples using our implemented functions can be found in 
[calculate.jl](doc/calculate.jl), which can be run with
```jl
# ] to enter package mode
activate .
# RETURN to exit package mode
include("doc/calculate.jl")
```
in the interactive Julia loop.

## Tests
To run the tests for the custom Gauss-Legendre nodes and weights implementation run the following command in the interactive loop:
```julia
# ] to enter package mode
activate .
test dn2
```

## Report PDF
To create the PDF report run the following command:
```julia
# ] to enter package mode
activate .
# RETURN to exit package mode
include("doc/make_report.jl")
```

