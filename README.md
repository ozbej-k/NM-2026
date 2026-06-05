# Mathematical Pendulum

Author: Ožbej Kresal

In this exercise we implemented the Dormand-Prince 5 (DOPRI5) method for solving ordinary differential equations and applied it to the mathematical pendulum. We first transformed the second-order pendulum equation into an equivalent first-order system and then used our implementation to compute the pendulum's motion for different initial conditions. Finally, we compared the obtained solutions with the harmonic pendulum approximation and investigated the dependence of the oscillation period on the pendulum's energy.

# Run instructions
Start the interactive Julia loop in the folder where this readme file is located.

Some examples using our implemented functions can be found in 
[plot_stuff.jl](doc/plot_stuff.jl), which can be run with
```jl
# ] to enter package mode
activate .
# RETURN to exit package mode
include("doc/plot_stuff.jl")
```
in the interactive Julia loop.

## Tests
To run the tests run the following command in the interactive loop:
```julia
# ] to enter package mode
activate .
test dn3
```

## Report PDF
To create the PDF report run the following command:
```julia
# ] to enter package mode
activate .
# RETURN to exit package mode
include("doc/make_report.jl")
```

