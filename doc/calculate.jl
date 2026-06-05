using dn2
using Plots

# Calculate the force between the cubes
n = 10
force_n = dn2.force(10)
rounded_force_n = round(force_n, digits=10)
println("Calculated force between two unit cubes one unit apart rounded to 10 digits: $rounded_force_n")

force_n1 = dn2.force(11)
error = abs(force_n1 - force_n)
println("\nCalculated error: $error")

# Plot the error dropping as we increase n
println("Plotting error drop off as we increase the number of nodes...")
ns = 2:n
vals = Float64[]
for _n in ns
  push!(vals, dn2.force(_n))
end

plot(ns, abs.(vals .- force_n1), ylabel="Error", xlabel="Number of nodes", yscale=:log10, marker=:circle, legend=false, yticks=10.0 .^ (0:-5:-15))

