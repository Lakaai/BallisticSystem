using Optim

# Define the objective function
f(v) = (v[1] - 3)^2 + (v[2] + 1)^2

# Initial guess
initial_x = [0.0, 0.0]

# Run BFGS optimization
result = optimize(f, initial_x, BFGS())

# Extract result
println("Optimal solution: ", Optim.minimizer(result))
println("Minimum value: ", Optim.minimum(result))
println("Converged: ", Optim.converged(result))

using ForwardDiff

result = optimize(f, initial_x, BFGS(); autodiff = :forward)

# Extract result
println("Optimal solution: ", Optim.minimizer(result))
println("Minimum value: ", Optim.minimum(result))
println("Converged: ", Optim.converged(result))