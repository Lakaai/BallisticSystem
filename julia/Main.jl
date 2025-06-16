using CSV
using Plots
using DataFrames
using LinearAlgebra
using Infiltrator

# TODO: Setup Unit Test for logpdf
# TODO: Setup toy problem for using Optim bfgs update method
# TODO: Gaussian Processes 

include("SystemEstimator.jl")
include("Gaussian.jl")

data = CSV.read("data/estimationdata.csv", DataFrame)


nx = 3  # Dimension of state vector
ny = 1  # Dimension of measurement vector

t_hist = data[:, 1][1:end-1] # Store time stamps
println(size(t_hist))
x_hist = [] # Store state vectors
y_hist = data[:, 6] # Store measurement data 

mu_hist = []
sigma_hist = []

mu0 = [14000.0; -450.0; 0.0005] # Initial state estimate
S0 = Matrix(Diagonal([2200.0, 100.0, 1e-3])) # Initial covariance estimate

# Create the Gaussian object 
density = Gaussian(mu0, S0)

nsteps = nrow(data) - 1 # Number of time steps
println("nsteps: ", nsteps)

for i = 1:nsteps

    time = t_hist[i]        # Get the current time
    println("time in main loop", typeof(time))
    y = y_hist[i]           # Get the current measurement

    # 1. Predict forward in time 
    predicted_density = predict(time, density) # Form the predicted density p(x[k] ∣ y[k]:y[k-1]) by propagating p(x[k-1] ∣ y[k]:y[k-1]) through the process model 
    println("predicted state density", predicted_density.mean)

    # 2. Process the measurement event
    posterior_density = measurement_update_bfgs(predicted_density, y) # Compute the filtered density p(x[k] ∣ y[1]:y[k])

    # 3. Store the data for plotting
    push!(mu_hist, posterior_density.mean)
    push!(sigma_hist, posterior_density.covariance)
    
end 

# Convert to matrices for easier plot handling
mu_matrix = hcat(mu_hist...)
sigma_matrix = hcat(sigma_hist...)
# println("mu_hist: ", mu_hist)
# println("sigma_hist: ",  sigma_hist)

println("mu_matrix: ", size(mu_matrix))
println("sigma_matrix: ", size(sigma_matrix))

# Plotting
println("size of t_hist: ", size(t_hist))
println("size of mu_matrix: ", size(mu_matrix))

# Example: mu_hist is a 3×N matrix, t_hist is a vector of length N
plot(t_hist, mu_matrix[1, :], label="h (Altitude)", xlabel="Time [s]", ylabel="States", lw=2)
plot!(t_hist, mu_matrix[2, :], label="v (Velocity)", lw=2)
plot!(t_hist, mu_matrix[3, :], label="c (Drag Coefficient)", lw=2)
