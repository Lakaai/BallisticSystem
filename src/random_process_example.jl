using Plots
using Distributions

include("Gaussian.jl")

v = 1.0           # Velocity (m/s)
Δt = 1.0          # Time step (s)
σ_u² = 0.1        # Motion noise variance
σ_z² = 0.5        # Measurement noise variance
n_steps = 20      # Number of time steps

# Initialize arrays
true_x = zeros(n_steps)        # True position
z = zeros(n_steps)             # Measurements
x_hat = zeros(n_steps)         # Estimated position
P = zeros(n_steps)             # State uncertainty
true_x[1] = 0.0                # Initial true position
x_hat[1] = 0.0                 # Initial estimated position
P[1] = 0.1                     # Initial uncertainty

# Simulate true motion and noisy measurements
motion_noise = Normal(0, sqrt(σ_u²))
measurement_noise = Normal(0, sqrt(σ_z²))

for t = 2:n_steps
    # True motion: x_t = x_{t-1} + v * Δt + noise
    true_x[t] = true_x[t-1] + v * Δt + rand(motion_noise)
    # Noisy measurement: z_t = x_t + noise
    z[t] = true_x[t] + rand(measurement_noise)
end

# Kalman Filter
for t = 2:n_steps
    # Prediction step
    x_hat_pred = x_hat[t-1] + v * Δt
    P_pred = P[t-1] + σ_u²
    # Update step
    K = P_pred / (P_pred + σ_z²)          # Kalman gain
    x_hat[t] = x_hat_pred + K * (z[t] - x_hat_pred)
    P[t] = (1 - K) * P_pred
    
end

# Plot results
plot(1:n_steps, true_x, label="True Position", linewidth=2, color=:blue)
plot!(1:n_steps, z, label="Measurements", marker=:circle, color=:red, linestyle=:dash)
plot!(1:n_steps, x_hat, label="Estimated Position", linewidth=2, color=:green)
xlabel!("Time Step")
ylabel!("Position (m)")
title!("1D Robot Localization with Kalman Filter")

# function process_model()

#     return xₜ = xₜ₋₁ + v * Δt + uₜ 
# end 

# function measurement_model()
#     σ_squared = 0.5
#     wₜ = from_moment(0, σ_squared)
#     zₜ = xₜ + wₜ
# end 