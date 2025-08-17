# Van der Pol oscillator state space equations
# ̇x₁ = x₂
# ̇x₂ = μ(1-x₁²) * x₂ - x₁

# TODO: Implement side by side comparison of additive noise vs inject noise into rk4
# TODO: Implement state augmentation 
# TODO: Apparently this effect becomes increasingly obvious when there are severe nonlinearities, maybe try this in such a system and compare results

using Plots
using Random, LinearAlgebra, Distributions

include("system_estimator.jl")

# Φ: [xₖ; Δwₖ] ↦ xₖ₊₁

# Parameters
μ = 1.0                      # Nonlinearity/damping coefficient
dt = 0.01                    # Time step
T = 10.0                     # Total time
N = Int(T / dt)              # Number of steps
t = range(0, T, length=N+1)

σ = 0.2                      # Standard deviation of process noise
Q = σ^2 * dt * I(2)          # Process noise covariance (scaled by dt)
noise_dist = MvNormal([0.0, 0.0], Q)

# Initial conditions
x0 = [2.0, 0.0]
x_add = zeros(2, N+1)
x_inj = zeros(2, N+1)
x_add[:, 1] = x0
x_inj[:, 1] = x0

# Dynamics function
function vdp(x)
    x1, x2 = x
    dx1 = x2
    dx2 = μ * (1 - x1^2) * x2 - x1
    return [dx1, dx2]
end

# RK4 integrator (no noise)
function rk4_step(f, x, dt)
    k1 = f(x)
    k2 = f(x .+ 0.5 * dt * k1)
    k3 = f(x .+ 0.5 * dt * k2)
    k4 = f(x .+ dt * k3)
    return x .+ dt / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
end

# Simulation loop
for k in 1:N

    # Method 1: Additive noise after RK4
    x_det = rk4_step(vdp, x_add[:, k], dt)
    x_add[:, k+1] = x_det + rand(noise_dist)

    # Method 2: Noise injected into each RK4 stage
    w = rand(noise_dist)

    function vdp_noisy(x)
        return vdp(x .+ w)
    end
    
    x_inj[:, k+1] = rk4_step(vdp_noisy, x_inj[:, k], dt)
end

# Plotting
plot(t, x_add[1, :], label="x₁ (Additive Noise)", lw=2, color=:blue)
plot!(t, x_inj[1, :], label="x₁ (Injected Noise)", lw=2, color=:green, linestyle=:dash)
xlabel!("Time (s)")
ylabel!("x₁")
title!("Van der Pol Oscillator with Additive vs Injected Noise (RK4)")