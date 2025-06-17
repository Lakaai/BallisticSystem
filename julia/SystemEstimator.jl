using LinearAlgebra
using ForwardDiff
using Optim
using Infiltrator

include("Measurement.jl")

@enum UpdateMethod AFFINE UNSCENTED NEWTONTRUSTEIG BFGSTRUSTSQRTINV

const p0 = 101.325e3       # Air pressure at sea level [Pa]
const M  = 0.0289644       # Molar mass of dry air [kg/mol]
const R  = 8.31447         # Gas constant [J/(mol·K)]
const L  = 0.0065          # Temperature gradient [K/m]
const T0 = 288.15          # Temperature at sea level [K]
const g  = 9.81            # Acceleration due to gravity [m/s²]
time_ = 0                  # Previous time stamp 

# Laplace Aproximation ℐ = ∫ₓf(x)dx ≈ f(x*) √det(2*π*P)
# x* = argmax f(x)
function predict(time::Any, density::Gaussian)
    global time_
    dt = time - time_   # Update the time increment
    time_= time         # Update the previous time
    @assert dt >= 0 "dt must be non-negative"
    if dt == 0  
        return density
    end 
    # Define process noise covariance Q (power spectral density)
    Q = Matrix(Diagonal([1e-20, 25e-12, 0.0])) # [velocity, drag coeff, altitude (no noise)]
    predicted_density = affineTransform(rk4_step, density, dt, Q) # TODO: See GROK 
    return predicted_density

    # TODO: Implement square root covariance

    """
    # 1.1 Augment state density with independent noise increment dw ~ N^{-1}(0, LambdaQ/dt)
    # [ x] ~ N^{-1}([ eta ], [ Lambda,          0 ])
    # [dw]         ([   0 ]  [      0, LambdaQ/dt ])
    # pdw = processNoiseDensity(dt) (cpp)

    # SQ is an upper triangular matrix such that SQ.'*SQ = Q is the power spectral density of the continuous time process noise
    SQ = Matrix(Diagonal([sqrt(1e-20), sqrt(25e-12)])) # Set the non-zero elements of SQ based on the given Q matrix

    # SQ(0, 0) = std::sqrt(1e-20);  // For velocity
    # SQ(1, 1) = std::sqrt(25e-12); // For drag coefficient

    # Distribution of noise increment dw ~ N(0, Q*dt) for time increment dt (cpp) Gaussian<double>::fromSqrtMoment(SQ*std::sqrt(dt));
    #  pdw = Gaussian(0, SQ*sqrt(dt)) # p(dw(idxQ)[k])

    #  pxdw = density * pdw # p(x[k], dw(idxQ)[k]) = p(x[k])*p(dw(idxQ)[k])
    # Combine x and dw: joint distribution
    # TODO pxdw = density * pdw  # p([x; dw]) = p(x) * p(dw)

    # // Map p(x[k], dw(idxQ)[k]) to p(x[k+1])
    # density = pxdw.affineTransform(func);
    """
end 

# Affine transform of Gaussian: y = f(x), linearised around μ
function affineTransform(process_model::Any, density::Gaussian, dt, Q)

    μx = density.mean
    Σx = density.covariance

    J = ForwardDiff.jacobian(x -> process_model(x, dt), μx)
    
    μy = process_model(μx, dt)
    Σy = J * Σx * J' + Q

    @assert isapprox(Σy, Σy', rtol=1e-6) "Covariance not symmetric"
    @assert all(eigvals(Σy) .>= -1e-10) "Covariance not positive semi-definite"

    return from_moment(μy, Σy)

    # TODO: Implement square root covariance

    # fx = func(μx, dt)
    # J = ForwardDiff.jacobian(x -> func(x, dt), μx)
    # μy = fx
    # Sy = J * Sx  # propagate covariance via Jacobian
    

    # Ensure Sy is upper triangular via QR decomposition
    # Q, R = qr(Sy)
    # Sy_upper = R 

    # return from_sqrt_moment(μy, Sy_upper) 

end 

    
# Map [x[k]; dw(idxQ)[k]] to x[k+1] using RK4
function rk4_step(x::Any, dt::Any)
    k1 = dynamics(x)
    k2 = dynamics(x .+ 0.5 .* dt .* k1)
    k3 = dynamics(x .+ 0.5 .* dt .* k2)
    k4 = dynamics(x .+ dt .* k3)

    # RK4 step for deterministic part
    x_drift = x .+ dt/6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)

    # Add noise increment (already scaled as dw ~ N(0, Q*dt))
    # x_next = x_drift .+ dw
    x_next = x_drift

    return x_next
end 
    
# Evaluate f(x) from the SDE dx = f(x)*dt + dw
function dynamics(x::Any)

    # Extract state variables
    h = x[1]
    v = x[2]
    c = x[3]
    # @show x
    f = similar(x)

    # Calculate temperature at altitude h
    T = T0 - L * h
    
    d = ((0.5 * M * p0) / R) * (1 / T) * (1 - L * h / T0)^(g * M / (R * L)) * v^2 * c

    # Set f according to the dynamics equations
    f[1] = v;                 # dh/dt = v
    f[2] = d - g;             # dv/dt = d - g
    f[3] = 0;                 # dc/dt = 0 (drag coefficient is constant)

    return f
end 

function dynamics(x::Vector{Float64}, J::Matrix{Float64}) # Here J is passed by reference 
    f = dynamics(x)  # call the single-argument dynamics function

     # Extract state variables
    h = x[1]
    v = x[2]
    c = x[3]

    # Calculate temperature at altitude h
    T = T0 - L * h

    # Calculate drag acceleration
    d = ((0.5 * M * p0) / R) * (1 / T) * (1 - L * h / T0)^(g * M / (R * L)) * v^2 * c

    # Calculate partial derivatives
    dd_dh = ((0.5 * M * p0) / R) * v^2 * c * (
        (L / (T^2)) * (1 - L * h / T0)^(g * M / (R * L)) -
        (g * M / (R * T0)) * (1 - L * h / T0)^((g * M / (R * L)) - 1) / T
    )
    dd_dv = 2 * d / v
    dd_dc = d / c

    # Resize J to the correct size and fill in values
    J = zeros(length(f), length(x))

    J[1, 2] = 1.0
    J[2, 1] = dd_dh
    J[2, 2] = dd_dv
    J[2, 3] = dd_dc
    # J[3, :] stays zero

    return f
end

# Evaluate F(X) from dX = F(X)*dt + dW
function augmentedDynamics(X::Any)
    @assert size(X, 1) > 0 "X must have at least one row"
    nx = size(X, 1)

    x = X[:, 1]
    f, J = dynamics(x, J)

    @assert size(f, 1) == nx "f must have nx rows"
    @assert size(J, 1) == nx && size(J, 2) == nx "J must be nx by nx"

    dX = hcat(f, J * X[:, 2:end])
    return dX
end 

"""
    cost_function_factory(density::Gaussian, measurement)

The mean and covariance values must be available at the time of evaluation of the cost function. However, the cost function 
passed to the optimiser requires a function signature of f(x), therefore we will create a function factory to create a cost function 
with the required signature whilst still having access to the mean and covaraince values. 

In general a `factory` is an object for creating other objects, in this case it is a function that returns another function.

# Arguments
- `density` The system density.
- `measurement`.

# Returns
- A function with signature f(x) that can be passed to the optimiser.

"""
function cost_function_factory(density::Gaussian, measurement)
    return function(x) # Returns a cost function f(x) which has the required signature for the optimiser
        logprior = log_sqrt_pdf(x, density)

        # You must define this based on your measurement model
        loglik, _ = logLikelihood(x, measurement; grad=true)

         # Return a scalar cost (negative log-likelihood (the measurement cost) + log-prior (the prediction cost))
        return -(logprior + loglik) # Return −logp(x∣z) = -(logp(x) + logp(z∣x))
    end
end


function measurement_update_bfgs(density::Gaussian, measurement::Any)
    x0 = density.mean
    S = density.covariance

    df = TwiceDifferentiable(cost_function_factory(density, measurement), x0, autodiff = :forward) # Store and reuse gradient and hessian 
    res = optimize(df, x0, BFGS())
    # @assert res.converged "res has not converged"

    x_map = Optim.minimizer(res)

    # Posterior sqrt covariance approximation (naive)
    H = ForwardDiff.hessian(cost_function_factory(density, measurement), x_map)
    @show H
    # H_inv = inv(H)
    # Q, R = qr(H_inv) # Perform QR decomposition
    F = cholesky(H)   # H = F'U F, F.U is upper triangular
    S = Matrix(inv(F.U))      # S * S' = H^{-1}

    # S = R
    
    return Gaussian(x_map, S)
    # ℐ = ∫ₓf(x)dx ≈ f(x*) √det(2*π*P)
end 

function measurement_update_unscented(density::Gaussian, measurement::Any)
    x0 = density.mean
    S = density.covariance
    return Gaussian(x0, S)
end 

function measurement_update_affine(density::Gaussian, measurement::Any)
    μ_pred = density.mean
    Σ_pred = density.covariance
 
    H = ForwardDiff.jacobian(predict_measurement, μ_pred)
    R = Matrix(Diagonal([50.0])) # Measurement noise covariance

    K = Σ_pred * H' / (H * Σ_pred * H' + R)
    μ = μ_pred + K * ([measurement] - predict_measurement(μ_pred))
    Σ = (I - K * H) * Σ_pred
    
    # return from_sqrt_moment(μ, S) # TODO: Implement square root covariance
    return from_moment(μ, Σ)
end 

function update(density::Gaussian, measurement::Any, update_method::UpdateMethod)
    if update_method == BFGSTRUSTSQRTINV
        density = measurement_update_bfgs(density, measurement)
    elseif update_method == UNSCENTED
        error("Not implemented yet")
        # density = measurement_update_unscented(density, measurement) 
    elseif update_method == AFFINE
        density = measurement_update_affine(density, measurement)
    elseif update_method == NEWTONTRUSTEIG
        error("Not implemented yet")
        # density = measurement_update_newtontrusteig(density, measurement)
    else
        error("Invalid update method: $update_method")
    end  
end  