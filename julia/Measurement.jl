include("Gaussian.jl")

const r1 = 5000; # Horizontal position of sensor [m]
const r2 = 5000; # Vertical position of sensor [m]

"""
    predict_measurement(x; grad=false)

Evaluate the nonlinear measurement function `h(x)` and, optionally, its Jacobian.

# Arguments
- `x::AbstractVector`: The state vector, assumed to be at least of length 1.
- `grad::Bool`: If `true`, also returns the Jacobian ∂h/∂x.

# Returns
- If `grad=false`: Returns the predicted measurement `h`.
- If `grad=true`: Returns a tuple `(h, dhdx)` where:
    - `h::Float64`: Predicted measurement.
    - `dhdx::Matrix`: Jacobian of the measurement function with respect to state `x`.

# Notes
- Assumes a measurement model of the form `y = h(x) + v` with noise `v`.
- Uses fixed parameters `r1`, `r2` assumed to be defined in the outer scope.
"""
function predict_measurement(x; grad=false)

    h1 = x[1]
    h = sqrt(r1^2 + (h1-r2)^2)

    if grad 
        dhdx = zeros(eltype(x), 1, 3) # Create matrix filled with zeros, where each zero has the same type as the elements of x (required for handling Dual types).
        dhdx[1, 1] = (h1 - r2)/h
        return h, dhdx
    else 
        return [h]
    end 
end 

"""
    predictDensity(x)

Computes the predicted measurement density `p(y | x)` as a Gaussian and returns its Jacobian.

# Arguments
- `x::AbstractVector`: The current state estimate.

# Returns
- `Gaussian([h], SR)`: A Gaussian distribution representing the predicted measurement.
- `dhdx::Matrix`: The Jacobian ∂h/∂x of the measurement function at `x`.

# Notes
- Assumes a fixed standard deviation (`sigma_rng = 50.0`) for the measurement noise.
- Measurement model: `y = h(x) + v` with v ∼ N(0, σ²).
"""
function predictDensity(x)
    sigma_rng = 50.0 # 50m standard deviation

    h, dhdx = predict_measurement(x; grad=true)
    
    SR = Matrix(Diagonal([sigma_rng]))
    # TODO: Fix (see lab5 cpp)

    return Gaussian([h], SR), dhdx
end 

"""
    logLikelihood(x, measurement; grad=false)

Computes the log-likelihood `log p(y | x)` of a measurement under the predicted measurement model,
optionally returning the gradient with respect to state `x`.

# Arguments
- `x::AbstractVector`: The current state estimate.
- `measurement::AbstractVector`: The actual observed measurement.
- `grad::Bool`: If `true`, also returns the gradient ∂/∂x log p(y | x).

# Returns
- If `grad=false`: Returns `log_likelihood::Float64`.
- If `grad=true`: Returns a tuple `(log_likelihood, gradient)` where:
    - `gradient::Vector`: Gradient of the log-likelihood with respect to state.

# Notes
- Assumes Gaussian measurement noise.
- Uses the chain rule: ∂/∂x log p(y | x) = -Jᵀ ∂/∂y log N(y; h(x), R)
- Evaluates log N(y; h(x), R) and d/dy log N(y; h(x), R)
"""
function logLikelihood(x, measurement; grad=false)
    
    likelihood, dhdx = predictDensity(x) # Compute the measurement likelihood p(x∣z) = N(y; h(x), R)
    log_likelihood, log_likelihood_gradient =  log_sqrt_pdf(measurement, likelihood; grad=true) # Compute the log likelihood logp(x∣z) = log N(y; h(x), R)
    if grad
        gradient = -dhdx' * log_likelihood_gradient;
        return log_likelihood, gradient;
    else 
        return log_likelihood;
    end 
end 