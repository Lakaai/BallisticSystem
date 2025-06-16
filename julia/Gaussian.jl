# using LinearAlgebra

""" 
    struct GaussianSqrt(mu, S) 
    Construct a Gaussian distribution with given mean and square root covariance matrix.

    # Arguments 
    - `μ` The mean vector of the Gaussian distribution.
    - `P` The covariance matrix of the Gaussian distribution.

"""
struct Gaussian{T}
    mean
    covariance::Matrix{T}
end

# Constructor with specified mean and covariance
function Gaussian(mean, covariance::Matrix{T}) where {T}
    @assert length(mean) == size(covariance, 1) "mean and covariance size mismatch"
    @assert size(covariance, 1) == size(covariance, 2) "covariance matrix must be square"
    Gaussian{T}(mean, Matrix(covariance))
end

""" 
    from_moment(μ, P) 
    Construct a Gaussian distribution with given mean and covariance matrix.

    # Arguments 
    - `μ` The mean vector of the Gaussian distribution.
    - `P` The covariance matrix of the Gaussian distribution.

"""
function from_moment(μ::Vector, P)
    return Gaussian(μ, Matrix(P))
end 

""" 
    from_sqrt_moment(μ, S) 
    Construct a Gaussian distribution with given mean and square root covariance matrix.

    # Arguments 
    - `μ` The mean vector of the Gaussian distribution.
    - `S` The square root covariance matrix (upper triangular) of the Gaussian distribution.

"""
function from_sqrt_moment(μ::Vector, S)
    @assert istriu(Matrix(S)) "S must be upper triangular"
    return Gaussian(μ, Matrix(S))
end 

""" 
    from_info(η, Λ) 
    Construct a Gaussian distribution with given information vector and information matrix.

    # Arguments 
    - `η` The information vector of the Gaussian distribution.
    - `Λ` The information matrix of the Gaussian distribution.

"""
function from_info(η::Vector, Λ)
    return Gaussian(η, Matrix(Λ))
end 

""" 
    from_sqrt_info(ν, Ξ) 
    Construct a Gaussian distribution with given information vector and square root information matrix.

    # Arguments 
    - `ν` The information vector of the Gaussian distribution.
    - `Ξ` The square root information matrix (upper triangular) of the Gaussian distribution.

"""
function from_sqrt_info(ν::Vector, Ξ)
    @assert istriu(Matrix(Ξ)) "Ξ must be upper triangular"
    return Gaussian(ν, Matrix(Ξ))
end 

"""
    log_sqrt_pdf(x, pdf) 

Compute the logarithm of a multivariate normal distribution in square-root form at the value `x`. 

# Arguments
- `x` The input vector at which to evaluate the log-likelihood.
- `pdf` A multivariate normal distribution with mean `μ` and square-root covariance matrix `S` such that SᵀS = P.

# Returns
- The log of the probability distribution function evaluated at `x`.
"""
function log_sqrt_pdf(x, pdf::Gaussian; grad=false)
    μ = pdf.mean
    S = pdf.covariance
    n = length(x)

    @assert istriu(S) "S is not upper triangular"
    @assert length(x) == length(μ) "Input x and mean μ must have same length"

    Δ = x .- μ # always use .- to support scalar/vector and AD types
    w = LowerTriangular(transpose(S)) \ Δ   
    logpdf = -(n/2)*log(2π)-sum(log.(abs.(diag(Matrix(S)))))-(1/2)*dot(w,w)

    if grad
        gradient = -UpperTriangular(S) \ w      # Gradient ∇logp = -S⁻¹ * w
        return logpdf, gradient
    else 
        return logpdf                           # Return log N(x; μ, S)
    end 
end 

"""
    log_pdf(x, pdf) 

Compute the logarithm of a multivariate normal distribution in standard form at the value `x`. 

# Arguments
- `x` The input vector at which to evaluate the log-likelihood.
- `pdf` A multivariate normal distribution with mean `μ` and covariance matrix `P`.

# Keyword Arguments
- `grad`: Whether to return the gradient of the log-likelihood.

# Returns
- The log of the probability distribution function evaluated at `x`.
- Optionally, the gradient ∇logp(x) if `grad=true`.
"""
# TODO: Implement logpdf to compute log of a pdf in standard form.
# function log_pdf(x, pdf::Gaussian; grad=false)
#     μ = pdf.mean
#     P = pdf.covariance
#     n = length(x)

#     @assert length(x) == length(μ) "Input x and mean μ must have same length"

#     Δ = x .- μ # always use .- to support scalar/vector and AD types
#     w = LowerTriangular(transpose(P)) \ Δ   
#     println(diag(S))
#     logpdf = -(n/2)*log(2*π)-(1/2)*sum(log.(abs.(diag(Matrix(S)))))-(1/2)*dot(w,w)

#     if grad
#         gradient = -UpperTriangular(S) \ w      # Gradient ∇logp = -S⁻¹ * w
#         return logpdf, gradient
#     else 
#         return logpdf                           # Return log N(x; μ, S)
#     end 
# end 
