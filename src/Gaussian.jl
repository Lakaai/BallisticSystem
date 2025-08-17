# using LinearAlgebra

using QuadGK

""" 
    struct GaussianSqrt(mu, S) 
    Construct a Gaussian distribution with given mean and square root covariance matrix.

    # Arguments 
    - `μ` The mean vector of the Gaussian distribution.
    - `P` The covariance matrix of the Gaussian distribution.
"""
# Custom gaussian struct because Distributions.jl does not support square root covariance
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

"""
    condition(density::Gaussian, idx_𝑥::Int, idx_𝑦::Int, 𝑦)

Given the joint Gaussian N(μ, Σ) and index sets for variables A and B,
return 𝑝(𝑥 | 𝑦) for 

# Arguments
- `density::Gaussian`: The joint probability distribution function as a Gaussian density 𝑝(𝑥, 𝑦).
- `idx_𝑥 ::Vector
- `idx_𝑦 ::Vector
- `𝑦 
# Returns
- μ_cond :: Vector         — conditional mean μ_A|B
- Σ_cond :: Matrix         — conditional covariance Σ_A|B
"""
function condition(density, idx_𝑥, idx_𝑦, 𝑦; sqrt=false)

    if sqrt
        error("Not implemented yet") # TODO
        return from_sqrt_moment(μ_cond, Σ_cond)
    else
        μ = density.mean
        println("μ inside condition: ", μ)
        Σ = density.covariance
        
        n𝑥 = length(idx_𝑥)
        n𝑦 = length(idx_𝑦)
        idx_𝑦 = [idx_𝑦]
        # if !(𝑦 isa AbstractVector)
        #     @show 𝑦
        #     𝑦 = [𝑦]
        # end 

        μ𝑥 = μ[idx_𝑥]
        println("μ𝑥 inside condition: ", μ𝑥)
        μ𝑦 = μ[idx_𝑦]
        println("μ𝑦 inside condition: ", μ𝑦)
    
        Σ𝑥𝑥 = Σ[1:n𝑥, 1:n𝑥]
        Σ𝑥𝑦 = Σ[1:n𝑥, n𝑥+1:end]
        Σ𝑦𝑥 = Σ[n𝑥+1:end, 1:n𝑥]
        Σ𝑦𝑦 = Σ[n𝑥+1:end, n𝑥+1:end]

        # Compute the new mean and covariance of the conditional distribution 𝑝(𝑥 | 𝑦)
        # Dont invert that matrix (Σ𝑦𝑦⁻¹) - https://www.johndcook.com/blog/2010/01/19/dont-invert-that-matrix/

        # Instead, solve the linear system Σ𝑦𝑦 * w = v to find w = Σ𝑦𝑦⁻¹ * v 
        println("Σ𝑦𝑦 inside condition: ", Σ𝑦𝑦)
        # @show typeof(Σ𝑦𝑦)
        # @show 𝑦
        # @show μ𝑦
        w = Σ𝑦𝑦 \ (𝑦 - μ𝑦)

        # Compute the conditional mean μ𝑥|𝑦 = μ𝑥 + Σ𝑥𝑦 * Σ𝑦𝑦⁻¹ * (𝑦 - μ𝑦)
        μ_cond = μ𝑥 + Σ𝑥𝑦 * w  

        # Again solve the linear system Σ𝑦𝑦 * w = Σ𝑦𝑥 to find w = Σ𝑦𝑦⁻¹ * Σ𝑦𝑥
        w = Σ𝑦𝑦 \ Σ𝑦𝑥

        # Compute the conditional covariance Σ𝑥|𝑦 = Σ𝑥𝑥 - Σ𝑥𝑦 * Σ𝑦𝑦⁻¹ * Σ𝑦𝑥
        Σ_cond = Σ𝑥𝑥 - Σ𝑥𝑦 * w  

        # Return the conditional distribution 𝑝(𝑥 | 𝑦)
        return from_moment(μ_cond, Σ_cond)
    end 
end 

"""
    marginal(density::Gaussian, idx::Vector{Int})

Compute the marginal distribution of a Gaussian density over a subset of variables.

# Arguments
- `density::Gaussian`: The Gaussian density to marginalize.
- `idx::Vector{Int}`: The indices of the variables to marginalize over.
"""
function marginal(density::Gaussian, idx::Vector{Int})

    return from_moment(density.mean[idx], density.covariance[idx, idx])
end



"""
join(density_𝑥, density_𝑦; sqrt=false)

Construct the joint distribution of two independent Gaussian densities.

# Arguments
- `density_𝑥`: A Gaussian distribution representing the first random variable.
- `density_𝑦`: A Gaussian distribution representing the second random variable.

# Keyword Arguments
- `sqrt`: If `true`, constructs the joint in square-root form (not yet implemented). Defaults to `false`.

# Returns
- A new Gaussian representing the joint distribution, with concatenated means and a block-diagonal covariance matrix.
"""
function join(density_𝑥, density_𝑦; sqrt=false)

    μ = vcat(density_𝑥.mean, density_𝑦.mean)
   
    n𝑥 = size(density_𝑥.covariance, 1)
    n𝑦 = size(density_𝑦.covariance, 1)

    if sqrt
        error("Not implemented yet") # TODO
        return from_sqrt_moment(μ, S)
    else 
        # Create block diagonal matrix from the two covariance matrices
        Σ = zeros(n𝑥 + n𝑦, n𝑥 + n𝑦)
        Σ[1:n𝑥, 1:n𝑥] = density_𝑥.covariance
        Σ[n𝑥+1:end, n𝑥+1:end] = density_𝑦.covariance
        return from_moment(μ, Σ)
    end 
end


# Not necessarily Gaussian, but can be used to compute the sum of two independent random variables.
# """
#     sum(𝐗, 𝐘)

# Compute the sum of two independent random variables 𝐗 + 𝐘 by computing z = (f ∗ g), that is the convolution of the two probability density functions f and g.

# # Arguments
# - `f`: A probability density function representing the random variable 𝐗.
# - `g`: A probability density function representing the random variable 𝐘. 

# # Returns
# - `z`: A probability density function that represents 𝐙, that is the sum of the two random variables 𝐗 + 𝐘.

# """
# function sum(f::Function, g::Function)
    
#     # Compute the function z(s) = ∫ f(x) g(z - x) dx

#     function z(s)
        
#         # Define the integrand, that is the function to be integrated
#         integrand = (x) -> f(x) * g(s - x)

#         # Integrate the integrand
#         value, error = quadgk(integrand, -Inf, Inf)

#         return value
#     end 

#     return z
# end 

"""
    sum(X::Gaussian, Y::Gaussian) -> Gaussian

Return the Gaussian distribution corresponding to the sum of two independent
Gaussian random variables X and Y.

# Arguments
- `X::Gaussian`: First Gaussian distribution with mean μ and covariance Σ.
- `Y::Gaussian`: Second Gaussian distribution with mean μ and covariance Σ.

# Returns
- `Z::Gaussian`: Gaussian distribution of the sum Z = X + Y.
"""
# function sum(X::Gaussian, Y::Gaussian; sqrt=false)
#     if sqrt
#         error("Not implemented yet") # TODO
#         return from_sqrt_moment(μ, S)
#     else 
#         μ_Z = X.mean + Y.mean
#         Σ_Z = X.covariance + Y.covariance
#         return from_moment(μ_Z, Σ_Z)
#     end 
# end