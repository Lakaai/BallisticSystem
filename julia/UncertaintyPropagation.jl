# using LinearAlgebra
# include("Gaussian.jl")

# μ = [1;2;3;4]
# Σ = Matrix(Diagonal([1,1,1,1]))

# px = Gaussian(μ, Σ)

# """
#     logpdf(x, pdf) 

# Compute the logarithm of a multivariate normal distribution in square-root form at the value `x`. 

# # Arguments
# - x .
# - `pdf` A multivariate normal distribution with mean `μ` and square-root covariance matrix `S` such that SᵀS = P.

# # Returns
# - The log of the probability distribution function evaluated at `x`.
# """
# function logpdf(x::AbstractVector, pdf::Gaussian)
#     μ = pdf.mu 
#     S = pdf.S
#     n = length(x)

#     @assert istriu(S) "S is not upper triangular"

#     Δ = x - μ
#     w = LowerTriangular(transpose(S)) \ Δ   
    
#     logpdf = -(n/2)*log(2*π)-(1/2)*sum(log(abs(diag(S))))-(1/2)*dot(w,w)

#     return logpdf::Real # Return log N(x; μ, S)
# end 