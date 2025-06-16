using Test
using Infiltrator
# using your module or file where the function is defined
# include("../src/your_code_file.jl") or using MyModule
# include("SystemEstimator.jl")
# include("Gaussian.jl")

# @testset "measurement_update_bfgs" begin
#     # Dummy Gaussian struct

#     # Dummy cost function and measurement
#     function cost_joint_density(x, density, measurement)
#         return sum((x .- measurement).^2)
#     end

#     # Dummy input
#     mu = [1.0, 2.0]
#     S = Matrix{Float64}(I, 2, 2)
#     density = Gaussian(mu, S)
#     measurement = [1.1, 2.2]

#     result = measurement_update_bfgs(g, measurement)
    
#    #  @infiltrate true  # will hit if test runs

#     @test isa(result, Gaussian)
# end

# Laplace approximation of a pdf
# Given a pdf p(x) where x is in the set X the laplace approximation at the maximum a posterio of p(x) is given by 
# p(x) ≈ N(x; μ, P) where μ = argmin V(x), P⁻¹ = d²V/dx² evaluted at x = μ and the cost function V(x) = -logp(x)
# 𝒩
#  (slide 10 gaus filtering 2 lecture) \


# Square root implementation of the cost function is implemented in slide 15 gausfilter 2 lecture

# BFGS is a quasi-newton method meaning that it only requires the gradient of the cost function and will approximate the hessian 
# So in quasi-Newton methos we start with an approximation to the Hessian and we update as we iterate