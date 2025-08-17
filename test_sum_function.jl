# Test file for the sum function in Gaussian.jl
# Run this file to test your sum function

# Include your Gaussian module
using LinearAlgebra
include("src/Gaussian.jl")

println("Testing sum function for Gaussian distributions...")
println("=" ^ 50)

# Test 1: Simple 1D case
println("\nTest 1: 1D Gaussian sum")
println("-" ^ 30)

# Create two 1D Gaussians
X1 = from_moment([1.0], reshape([2.0], 1, 1))  # mean=1, variance=2 (1x1 matrix)
Y1 = from_moment([3.0], reshape([1.5], 1, 1))  # mean=3, variance=1.5 (1x1 matrix)

println("X1: mean = $(X1.mean), variance = $(X1.covariance)")
println("Y1: mean = $(Y1.mean), variance = $(Y1.covariance)")

# Compute sum
Z1 = sum(X1, Y1)
println("Z1 = X1 + Y1: mean = $(Z1.mean), variance = $(Z1.covariance)")

# Expected result: mean = 1 + 3 = 4, variance = 2 + 1.5 = 3.5
expected_mean = 4.0
expected_variance = 3.5
println("Expected: mean = $expected_mean, variance = $expected_variance")
println("Test passed: $(abs(Z1.mean[1] - expected_mean) < 1e-10 && abs(Z1.covariance[1,1] - expected_variance) < 1e-10)")

# Test 2: 2D case
println("\nTest 2: 2D Gaussian sum")
println("-" ^ 30)

# Create two 2D Gaussians
X2 = from_moment([1.0, 2.0], [2.0 0.5; 0.5 1.0])
Y2 = from_moment([3.0, -1.0], [1.0 0.2; 0.2 1.5])

println("X2: mean = $(X2.mean)")
println("X2: covariance = $(X2.covariance)")
println("Y2: mean = $(Y2.mean)")
println("Y2: covariance = $(Y2.covariance)")

# Compute sum
Z2 = sum(X2, Y2)
println("Z2 = X2 + Y2: mean = $(Z2.mean)")
println("Z2: covariance = $(Z2.covariance)")

# Expected result: mean = [1+3, 2+(-1)] = [4, 1]
expected_mean_2d = [4.0, 1.0]
expected_cov_2d = [3.0 0.7; 0.7 2.5]
println("Expected mean: $expected_mean_2d")
println("Expected covariance: $expected_cov_2d")

# Check if results match
mean_correct = all(abs.(Z2.mean - expected_mean_2d) .< 1e-10)
cov_correct = all(abs.(Z2.covariance - expected_cov_2d) .< 1e-10)
println("Test passed: $(mean_correct && cov_correct)")

# Test 3: Edge case - zero mean
println("\nTest 3: Zero mean case")
println("-" ^ 30)

X3 = from_moment([0.0, 0.0], [1.0 0.0; 0.0 1.0])
Y3 = from_moment([0.0, 0.0], [1.0 0.0; 0.0 1.0])

Z3 = sum(X3, Y3)
println("X3 + Y3: mean = $(Z3.mean), covariance = $(Z3.covariance)")
println("Expected: mean = [0, 0], covariance = [2 0; 0 2]")
println("Test passed: $(all(Z3.mean .== 0) && all(Z3.covariance .== 2*I(2)))")


println("\n" ^ 2)
println("All tests completed!") 