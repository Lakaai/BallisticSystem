using Test
using Distributions
using LinearAlgebra

include("../system_estimator.jl")

@testset "system_estimator.jl" begin
            @testset "sqrt_affine_transform" begin
                
                # TODO: Implement test

                μ = [0.0, 0.0, 0.0]
                S = [1.0 0.5 1.0; 0.0 1.0 0.5; 0.0 0.0 1.0]
                x = [0.1, -0.2, 0.3]          
                Δₜ = 0.1
                Q = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
                pdf = Gaussian(μ, Matrix(S))
                py = sqrt_affine_transform(rk4_step, pdf, Δₜ, Q) # Proability distribution y = Ν(μy, Sy)
                
                expected_mean = [0.0, 0.0, 0.0] # TODO
                expected_sqrt_covariance = [1.0 0.5 1.0; 0.0 1.0 0.5; 0.0 0.0 1.0] # TODO

                @test isapprox(py.mean, expected_mean; atol=1e-8)
                @test isapprox(py.covariance, expected_sqrt_covariance; atol=1e-8)
            end
            @testset "placeholder" begin
               @test skip=true
               @test skip=true
            end
        end;