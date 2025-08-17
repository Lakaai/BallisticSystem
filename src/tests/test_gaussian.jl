using Test
using Distributions
using LinearAlgebra

include("../Gaussian.jl")

@testset "Gassian.jl" begin
            @testset "join" begin

                μ𝑥 = [1, 2, 3]
                Σ𝑥 = [4.0 1.5 -0.8; 1.5 3.2 0.6; -0.8 0.6 2.1]
                denisty_𝑥 = Gaussian(μ𝑥, Σ𝑥)

                μ𝑦 = [4, 5, 6]
                Σ𝑦 = [2.5 0.9 1.2; 0.9 5.0 -1.8; 1.2 -1.8 3.6] 
               
                density_𝑦 = Gaussian(μ𝑦, Σ𝑦)
                
                joint_density = join(denisty_𝑥, density_𝑦)

                expected_mean = [1, 2, 3, 4, 5, 6] 
                expected_covariance = [4.0  1.5 -0.8  0.0  0.0  0.0;
                                       1.5  3.2  0.6  0.0  0.0  0.0;
                                       -0.8  0.6  2.1  0.0  0.0  0.0;
                                       0.0  0.0  0.0  2.5  0.9  1.2;
                                       0.0  0.0  0.0  0.9  5.0 -1.8;
                                       0.0  0.0  0.0  1.2 -1.8  3.6]

                @test isapprox(joint_density.mean, expected_mean; atol=1e-8)
                @test isapprox(joint_density.covariance, expected_covariance; atol=1e-8)
            end
            
            @testset "placeholder" begin
               @test skip=true
               @test skip=true
            end
        end;