using Test
using Distributions
using LinearAlgebra

include("../src/gaussian.jl")

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

            @testset "add" begin

                # Case 1
                μ₁ = [1; 1; 1]
                S₁ = Matrix(I, 3, 3)
                𝑥₁ = from_sqrt_moment(μ₁, S₁)

                μ₂ = [2; 2; 2]
                S₂ = zeros(3, 3)
                𝑥₂ = from_sqrt_moment(μ₂, S₂)

                res = add(𝑥₁, 𝑥₂; sqrt=true)
                
                @test isapprox(res.mean, [3, 3, 3]; atol=1e-8)
                @test isapprox(res.covariance, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]; atol=1e-8)

                # Case 2
                μ₁ = [1; 1; 1]
                S₁ = zeros(3, 3)
                𝑥₁ = from_sqrt_moment(μ₁, S₁)

                μ₂ = [2; 2; 2]
                S₂ = Matrix(I, 3, 3)
                𝑥₂ = from_sqrt_moment(μ₂, S₂)

                res = add(𝑥₁, 𝑥₂; sqrt=true)

                @test isapprox(res.mean, [3, 3, 3]; atol=1e-8)
                @test isapprox(res.covariance, [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]; atol=1e-8)


                # Case 3
                μ₁ = [1; 1; 1]
                S₁ = Matrix(I, 3, 3) * 3.0
                𝑥₁ = from_sqrt_moment(μ₁, S₁)

                μ₂ = [2; 2; 2]
                S₂ = Matrix(I, 3, 3) * 4.0
                𝑥₂ = from_sqrt_moment(μ₂, S₂)

                res = add(𝑥₁, 𝑥₂; sqrt=true)

                @test isapprox(res.mean, [3, 3, 3]; atol=1e-8)
                @test isapprox(res.covariance, [-5.0 0.0 0.0; 0.0 -5.0 0.0; 0.0 0.0 -5.0]; atol=1e-8)

                # Case 4
                μ₁ = [1; 1; 1; 1]
                S₁ = [1 2 3 4; 0 5 6 7; 0 0 8 9; 0 0 0 16]
                𝑥₁ = from_moment(μ₁, S₁)

                μ₂ = [2; 2; 2; 2]
                S₂ = [0 10 -1 -3; 0 0 0 0; 0 0 0 0; 0 0 0 0]
                𝑥₂ = from_moment(μ₂, S₂)

                res = add(𝑥₁, 𝑥₂; sqrt=true)

                @test isapprox(res.mean, [3, 3, 3, 3]; atol=1e-5)
                @test isapprox(res.covariance, [1.0 2.0 3.0 4.0; 0.0 -11.180339 -1.788854 -0.447213; 0.0 0.0 -9.889388 -11.749968; 0.0 0.0 0.0 -16.023053]; atol=1e-5)

            end 

            @testset "log_pdf" begin
                        
                μ = [0.0, 0.0]
                Σ = [1.0 0.5; 0.5 1.5]
                x = [0.1, -0.2]
                
                S = cholesky(Σ).U

                pdf = Gaussian(μ, Matrix(S))
                
                logp_expected = logpdf(MvNormal(μ, Σ), x) # Compare to Distributions.jl
                logp = log_pdf(x, pdf; grad=false, sqrt=true)

                @test isapprox(logp, logp_expected; atol=1e-8)

                μ = [0.0, 0.0]
                Σ = [1e-10 0.0; 0.0 1e-10]
                x = [0.0, 0.0]

                S = cholesky(Σ).U
                pdf = Gaussian(μ, Matrix(S))

                logp_expected = logpdf(MvNormal(μ, Σ), x)  # Compare to Distributions.jl
                logp = log_pdf(x, pdf; grad=false, sqrt=true)
                @test isapprox(logp, logp_expected; atol=1e-8)

                μ = [0.0, 0.0]
                Σ = [2.0 0.3; 0.3 1.0]
                x = [10.0, -10.0]

                S = cholesky(Σ).U
                pdf = Gaussian(μ, Matrix(S))

                logp_expected = logpdf(MvNormal(μ, Σ), x)
                logp = log_pdf(x, pdf; grad=false, sqrt=true)
                @test isapprox(logp, logp_expected; atol=1e-8)

           end

            @testset "placeholder" begin
               @test skip=true
               @test skip=true
            end
        end;


        
