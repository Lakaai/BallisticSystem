using Test
using Distributions
using LinearAlgebra

include("Gaussian.jl")

@testset "Gaussian.jl" begin
           @testset "log_sqrt_pdf" begin
                
                μ = [0.0, 0.0]
                Σ = [1.0 0.5; 0.5 1.5]
                x = [0.1, -0.2]
                
                S = cholesky(Σ).U

                pdf = Gaussian(μ, Matrix(S))
                
                logp_expected = logpdf(MvNormal(μ, Σ), x) # Compare to Distributions.jl
                logp = log_sqrt_pdf(x, pdf)

                @test isapprox(logp, logp_expected; atol=1e-8)

                μ = [0.0, 0.0]
                Σ = [1e-10 0.0; 0.0 1e-10]
                x = [0.0, 0.0]

                S = cholesky(Σ).U
                pdf = Gaussian(μ, Matrix(S))

                logp_expected = logpdf(MvNormal(μ, Σ), x)
                logp = log_sqrt_pdf(x, pdf)
                @test isapprox(logp, logp_expected; atol=1e-8)

                μ = [0.0, 0.0]
                Σ = [2.0 0.3; 0.3 1.0]
                x = [10.0, -10.0]

                S = cholesky(Σ).U
                pdf = Gaussian(μ, Matrix(S))

                logp_expected = logpdf(MvNormal(μ, Σ), x)
                logp = log_sqrt_pdf(x, pdf)
                @test isapprox(logp, logp_expected; atol=1e-8)

           end
           @testset "log_pdf" begin
               @test skip=true
               @test skip=true
           end
       end;