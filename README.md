## Resources

### Stochastic Differential Equations

MIT Lecture Series: https://www.youtube.com/watch?v=qdbkvD4N-us&t=21s&ab_channel=MITOpenCourseWare
Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations
Numerical Solution of Stochastic Differential Equations

Seems by far the best:

Applied Stochastic Differential Equations
Book by Arno Solin and Simo Särkkä

 Discrete-time integration of continuous-time stochastic differential equations (SDEs) with sampled noise increments, i.e., stochastic Runge-Kutta integration.

 https://www.cs.unc.edu/~welch/kalman/media/pdf/maybeck_ch1.pdf

 https://z-library.sk/book/2329291/72ed72/introduction-to-stochastic-differential-equations.html


### Gaussian Processes 

ML Tutorial Gaussian Process: https://www.youtube.com/watch?v=92-98SYOdlY&ab_channel=MarcDeisenroth
Visual Exploration of Gaussian Processes: https://distill.pub/2019/visual-exploration-gaussian-processes/

# GLMakie Rendering Issue Fix

Unsure which one of these fixed the rednering issue of fading lines but it was one of them. 

See here for example code: 

https://docs.makie.org/dev/explanations/cameras
 
# Set the NVIDIA adapter as default
export MESA_D3D12_DEFAULT_ADAPTER_NAME=NVIDIA

# Force the D3D12 driver
export GALLIUM_DRIVER=d3d12

# Disable software fallback
export LIBGL_ALWAYS_SOFTWARE=0

# ALL gausssian filtering on slide 48 at time 1:08:56 in video below 
https://www.youtube.com/watch?v=JUDqgLMqrgk&list=PLYzdsV9oalhtLVRx7Se1w-TvIxKInK97g&index=9&ab_channel=LukeThompson


old 
Unscented trasnform code: 

    # else   
    #     Σ𝑥 = density.covariance
    #     S𝑥 = cholesky((L + λ) * Σ𝑥).L  
    #     # Weights Wᵢ 𝒘 𝑾
        
    #     # Plus one because non-zero based indexing
    #     𝛘 = zeros(L, 2L + 1)
    #     𝛘[:, 1] = μ𝑥

    #     for i in 1:L
    #         𝛘[:, i+1] = μ𝑥 + S𝑥[:, i]
    #         𝛘[:, i+1+L] = μ𝑥 - S𝑥[:, i]
    #     end

    #     # Construct the weights matrix
    #     𝑾ᵐ = zeros(2L + 1)
    #     𝑾ᵐ[1] = λ / (L + λ)
    #     𝑾ᵐ[2:end] .= 1 / (2 * (L + λ))
        
    #     𝑾ᶜ = zeros(2L + 1)
    #     𝑾ᶜ[1] = λ / (L + λ) + (1 - α^2 + β)
    #     𝑾ᶜ[2:end] .= 1 / (2 * (L + λ))

    #     @assert abs(sum(𝑾ᵐ) - 1.0) < 1e-12 "Mean weights don't sum to 1"
        
    #     # for i in 1:(2L + 1)
    #     #     @show 𝛘[:, i]
    #     #     𝒴[i] = func(𝛘[:, i])
    #     # end 

    #     # for i in 1:(2 * L + 1)
    #     #     push!(𝒴, predict_measurement(𝛘[i])) # TODO: Compare the two methods of building Y
    #     # end

    #     # 𝒴 = hcat[predict_measurement(𝛘[:, i]) for i in 1:(2L + 1)]
    #     # μ𝑦 = sum(𝑾ᵐ[i] * 𝒴[:, i] for i in 1:(2L + 1))

    #     # Σ𝑦 = sum(𝑾ᶜ[i] * (𝒴[:, i] - μ𝑦) * (𝒴[:, i] - μ𝑦)' for i in 1:(2L + 1)) 

    #     # # Then compute mean as:
    #     # μ𝑦 = sum(𝑾ᵐ[i] * 𝒴[i] for i in 1:(2L + 1))

    #     # # And covariance as:
    #     # Σ𝑦 = sum(𝑾ᶜ[i] * (𝒴[i] - μ𝑦)^2 for i in 1:(2L + 1))

    #     # CLAUDE CODE
    #     # Propagate sigma points
    #     𝒴 = zeros(Float64, L, 2L + 1)

    #     for i in 1:(2L + 1)
    #         𝒴[:, i] = func(𝛘[:, i])
    #     end
        
    #     # Compute predicted mean and covariance (no Q added here)
    #     μ𝑦 = sum(𝑾ᵐ[i] * 𝒴[:, i] for i in 1:(2L + 1))
    #     Σ𝑦 = sum(𝑾ᶜ[i] * (𝒴[:, i] - μ𝑦) * (𝒴[:, i] - μ𝑦)' for i in 1:(2L + 1))
        
    #     # Ensure symmetry
    #     Σ𝑦 = 0.5 * (Σ𝑦 + Σ𝑦') # TODO: Review thoery 


 AND MORE 


        # # μ𝑦 = transformed_density.mean
        # # Σ𝑦 = transformed_density.covariance
        # μ𝑥 = density.mean
        # Σ𝑥 = density.covariance
        # L = length(μ𝑥)
        
        # # UKF parameters
        # κ = 0
        # α = 0.2
        # β = 2
        # λ = α^2 * (L + κ) - L
        
        # # Add regularization for numerical stability
        # ε = 1e-9
        # Σ𝑥_reg = Σ𝑥 + ε * I
        
        # # Generate sigma points from STATE distribution
        # Sₓ = cholesky((L + λ) * Σ𝑥_reg).L
        # 𝛘 = zeros(Float64, L, 2L + 1)
        # 𝛘[:, 1] = μ𝑥
        
        # for i in 1:L
        #     𝛘[:, i+1] = μ𝑥 + Sₓ[:, i]
        #     𝛘[:, i+1+L] = μ𝑥 - Sₓ[:, i]
        # end
        
        # # Weights
        # 𝑾ᵐ = zeros(2L + 1)
        # 𝑾ᶜ = zeros(2L + 1)
        # 𝑾ᵐ[1] = λ / (L + λ)
        # 𝑾ᶜ[1] = λ / (L + λ) + (1 - α^2 + β)
        # 𝑾ᵐ[2:end] .= 1 / (2 * (L + λ))
        # 𝑾ᶜ[2:end] .= 1 / (2 * (L + λ))

        # # Transform sigma points through measurement model
        # 𝒴 = zeros(2L + 1)  # Assuming scalar measurements

        # for i in 1:(2L + 1)
        #     𝒴[i] = predict_measurement(𝛘[:, i])
        # end
        
        # # Compute measurement statistics
        # μ𝑦 = sum(𝑾ᵐ[i] * 𝒴[i] for i in 1:(2L + 1))
        # Σ𝑦 = sum(𝑾ᶜ[i] * (𝒴[i] - μ𝑦)^2 for i in 1:(2L + 1)) + R
        
        # # Compute cross-covariance (state-measurement)
        # Σ𝑥𝑦 = sum(𝑾ᶜ[i] * (𝛘[:, i] - μ𝑥) * (𝒴[i] - μ𝑦) for i in 1:(2L + 1))
        
        # # Kalman update
        # K = Σ𝑥𝑦 / Σ𝑦
        # innovation = measurement - μ𝑦
        # μ_updated = μ𝑥 + K * innovation
        # Σ_updated = Σ𝑥 - K * Σ𝑥𝑦'  # More numerically stable than K * Σ𝑦 * K'
        
        # # Ensure symmetry
        # Σ_updated = 0.5 * (Σ_updated + Σ_updated')