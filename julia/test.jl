# using Optim
# using LinearAlgebra

# # # Define prediction density
# # struct Gaussian
# #     μ
# #     Σ
# # end

# function logLikelihood(x, density, measurement)
    
#     # Eigen::MatrixXd dhdx;
#     dhdx = zeros(1, 3)
#     # auto likelihood = predictDensity(x, system, dhdx);
#     likelihood, dhdx = predictDensity(x, density, dhdx)

#     log_likelihood, log_likelihood_grad =  logpdf(measurement, likelihood, grad=true)


#     # // Evaluate log N(y; h(x), R) and d/dy log N(y; h(x), R)
#     # Eigen::VectorXd loglikGrad;
#     # double logLik = likelihood.log(y_, loglikGrad);
#     # // Note:
#     # //  d                        d     
#     # // -- log N(y; h(x), R) = - -- log N(y; h(x), R)
#     # // dh                       dy

#     # // Gradient of log likelihood:
#     # //
#     # //         d 
#     # // g_i = ---- log N(y; h(x), R)
#     # //       dx_i
#     # //
#     # //             dh_k     d
#     # // g_i = sum_k ---- * ---- log N(y; h(x), R)
#     # //             dx_i   dh_k
#     # //
#     # //               dh_k     d
#     # // g_i = - sum_k ---- * ---- log N(y; h(x), R)
#     # //               dx_i   dy_k
#     # //

#     g = -dhdx' * log_likelihood_grad;
#     return log_likelihood, g;

# end 

# function logpdf(x, pdf::Gaussian; grad=false)
#     μ = pdf.μ isa AbstractVector ? pdf.μ : [pdf.μ]
#     S = pdf.Σ
#     @show S
#     n = length(x)
#     @show x
#     @show μ
#     @show S

#     @assert istriu(S) "S is not upper triangular"
#     @assert length(x) == length(μ) "Input x and mean μ must have same length"

#     Δ = x - μ
#     w = LowerTriangular(transpose(S)) \ Δ   
#     println(diag(S))
#     logpdf = -(n/2)*log(2*π)-(1/2)*sum(log.(abs.(diag(Matrix(S)))))-(1/2)*dot(w,w)

#     if grad
#         gradient = -UpperTriangular(S) \ w      # Gradient ∇logp = -S⁻¹ * w
#         return logpdf, gradient
#     else 
#         return logpdf                           # Return log N(x; μ, S)
#     end 
# end 


# # Measurement model
# function h(x)
#     return x  # Identity model for simplicity
# end

# R = 0.1I  # Measurement noise covariance


# # Evaluate h(x) from the measurment model y = h(x) + v
# function predict_measurement(x, density)
#     r1 = 5000; # Horizontal position of sensor [m]
#     r2 = 5000; # Vertical position of sensor [m]
#     h1 = x[1]
#     h = sqrt(r1^2 + (h1-r2)^2)
#     return h
# end 

# # Evaluate h(x) and its Jacobian J = dh/fx from the measurement model y = h(x) + v
# function predict_measurement(x, density, dhdx)
#     r1 = 5000; # Horizontal position of sensor [m]
#     r2 = 5000; # Vertical position of sensor [m]
#     h = predict_measurement(x, density)
#     dhdx = zeros(1, 3)
#     @show x
#     println("inside predict_measurement")
#     h1 = x[1]
#     range = sqrt(r1^2 + (h1-r2)^2)

#     dhdx[1, 1] = (h1 - r2)/range
#     return h, dhdx
# end

# function predictDensity(x, density, dhdx)
#     sigma_rng = 50.0 # 50m standard deviation

#     h, dhdx = predict_measurement(x, density, dhdx)
#     # SR = Matrix(Diagonal(1.0)) * sigma_rng
#     SR = Diagonal([50.0])
#     # TODO: Fix (see lab5 cpp)
#     @show length(x)
#     @show h
#     @show SR
#     return Gaussian(h, SR), dhdx
# end 


# # function cost_function_factory(z, prediction_density::Gaussian)

# #     μ_pred = prediction_density.μ
# #     Σ_pred_inv = inv(prediction_density.Σ)
# #     R_inv = inv(0.1I)  # measurement noise

# #     return function(x)
# #         # Negative log-prior
# #         diff = x .- μ_pred
# #         prior = 0.5 * dot(diff, Σ_pred_inv * diff)

# #         # Negative log-likelihood
# #         h_x = x  # identity measurement function for simplicity
# #         residual = h_x .- z
# #         likelihood = 0.5 * dot(residual, R_inv * residual)

# #         return prior + likelihood
# #     end
# # end

# function cost_function_factory(density::Gaussian, measurement)
#     return function(x) # Returns a cost function f(x) which has the required signature for the optimiser
#         logprior = logpdf(x, density)

#         # You must define this based on your measurement model
#         loglik, _ = logLikelihood(x, density, measurement)
#         @show logprior
#         @show loglik
#          # Return a scalar cost (negative log-likelihood (the measurement cost) + log-prior (the prediction cost))
#         return -(logprior + loglik) # Return −logp(x∣z) = -(logp(x) + logp(z∣x))
#     end
# end

# measurement = [0.9]
# density = Gaussian([1.0, 2.0], [1.0 0.0; 0.0 1.0])

# cost = cost_function_factory(density, measurement)

# # Now this works:
# result = optimize(cost, [0.0, 0.0], BFGS())




# # If you're not explicitly passing them into the cost function, they must be captured by closure or global scope, 
# # which is often bad practice or unintentional.