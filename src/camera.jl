# Define the struct using kwdef to allow for default values and keyword based construction. 
@kwdef mutable struct Camera
    camera_matrix::Matrix{T} = zeros(3, 3)
    distortion_coefficients::Vector{T} = zeros(5)
end 

# function Camera(camera_matrix::Matrix{T}) where {T}
#     @assert size(camera_matrix, 1) == 3 "camera_matrix must be 3x3"
#     @assert size(camera_matrix, 2) == 3 "camera_matrix must be 3x3"
#     Camera{T}(camera_matrix)
# end 

# function Camera(f::T, cx::T, cy::T) where {T}