# include("camera.jl")

camera_matrix = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
distortion_coefficients = [0.0, 0.0, 0.0, 0.0, 0.0]

camera = Camera(camera_matrix, distortion_coefficients)

