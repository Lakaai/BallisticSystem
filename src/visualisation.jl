using GLMakie
using GeometryBasics

figure = Figure()
ax = Axis3(figure[1, 1], 
    title = "3D Plot",
    limits = (-5, 5, -5, 5, -5, 5)  # (xmin, xmax, ymin, ymax, zmin, zmax)
)

# Add a scatter point
scatter!(ax, Point3f(1, 1, 5), color = :red, markersize = 20)

# Add coordinate axes for reference
lines!(ax, [Point3f(0.0,0.0,0.0), Point3f(2.0,0.0,0.0)], color = :red, linewidth = 3)   # X-axis
lines!(ax, [Point3f(0,0,0), Point3f(0,2,0)], color = :green, linewidth = 3) # Y-axis  
lines!(ax, [Point3f(0,0,0), Point3f(0,0,2)], color = :blue, linewidth = 3)  # Z-axis

# Add a frustum at position (1, 1, 1)
frustum_center = Point3f(1, 1, 1)
frustum_size = 0.5

# Create frustum edges (a truncated pyramid)
# Base rectangle
base_points = [
    Point3f(-frustum_size, -frustum_size, 0),
    Point3f(frustum_size, -frustum_size, 0),
    Point3f(frustum_size, frustum_size, 0),
    Point3f(-frustum_size, frustum_size, 0)
]

# Top rectangle (smaller)
top_points = [
    Point3f(-frustum_size/5, -frustum_size/5, frustum_size),
    Point3f(frustum_size/5, -frustum_size/5, frustum_size),
    Point3f(frustum_size/5, frustum_size/5, frustum_size),
    Point3f(-frustum_size/5, frustum_size/5, frustum_size)
]

# Move frustum to desired position
base_points = [p + frustum_center for p in base_points]
top_points = [p + frustum_center for p in top_points]

# Draw base rectangle
for i in 1:4
    next_i = mod(i, 4) + 1
    lines!(ax, [base_points[i], base_points[next_i]], color = :blue, linewidth = 2, linestyle = :dot)
end

# Draw top rectangle
for i in 1:4
    next_i = mod(i, 4) + 1
    lines!(ax, [top_points[i], top_points[next_i]], color = :blue, linewidth = 2, linestyle = :dot)
end

# Draw connecting edges
for i in 1:4
    lines!(ax, [base_points[i], top_points[i]], color = :blue, linewidth = 2, linestyle = :dot)
end

display(figure)