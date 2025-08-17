function rotx(ϕ)
    return [1 0 0; 0 cos(ϕ) -sin(ϕ); 0 sin(ϕ) cos(ϕ)]
end

function roty(θ)
    return [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]
end

function rotz(ψ)
    return [cos(ψ) -sin(ψ) 0; sin(ψ) cos(ψ) 0; 0 0 1]
end

"""
    rotxyz(ϕ, θ, ψ)

Compute the rotation matrix for a rotation around the x, y, and z axes as a passive rotation.

# Arguments
- `ϕ`: Rotation around the x-axis.
- `θ`: Rotation around the y-axis.
- `ψ`: Rotation around the z-axis.
"""
function rotxyz(ϕ, θ, ψ)
    return rotz(ψ) * roty(θ) * rotx(ϕ)
end