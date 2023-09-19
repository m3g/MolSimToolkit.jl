"""
    procrustes(x,y)

Aligns two structures [sets of points in 3D space]. Solves
the "Procrustes" problem. Structures are expected to be of the same size, and the 
correspondence is assumed from the vector indices. 

Returns x aligned, by performing the rigid body transformation [rotation
and translation that minimizes the RMSD between x and y].

x, y, and xnew (return) are matrices of dimensions (n,3) 
(n is the number of points, 3 is the dimension of the space).

"""
using LinearAlgebra

function procrustes(x::AbstractMatrix{T}, y::AbstractMatrix{T})
    N = size(x,1)
    xvec = reinterpret(reshape, SVector{N,T}, x)
    yvec = reinterpret(reshape, SVector{N,T}, y)
    return procrustes(xvec, yvec)
end

function procrustes(x::AbstractVector{<:AbstractVector}, y::AbstractVector{<:AbtractVector}; mass = nothing)
    # to be implemented
end

function procrustes(x::AbstractVector{<:AbstractVector}, y::AbstractVector{<:AbtractVector}; mass = nothing)
    length(x) == length(y) || throw(DimensionMismatch("x and y must have the same length"))

    cmx = center_of_mass(x, mass)
    cmy = center_of_mass(y, mass)
    x .= x .- Ref(cmx)
    y .= y .- Ref(cmy)

    # These arrays can be of the size of the bijection, not necessarily of the
    # sizes of the complete input arrays. Let us see what interface we will provide
    # for the cases where the user wants to move other atoms than the ones 
    # explicitly involved in the bijection.
    xm = zeros(3,length(x))
    xp = zeros(3,length(x))

    for i in eachindex(x,y,xm,ym)
        xm[1:3,i] .= y[i] .- x[i]
        xp[1:3,i] .= y[i] .+ x[i]
    end

    q = zeros(MMatrix{4,4,eltype(xm),16})
    for i in eachindex(xm, xp)
      q[1,1] = q[1,1] + sum(abs2, @view(xm[1:3,i]))
      q[1,2] = q[1,2] + xp[2,i]*xm[3,i] - xm[2,i]*xp[3,i]
      q[1,3] = q[1,3] + xm[1,i]*xp[3,i] - xp[1,i]*xm[3,i]
      q[1,4] = q[1,4] + xp[1,i]*xm[2,i] - xm[1,i]*xp[2,i]
      q[2,2] = q[2,2] + xp[2,i]^2 + xp[3,i]^2 + xm[1,i]^2
      q[2,3] = q[2,3] + xm[1,i]*xm[2,i] - xp[1,i]*xp[2,i]
      q[2,4] = q[2,4] + xm[1,i]*xm[3,i] - xp[1,i]*xp[3,i]
      q[3,3] = q[3,3] + xp[1,i]^2 + xp[3,i]^2 + xm[2,i]^2
      q[3,4] = q[3,4] + xm[2,i]*xm[3,i] - xp[2,i]*xp[3,i]
      q[4,4] = q[4,4] + xp[1,i]^2 + xp[2,i]^2 + xm[3,i]^2
    end
    q[2,1] = q[1,2]
    q[3,1] = q[1,3]
    q[3,2] = q[2,3]
    q[4,1] = q[1,4]
    q[4,2] = q[2,4]
    q[4,3] = q[3,4]          
    q = SMatrix(q)

  # Computing the eigenvectors 'v' of the q matrix

  v = LinearAlgebra.eigvecs(q)

  # Compute rotation matrix
  
  u = zeros(MMatrix{3,3,Float64,9})
  u[1,1] = v[1,1]^2 + v[2,1]^2 - v[3,1]^2 - v[4,1]^2
  u[1,2] = 2. * ( v[2,1]*v[3,1] + v[1,1]*v[4,1] )
  u[1,3] = 2. * ( v[2,1]*v[4,1] - v[1,1]*v[3,1] )
  u[2,1] = 2. * ( v[2,1]*v[3,1] - v[1,1]*v[4,1] )
  u[2,2] = v[1,1]^2 + v[3,1]^2 - v[2,1]^2 - v[4,1]^2
  u[2,3] = 2. * ( v[3,1]*v[4,1] + v[1,1]*v[2,1] )
  u[3,1] = 2. * ( v[2,1]*v[4,1] + v[1,1]*v[3,1] )
  u[3,2] = 2. * ( v[3,1]*v[4,1] - v[1,1]*v[2,1] )
  u[3,3] = v[1,1]^2 + v[4,1]^2 - v[2,1]^2 - v[3,1]^2      
  u = SMatrix(u)

  # Rotate and translate vector x aligning it to y
  xnew = Ref(u) .* x .+ Ref(cmy)

  # Restore input arrays
  x .= x .+ Ref(cmx)
  y .= y .+ Ref(cmy)

  return xnew
end



