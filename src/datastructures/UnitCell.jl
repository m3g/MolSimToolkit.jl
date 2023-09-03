
export UnitCell

"""
    UnitCell{T}

Container for the unit cell matrix. The container is used such that using the unit cell 
from a `Chemfiles.Frame` is transparent to the user, without any further conversions.

# Example

```julia-repl # to be doctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> trajectory = Trajectory(Testing.namd_traj);

julia> frame = currentframe(trajectory);

julia> unitcell = UnitCell(frame)
UnitCell{Float64}
  47.411   0.000   0.000
   0.000  47.411   0.000
   0.000   0.000  87.798

julia> positions = Positions(frame);

julia> wrap(positions[1], positions[2], unitcell)
3-element Point3D{Float64} with indices SOneTo(3):
  5.912472724914552
 10.768872261047363
 28.27700805664062
```

"""
struct UnitCell{T}
    matrix::SMatrix{3,3,T,9}
end
UnitCell(u::Chemfiles.UnitCell) = UnitCell(SMatrix{3,3,Float64,9}(transpose(Chemfiles.matrix(u))))
UnitCell(f::Chemfiles.Frame) = UnitCell(Chemfiles.UnitCell(f))

function Base.show(io::IO, u::UnitCell) 
    print(io, "$(typeof(u))")
    for i in 1:3
        print(io, "\n")
        for j in 1:3
            print(io, @sprintf("%8.3f", u.matrix[i,j]))
        end
    end
end

@testitem "UnitCell" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    using StaticArrays
    import Chemfiles
    traj = Chemfiles.Trajectory(Testing.namd_traj)
    frame = Chemfiles.read(traj)
    u = UnitCell(frame)
    @test u.matrix ≈ transpose(Chemfiles.matrix(Chemfiles.UnitCell(frame)))
    pos = Positions(frame)
    x = pos[1]
    y = pos[2]
    @test wrap(x, y, u) ≈ wrap(x, y, u.matrix)
    @test wrap_to_first(x, u) ≈ wrap_to_first(x, u.matrix)
end


