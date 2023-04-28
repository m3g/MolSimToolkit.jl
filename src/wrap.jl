"""
    wrap(x::SVector{3}, x_ref::SVector{3}, unitcell::Chemfiles.UnitCell)

Function that wraps a vector `x` according to the periodic boundary conditions of the unit cell `unitcell` and the reference position `x_ref`.

```

"""
function wrap(x::SVector{3}, x_ref::SVector{3}, unitcell::Chemfiles.UnitCell)
    unit_cell_matrix = SMatrix{3,3}(transpose(Chemfiles.matrix(unitcell)))
    x = wrap_relative_to(x,x_ref,unit_cell_matrix)
    return x
end

@testitem "wrap" begin
    # Test wrapping of a vector
    x = SVector{3}([0.5, 0.5, 0.5])
    x_ref = SVector{3}([0.0, 0.0, 0.0])
    unitcell = Chemfiles.UnitCell([ 10.0 0.0 0.0 
                                    0.0 10.0 0.0
                                    0.0 0.0 10.0 ])    
    @test wrap(x, x_ref, unitcell) ≈ x
    # Test if vector is not inside unit cell
    x = SVector{3}([10.5, 10.5, 10.5])
    x_ref = SVector{3}([0.0, 0.0, 0.0])
    @test wrap(x, x_ref, unitcell) ≈ SVector{3}([0.5, 0.5, 0.5])
    # Test with a triclinic cell
    unitcell = Chemfiles.UnitCell([ 10.0 5.0 5.0 
                                    0.0 10.0 5.0
                                    0.0 0.0 10.0])
    x = SVector{3}([0.5, 0.5, 0.5])
    x_ref = SVector{3}([0.0, 0.0, 0.0])
    @test wrap(x, x_ref, unitcell) ≈ SVector{3}([0.5, 0.5, 0.5])
    # Test with atoms far from each other
    x = SVector{3}([10.5, 10.5, 10.5])
    x_ref = SVector{3}([0.0, 0.0, 0.0])
    @test wrap(x, x_ref, unitcell) ≈ SVector{3}([0.5, 0.5, 0.5])
end
