"""
    distances(pdbfile::String, trajectory::Chemfiles.Trajectory, sel1, sel2; first=1, last=nothing, stride=1)

Function that calculates the distance between the center of mass of the selections in a trajectory.

"""
function distances(
    pdb_file::String, 
    trajectory::Trajectory, 
    indexes1::AbstractVector{Int}, 
    indexes2::AbstractVector{Int} 
)
    distances = zeros(nframes(trajectory))
    for frame in trajectory
        unitcell = UnitCell(frame)
        coordinates = Positions(frame)
        cm1 = centerofmass(coordinates, indexes1)
        cm2 = centerofmass(coordinates, indexes2)
        # wrap coordinates according to PBC. Note that here we use an
        # internal function of `CellListMap`. 
        y_wrapped = wrap(y,x,unitcell)
        # compute distance
        d = norm(y_wrapped-x)
        # add data to distance array
        distances[iframe] = d
        # read next frame from the trajectory file
    end
    return distances
end

@testitem "distances" begin
    using PDBTools
    using MolSimToolkit.Testing

    traj = Trajectory(Testing.namd_traj)
    pdb = Testing.namd_pdb
    i1 = selindex(pdb, "name CA and residue 1")
    i2 = selindex(pdb, "name CA and residue 2")
    @test distances(pdb, traj, i1, i2) ≈ 
        [3.85915490559465, 3.8968489296128737, 3.8413137354273266, 3.924179031288244, 3.9073096020456592]

end

