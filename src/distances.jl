"""
    distances(pdbfile::String, trajectory::Chemfiles.Trajectory, sel1, sel2; first=1, last=nothing, stride=1)

Function that calculates the distance between the center of mass of the selections in a trajectory.

"""
function distances(
    pdb_file::String, 
    trajectory_file::String, 
    selection1, 
    selection2; 
    first=1, last=nothing, stride=1
)
    trajectory = Chemfiles.Trajectory(trajectory_file)
    atoms1 = PDBTools.readPDB(pdb_file, selection1)
    atoms2 = PDBTools.readPDB(pdb_file, selection2)
    last = isnothing(last) ? Chemfiles.length(trajectory) : last
    nframes = floor(Int,(last - first + 1) / stride)
    @show nframes, last, first, stride
    distances = zeros(nframes)
    # read first frame
    frame = Chemfiles.read(trajectory)
    for iframe in 1:min(Chemfiles.length(trajectory), last)
        if (iframe < first) || ((iframe % stride) != 0)
            Chemfiles.read!(trajectory, frame)
            continue
        end
        unitcell_read = Chemfiles.matrix(Chemfiles.UnitCell(frame))
        unit_cell = SMatrix{3,3}(transpose(unitcell_read))
        # read coordinates, convert them to small static vectors
        coordinates = Chemfiles.positions(frame)
        # TODO: below we compute the distance between the first atom of each selection
        # this has to be replaced by the computation of their centers of mass
        x = SVector{3}(@view(coordinates[:,atoms1[1].index]))
        y = SVector{3}(@view(coordinates[:,atoms2[1].index]))
        # wrap coordinates according to PBC. Note that here we use an
        # internal function of `CellListMap`. 
        y_wrapped = wrap(y,x,unit_cell)
        # compute distance
        d = norm(y_wrapped-x)
        # add data to distance array
        distances[iframe] = d
        # read next frame from the trajectory file
        if iframe < last  
            Chemfiles.read!(trajectory, frame)
        end
    end
    Chemfiles.close(trajectory)
    return distances
end

@testitem "distances" begin
    using Chemfiles
    using PDBTools
    using MolSimToolkit.Testing

    traj = Testing.namd_traj
    pdb = Testing.namd_pdb
    @test distances(pdb, traj, "name CA and residue 1", "name CA and residue 2") ≈ 
        [3.85915490559465, 3.8968489296128737, 3.8413137354273266, 3.924179031288244, 3.9073096020456592]

end

