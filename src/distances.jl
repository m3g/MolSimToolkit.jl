"""
    distances(trajectory, sel1, sel2; first=1, last=nothing, stride=1)

Function that calculates the distance between the center of mass of the selections in a trajectory.

"""
function distances(trajectory::Chemfiles.Trajectory, sel1, sel2; first=1, last=nothing, stride=1)
    atoms1 = select(pdb, selection1)
    atoms2 = select(pdb, selection2)
    last = isnothing(last) ? Chemfiles.length(trajectory) : last
    nframes = floor(Int,(last - first + 1) / stride)
    distances = zeros(nframes)
    for iframe in 1:nframes
        # skip strided frames
        if (iframe % stride) != 0
            continue
        end
        # If last frame was reached, stop
        if iframe > last
            break
        end
        # read unit cell matrix from trajectory
        frame = Chemfiles.frame(trajectory, iframe)
        unitcell = Chemfiles.UnitCell(frame)
        unit_cell = SMatrix{3,3}(transpose(matrix_read))
        # read coordinates, convert them to small static vectors
        coordinates = Chemfiles.positions(frame)
        x = SVector{3}(@view(coordinates[:,atom1[1].index]))
        y = SVector{3}(@view(coordinates[:,atom2[1].index]))
        # wrap coordinates according to PBC. Note that here we use an
        # internal function of `CellListMap`. 
        y_wrapped = wrap(y,x,unit_cell)
        # compute distance
        d = norm(y_wrapped-x)
        # add data to distance array
        push!(distances, d)
    end
    return distances
end

