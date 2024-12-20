#
# Testing routines: computing the lists with naive algorithms
#

import CellListMap: Box, wrap_relative_to

# For a single set of molecules
function naive_md(x, xn_atoms_per_molecule::Integer, unitcell::AbstractVecOrMat, cutoff::Real)
    box = Box(unitcell, cutoff)
    naive_md(x, xn_atoms_per_molecule, box)
end
function naive_md(x, xn_atoms_per_molecule::Integer, box::Box)
    x_list = init_list(x, i -> _mol_indices(i, xn_atoms_per_molecule))
    for i in 1:length(x)-1
        vx = x[i]
        for j in i+1:length(x)
            vy = x[j]
            imol = _mol_indices(i, xn_atoms_per_molecule)
            jmol = _mol_indices(j, xn_atoms_per_molecule)
            if imol == jmol
                continue
            end
            vy_wrapped = wrap_relative_to(vy, vx, box.input_unit_cell.matrix)
            d = sqrt(sum(abs2, vy_wrapped - vx))
            if d > box.cutoff
                continue
            end
            if d < x_list[imol].d
                x_list[imol] = MinimumDistance(true, i, j, d)
            end
            if d < x_list[jmol].d
                x_list[jmol] = MinimumDistance(true, j, i, d)
            end
        end
    end
    return x_list
end

# For disjoint sets, returning only one list
function naive_md(x, y, xn_atoms_per_molecule::Integer, unitcell::AbstractVecOrMat, cutoff::Real)
    box = Box(unitcell, cutoff)
    naive_md(x, y, xn_atoms_per_molecule, box)
end
function naive_md(x, y, xn_atoms_per_molecule, box)
    x_list = init_list(x, i -> _mol_indices(i, xn_atoms_per_molecule))
    for (i, vx) in pairs(x)
        for (j, vy) in pairs(y)
            vy_wrapped = wrap_relative_to(vy, vx, box.input_unit_cell.matrix)
            d = sqrt(sum(abs2, vy_wrapped - vx))
            if d > box.cutoff
                continue
            end
            imol = _mol_indices(i, xn_atoms_per_molecule)
            if d < x_list[imol].d
                x_list[imol] = MinimumDistance(true, i, j, d)
            end
        end
    end
    return x_list
end

# For disjoint sets, returning both lists
function naive_md(x, y, xn_atoms_per_molecule::Integer, yn_atoms_per_molecule::Integer, unitcell::AbstractVecOrMat, cutoff::Real)
    box = Box(unitcell, cutoff)
    naive_md(x, y, xn_atoms_per_molecule, yn_atoms_per_molecule, box)
end
function naive_md(x, y, xn_atoms_per_molecule, yn_atoms_per_molecule, box)
    x_list = init_list(x, i -> _mol_indices(i, xn_atoms_per_molecule))
    y_list = init_list(y, i -> _mol_indices(i, yn_atoms_per_molecule))
    for (i, vx) in pairs(x)
        for (j, vy) in pairs(y)
            vy_wrapped = wrap_relative_to(vy, vx, box.input_unit_cell.matrix)
            d = sqrt(sum(abs2, vy_wrapped - vx))
            if d > box.cutoff
                continue
            end
            imol = _mol_indices(i, xn_atoms_per_molecule)
            if d < x_list[imol].d
                x_list[imol] = MinimumDistance(true, i, j, d)
            end
            jmol = _mol_indices(j, yn_atoms_per_molecule)
            if d < y_list[jmol].d
                y_list[jmol] = MinimumDistance(true, j, i, d)
            end
        end
    end
    return x_list, y_list
end

#
# Comparison functions
#
import Base.isapprox
function isapprox(
    list1::List,
    list2::List,
)
    length(list1) != length(list2) && return false
    for i in eachindex(list1)
        list1[i].within_cutoff != list2[i].within_cutoff && return false
        list1[i].i != list2[i].i && return false
        list1[i].j != list2[i].j && return false
        !(list1[i].d ≈ list2[i].d) && return false
    end
    return true
end
function isapprox(
    lists1::Tuple{T,T},
    lists2::Tuple{T,T},
) where {T<:List}
    !(lists1[1] ≈ lists2[1]) && return false
    !(lists1[2] ≈ lists2[2]) && return false
    return true
end

#
# Auxiliary functions to plot a 2D molecular arrangement for 
# inspection and visualization
#
function plot_mol!(
    p,
    x;
    cycle=false,
    markercolor=nothing
)
    for i in 1:length(x)-1
        xat = [x[i][1], x[i+1][1]]
        yat = [x[i][2], x[i+1][2]]
        Main.plot!(p, xat, yat, label=:none, linewidth=2, color=:black)
    end
    if cycle
        xat = [x[end][1], x[1][1]]
        yat = [x[end][2], x[1][2]]
        Main.plot!(p, xat, yat, label=:none, linewidth=2, color=:black)
    end
    if isnothing(markercolor)
        Main.scatter!(p, Tuple.(x), label=:none, markersize=4.0)
    else
        Main.scatter!(p, Tuple.(x), label=:none, markersize=4.0, markercolor=markercolor)
    end
    return p
end

function plot_mol!(
    p,
    x,
    n_atoms_per_molecule::Integer;
    cycle=false,
    markercolor=nothing
)
    for i in 1:n_atoms_per_molecule:length(x)
        plot_mol!(
            p,
            @view(x[i:i+n_atoms_per_molecule-1]);
            cycle=cycle,
            markercolor=markercolor
        )
        if isnothing(markercolor)
            markercolor = p.series_list[end][:markercolor]
        end
    end
    return p
end

function plot_md!(
    p,
    x,
    xn_atoms_per_molecule::Integer,
    y,
    yn_atoms_per_molecule::Integer,
    md::List;
    x_cycle=false,
    y_cycle=false
)
    for pair in md
        !pair.within_cutoff && continue
        xp = [x[pair.i][1], y[pair.j][1]]
        yp = [x[pair.i][2], y[pair.j][2]]
        Main.plot!(p, xp, yp, label=:none, linewidth=2, color=:black, alpha=0.3)
    end
    plot_mol!(p, x, xn_atoms_per_molecule; cycle=x_cycle, markercolor=:red)
    plot_mol!(p, y, yn_atoms_per_molecule; cycle=y_cycle, markercolor=:blue)
end

download_example() =
    Main.PDBTools.readPDB(download("https://raw.githubusercontent.com/m3g/ComplexMixtures.jl/master/test/data/NAMD/structure.pdb"))
