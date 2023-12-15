module PackmolInputCreator

using ..MolSimToolkit: center_of_mass
using PDBTools
using TestItems

export convert_concentration, convert_density_table
export density_pure_solvent, density_pure_cossolvent
export write_packmol_input
export SolutionBoxUSC

abstract type SolutionBox end
density_pure_solvent(system::SolutionBox) = system.density_table[begin, 2]
density_pure_cossolvent(system::SolutionBox) = system.density_table[end, 2]

const PackmolInputCreatorDirectory = @__DIR__

@kwdef mutable struct SolutionBoxUSC <: SolutionBox
    concentration_units::String = "x"
    solute_pdbfile::String
    solvent_pdbfile::String
    cossolvent_pdbfile::String
    density_table::Matrix{Float64}
    solute_molar_mass::Float64 = mass(readPDB(solute_pdbfile))
    solvent_molar_mass::Float64 = mass(readPDB(solvent_pdbfile))
    cossolvent_molar_mass::Float64 = mass(readPDB(cossolvent_pdbfile))
end

"""
    SolutionBoxUSC(; 
        solute_pdbfile::String, 
        solvent_pdbfile::String,
        cossolvent_pdbfile::String,
        density_table::Matrix{Float64},
        concentration_units = "x",
    )

Setup a system composed of a solute (U) a solvent (S) and a cossolvent (C). 

The concentration units of the density table can be provided explicitly and
are assumed by default to be the molar fraction, `x`, of the cossolvent.

"""
function SolutionBoxUSC end

function unit_name(u::String)
    u == "mol/L" && return "molarity"
    u == "x" && return "molar fraction"
    u == "vv" && return "volume fraction"
    u == "mm" && return "mass fraction"
end

function Base.show(io::IO, ::MIME"text/plain", system::SolutionBoxUSC)
    print(io, chomp("""
    ==================================================================
    SolutionBoxUSC properties (Solute + Solvent + Cossolvent):
    ==================================================================
        Solute pdb file: $(basename(system.solute_pdbfile))
        Solvent pdb file: $(basename(system.solvent_pdbfile))
        Cossolvent pdb file: $(basename(system.cossolvent_pdbfile))
        Density of pure solvent: $(density_pure_solvent(system)) g/mL
        Density of pure cossolvent: $(density_pure_cossolvent(system)) g/mL
        Molar masses: 
            solute: $(system.solute_molar_mass) g/mol
            solvent: $(system.solvent_molar_mass) g/mol
            cossolvent: $(system.cossolvent_molar_mass) g/mol
        Concentration units: $(system.concentration_units) ($(unit_name(system.concentration_units)))
        Cocentration range: $(first(system.density_table[:, 1])) - $(last(system.density_table[:, 1]))
    ==================================================================
    """))
end

# conversion factor from mL/mol to Å^3/molecule
const CMV = 1e24 / 6.02e23

# conversion factor from mol/L to molecules/Å^3
const CMC = 6.02e23 / 1e27

#=
    interpolate_concentration(system, x)

Obtain by interpolation the value of of ρ[:,2] that corresponds
to a the best estimate given a value of x corresponding to 
to the domain ρ[:,1].

=#
function interpolate_concentration(system, x)
    ρ = system.density_table
    i = findfirst(d -> d >= x, @view(ρ[:, 1]))
    i == firstindex(@view(ρ[:, 1])) && return ρ[i, 2]
    ρ[i] == x && return ρ[i, 2]
    dρdx = (ρ[i, 2] - ρ[i-1, 2]) / (ρ[i, 1] - ρ[i-1, 1])
    d = ρ[i-1, 2] + dρdx * (x - ρ[i-1, 1])
    return d
end

# Adjust x to be in the [0,1] range, to avoid problems with numerical precision
fixrange(x) = x < 0 ? 0 : (x > 1 ? 1 : x)

"""
    convert_concentration(
        system::SolutionBoxUSC,
        input_concentration, 
        units;
        density, # of the solution, required when converting from/to mol/L
    )

Convert concentration from one unit to another. The input
concentration is given in `input_concentration`, and the unit conversion 
is given by `units` keyword, that can be one of the following pairs:

The supported concentration units are:

- `"mol/L"`: molarity
- `"x"`: molar fraction
- `"vv"`: volume fraction
- `"mm"`: mass fraction

Conversion among types consists in passing the `units` keyword argument,
which is a pair of the form `"from" => "to"`, where `"from"` and `"to"`
are one of the supported units.

# Example

For example, to convert from molarity to molar fraction, use:

```julia
convert_concentration(system, 55.5, "mol/L" => "x")
```
where `system` is a `SolutionBoxUSC` object, and `55.5` is the molarity.

"""
function convert_concentration(
    system::SolutionBoxUSC,
    input_concentration::Real, 
    units::Pair{String,String};
    density::Union{Real,Nothing} = nothing, # of the solution
)

    (; solvent_molar_mass, cossolvent_molar_mass) = system

    # If the units didn't change, just return the input concentrations
    first(units) == last(units) && return input_concentration

    ρ = density # density of the solution
    ρc = density_pure_cossolvent(system) # density of the pure cossolvent
    ρw = density_pure_solvent(system) # density of the pure solvent
    Mw = solvent_molar_mass # molar mass of solvent 
    Mc = cossolvent_molar_mass # molar mass of the cossolvent

    # Check if density of solution was provided, in the case of molar
    # fraction conversions. Also, to converto to other units from molarity, we need ρ
    if isnothing(ρ)
        if last(units) == "mol/L" || last(units) == "mol/L" 
            throw(ArgumentError(
                """
                Density of solution is required to convert from/to molarity.
                Use the optional keyword argument `density` to provide it.
                """))
        end
    end

    # nc and nw are the molar concentrations
    if first(units) == "vv"
        if !(0 <= input_concentration <= 1)
            throw(ArgumentError("Volume fraction must be in the [0,1] range."))
        end
        # mL * g/ml / g/mol = mol 
        vv = input_concentration
        Nc = (ρc * vv) / Mc # number of cossolvent molecules in 1mL of ideal solution
        Nw = ρw * (1 - vv) / Mw # number of solvent molecules in 1mL of ideal solution
        if last(units) == "x"
            x = Nc / (Nc + Nw) # molar fraction
            return fixrange(x)
        end
        if last(units) == "mol/L"
            mt = Nc * Mc + Nw * Mw # mass of the solution (for 1 mol total)
            v = mt / (1000*ρ) # volume of the solution (for 1 mol total)
            return Nc / v # mol/L
        end
        if last(units) == "mm"
            return fixrange(Nc * Mc / (Nc * Mc + Nw * Mw))
        end
    end

    if first(units) == "x"
        if !(0 <= input_concentration <= 1)
            throw(ArgumentError("Molar fraction must be in the [0,1] range."))
        end
        x = input_concentration # molar fraction
        if last(units) == "mol/L"
            m = x * Mc + (1 - x) * Mw # g/mol: mass of the solution
            v = m / (ρ*1000) # L/mol: volume
            return x / v # mol/L: molarity of the solution 
        end
        if last(units) == "vv"
            vc = Mc *  x / ρc # Volume of x mols of pure cossolvent 
            vw = Mw * (1 - x) / ρw # Volume of (1-x) mols of pure solvent 
            vv = vc / (vc + vw) # volume fraction of cossolvent in ideal solution
            return fixrange(vv)
        end
        if last(units) == "mm"
            return fixrange(x * Mc / (x * Mc + (1 - x) * Mw))
        end
    end

    if first(units) == "mol/L"
        pure_c = 1000 * ρc / Mc  
        if !(0 <= input_concentration <= pure_c)
            throw(ArgumentError("Cossolvent molarity must be in the [0,$pure_c] range."))
        end
        nc = input_concentration / 1000
        nw = (ρ - nc * Mc) / Mw
        if last(units) == "x"
            return fixrange(nc / (nc + nw))
        end
        if last(units) == "vv"
            vc = nc * Mc / ρc
            vw = nw * Mw / ρw
            return fixrange(vc / (vc + vw))
        end
        if last(units) == "mm"
            return fixrange(nc * Mc / (nc * Mc + nw * Mw))
        end
    end

    if first(units) == "mm"
        if !(0 <= input_concentration <= 1)
            throw(ArgumentError("Mass fraction must be in the [0,1] range."))
        end
        mm = input_concentration # mass fraction
        Nc = mm * 1 / Mc # mol of ethanol in 1g
        Nw = (1 - mm) / Mw # mol of water in 1g
        if last(units) == "x"
            return fixrange(Nc / (Nc + Nw)) # molar fraction
        end
        if last(units) == "vv"
            vc = Nc * Mc  / ρc 
            vw = Nw * Mw / ρw
            return fixrange(vc / (vc + vw)) # volume fraction
        end
        if last(units) == "mol/L"
            v = 1 / ρ # volume of 1g of solution
            return 1000 * Nc / v # mol/L
        end
    end
end

"""
    convert_density_table(system::SolutionBoxUSC, target_units)

Converts the density table of the system from one unit to another. Returns the 
input `system` with the density table converted to the new units.

The target units may be one of: `"mol/L"`, `"x"`, `"vv"`, `"mm"`.

## Example

```julia
convert_density_table(system, "mol/L")
```

"""
function convert_density_table(
    system::SolutionBoxUSC,
    target_units::String;
)
    current_units = system.concentration_units
    density_table = system.density_table
    for irow in eachindex(eachrow(density_table))
        cin = density_table[irow, 1]
        ρ = density_table[irow, 2]
        cout = convert_concentration(system, cin, current_units => target_units; density = ρ)
        density_table[irow, 1] = cout
    end
    system.concentration_units = target_units
    return system
end

"""
    write_packmol_input(
        system::SolutionBoxUSC;
        concentration::Real, 
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Real}, # or
        margin::Real
    )

Function that generates an input file for Packmol. 

The box sides are given in Ångströms, and can be provided as a vector of 3 elements.
Alternativelly, the margin can be provided, and the box sides will be calculated as
the maximum and minimum coordinates of the solute plus the margin in all 3 dimensions.

"""
function write_packmol_input(
    system::SolutionBoxUSC;
    concentration::Real, 
    input="box.inp",
    output="system.pdb",
    box_sides::Union{AbstractVector{<:Real},Nothing} = nothing,
    margin::Union{Real,Nothing} = nothing, 
    # testing option
    debug = false,
)

    (; solute_pdbfile,
       solvent_pdbfile,
       cossolvent_pdbfile,
       solute_molar_mass,
       solvent_molar_mass,
       cossolvent_molar_mass,
    ) = system

    # molar masses (g/mol)
    Mp = solute_molar_mass
    Mc = cossolvent_molar_mass
    Mw = solvent_molar_mass

    # Check consistency of the concentrations given
    cunit = system.concentration_units
    ρs = @view(system.density_table[:, 1])
    if cunit == "mol/L" && (first(ρs) < 1e-3 && last(ρs) ≈ 1)
        throw(ArgumentError("Concentrations in density table do not appear to be in mol/L."))
    end
    if (cunit in ("x", "mm", "vv")) && (any(x -> !(0 <= x <= 1), ρs)) 
        throw(ArgumentError("Concentrations in density table outside [0,1] range, and units are: $cunit"))
    end

    # Find the density corresponding to the target concentration
    ρ = interpolate_concentration(system, concentration)

    # Obtain the concentration in all units, for testing
    c_x = convert_concentration(system, concentration, cunit => "x"; density = ρ)
    c_vv = convert_concentration(system, concentration, cunit => "vv"; density = ρ)
    cc_mol = convert_concentration(system, concentration, cunit => "mol/L"; density = ρ)
    cc_mm = convert_concentration(system, concentration, cunit => "mm"; density = ρ)

    # aliases for clearer formulas
    ρc = density_pure_cossolvent(system) 
    ρw = density_pure_solvent(system) 

    # Convert cossolvent concentration in molecules/Å³
    cc = CMC * cc_mol

    # Set box side
    if isnothing(box_sides) && isnothing(margin)
        throw(ArgumentError("Either box_sides or margin must be provided."))
    elseif !isnothing(box_sides) && !isnothing(margin)
        throw(ArgumentError("Either box_sides or margin must be provided, but not both."))
    end
    solute_atoms = readPDB(system.solute_pdbfile)
    solute_extrema = round.(maxmin(solute_atoms).xlength; digits=3)
    if !isnothing(margin)
        box_sides = solute_extrema .+ 2 .* margin
    end

    # Box volume (Å³)
    vbox = box_sides[1] * box_sides[2] * box_sides[3]

    # Solution volume (vbox - vsolute) - vsolute is estimated
    # as if it had the same mass density of the solution
    vs = vbox - CMV * Mp / ρ

    # number of cossolvent molecules: cossolvent concentration × volume of the solution
    nc = round(Int, cc * vs)

    # Number of solvent molecules
    nw = round(Int, nc * (1 - c_x) / c_x)

    # Final density of the solution (not inclusing solute volume)
    ρ = CMV * (Mc * nc + Mw * nw) / vs

    # Final cossolvent concentration (mol/L)
    cc_f = 1000 * (nc / vs) * CMV

    # Final solvent concentration (mol/L)
    cw_f = 1000 * (nw / vs) * CMV

    # Half of box sides, to center the solute at the origin
    l = round.(box_sides ./ 2; digits=3)

    summary = """
        ==================================================================
        Summary:
        ==================================================================

        Target concentration = $cc_mol mol/L
        (of cossolvent)      = $c_x molar fraction
                             = $c_vv volume fraction
                             = $cc_mm mass fraction 
                             = $cc molecules/Å³

        Box volume = $vbox Å³
        Solution volume = $vs Å³   
        Solute extrema = [ $(join(-0.5*solute_extrema, ", ")), $(join(0.5*solute_extrema, ", ")) ] Å
        Box sides = [ $(join( -1.0*l, ", ")), $(join( l, ", ")) ] Å 

        Density of solution = $ρ g/mL
        Solute molar mass = $Mp g/mol
        Cossolvent molar mass = $Mc g/mol
        Solvent molar mass = $Mw g/mol

        Number of cossolvent ($(basename(cossolvent_pdbfile))) molecules = $nc 
        Number of solvent ($(basename(solvent_pdbfile))) molecules = $nw 

        Final cossolvent concentration = $cc_f mol/L
                                    = $(CMC*cc_f) molecules/Å³
        Final solvent concentration = $cw_f mol/L
                                    = $(CMC*cw_f) molecules/Å³
                                    
        Final volume fraction = $((nc * Mc / ρc)/((nc * Mc / ρc) + (nw * Mw / ρw)))
        Final molar fraction = $(nc/(nc+nw))

        ==================================================================
        """
    println(summary)

    open(input, "w") do io
        print(io,
            """
            # 
            # Packmol input file
            # 
            # Generated by MolSimToolkit.jl
            #
            """
        )
        for line in split(summary, "\n")
            println(io, "# $line")
        end
        println(io,
            """
            #
            tolerance 2.0
            output $output
            add_box_sides 1.0
            filetype pdb
            seed -1

            structure $solute_pdbfile
              number 1
              center
              fixed 0. 0. 0. 0. 0. 0.
            end structure

            structure $solvent_pdbfile
              number $nw
              inside box $(join( -1.0*l, " ")), $(join( l, " "))  
            end structure
            """)
        if nc > 0
            println(io,
                """
                structure $cossolvent_pdbfile
                  number $nc
                  inside box $(join( -1.0*l, " ")), $(join( l, " "))  
                end structure
                """)
        end
    end
    print(chomp(
        """
        Wrote file: $input

        ==================================================================
        """))
    
    if debug 
        return nw, nc, l
    else
        return nothing
    end
end # function write_packmol_input

# Tests
include("./test/runtests.jl")

end # module PackmolInputCreator

