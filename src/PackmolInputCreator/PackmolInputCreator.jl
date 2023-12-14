module PackmolInputCreator

using ..MolSimToolkit: center_of_mass
using PDBTools
using TestItems

export convert_concentration
export density_pure_solvent, density_pure_cossolvent
export write_packmol_input
export SolutionBoxUSC

abstract type SolutionBox end
density_pure_solvent(system::SolutionBox) = system.density_table[begin, 2]
density_pure_cossolvent(system::SolutionBox) = system.density_table[end, 2]

const PackmolInputCreatorDirectory = @__DIR__

@kwdef struct SolutionBoxUSC <: SolutionBox
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
    )

Setup a system composed of a solute (U) a solvent (S) and a cossolvent (C). 

"""
function SolutionBoxUSC end

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

"""
    convert_concentration(
        system::SolutionBoxUSC,
        input_concentration, 
        units;
        density, # of the solution, required if the input concentration is mol/L
    )

Convert concentration from one unit to another. The input
concentration is given in `input_concentration`, and the unit conversion 
is given by `units` keyword, that can be one of the following pairs:

- `"mol/L"` => `"x"`: from molarity to molar fraction
- `"mol/L"` => `"vv"`: from molarity to volume fraction
- `"x"` => `"mol/L"`: from molar fraction to molarity
- `"x"` => `"vv"`: from molar fraction to volume fraction
- `"vv"` => `"mol/L"`: from volume fraction to molarity
- `"vv"` => `"x"`: from volume fraction to molar fraction

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
            return x
        end
        if last(units) == "mol/L"
            mt = Nc * Mc + Nw * Mw # mass of the solution (for 1 mol total)
            v = mt / (1000*ρ) # volume of the solution (for 1 mol total)
            return Nc / v # mol/L
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
            return vv
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
            return nc / (nc + nw)
        end
        if last(units) == "vv"
            vc = nc * Mc / ρc
            vw = nw * Mw / ρw
            return vc / (vc + vw)
        end
    end

end

"""
    write_packmol_input(
        system::SolutionBoxUSC;
        concentration::Real, 
        cunit="mol/L",
        input="box.inp",
        output="system.pdb",
        # box size
        box_sides::AbstractVector{<:Real}, # or
        margin::Real
    )

Function that generates an input file for Packmol. By default, 
the concentrations is given in `mol/L`, but it can also be 
given in molar fraction "x" or volume fraction "vv", using `cunit="x"` or `cunit="vv"`. 

The box sides are given in Ångströms, and can be provided as a vector of 3 elements.
Alternativelly, the margin can be provided, and the box sides will be calculated as
the maximum and minimum coordinates of the solute plus the margin in all 3 dimensions.

"""
function write_packmol_input(
    system::SolutionBoxUSC;
    concentration::Real, 
    cunit="mol/L",
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

    # If the concentration was given in mol/L we need to find to which
    # this corresponds in molar fraction
    c_x = if cunit == "mol/L"
        xl = 0.0
        xr = 1.0
        it = 0
        while xr - xl > 1e-6 
            xm = (xl + xr) / 2
            ρxm = interpolate_concentration(system, xm)
            if convert_concentration(system, xm, "x" => "mol/L"; density = ρxm) > concentration
                xr = xm
            else
                xl = xm
            end
            it += 1
            if it > 1000
                throw(ArgumentError("Could not find the molar fraction corresponding to the given concentration."))
            end
        end
        xl
    else
        convert_concentration(system, concentration, cunit => "x")
    end

    # Convert concentrations, for checking
    ρ = interpolate_concentration(system, c_x)
    c_vv = convert_concentration(system, concentration, cunit => "vv"; density = ρ)
    cc_mol = convert_concentration(system, concentration, cunit => "mol/L"; density = ρ)

    # aliases for clearer formulas
    ρc = density_pure_cossolvent(system) # of the pure cossolvent

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

    # Solution volume (vbox - vprotein)
    vs = vbox - CMV * Mp / ρ

    # number of cossolvent molecules: cossolvent concentration × volume of the solution
    nc = round(Int, cc * vs)

    #
    # number of water molecules, obtained from the mass difference
    #
    # nc (molecules) / NA (mol) * Mc (g/mol) is the mass of the cossolvent molecules
    # ρ*vs is the total mass of the solution (density (g/L) × volume (L) 
    # ρ*vs - (nc/NA)*Mc is the mass of water in the solution (g)
    # (ρ*vs - (nc/NA)*Mc)/Mw is the number of mols of water in the solution (mol)
    # NA*(ρ*vs - (nc/NA)*Mc)/Mw is the number of water molcules
    # rearranging: nw = (NA*ρ*vs - nc*Mc)/Mw 
    # But since we have vs in Å³, we need the conversion vs = vs / 1e24, and we have
    # nw = (NA*ρ*vs/1e24 - nc*Mc)/Mw 
    # given that CMV = 1e24/NA, we have nw = (ρ*vs/CMV - nc*Mc)/Mw
    nw = round(Int, (ρ * vs / CMV - nc * Mc) / Mw)

    # Final density of the solution
    ρ = CMV * (Mc * nc + Mw * nw) / vs

    # Final cossolvent concentration (mol/L)
    cc_f = 1000 * (nc / vs) * CMV

    # Final water concentration (mol/L)
    cw_f = 1000 * (nw / vs) * CMV

    # Final recovered concentration in vv
    vv = CMV * (nc * Mc / ρc) / vs

    # Half of box sides, to center the solute at the origin
    l = round.(box_sides ./ 2; digits=3)

    println(
        """
        ==================================================================
        Summary:
        ==================================================================

        Target concentration = $cc_mol mol/L
                             = $c_vv volume fraction
                             = $c_x molar fraction
                             = $cc molecules/Å³

        Box volume = $vbox Å³
        Solution volume = $vs Å³   
        Solute extrema = [ $(join(-0.5*solute_extrema, ", ")), $(join(0.5*solute_extrema, ", ")) ] Å
        Box sides = [ $(join( -1.0*l, ", ")), $(join( l, ", ")) ] Å 

        Density = $ρ g/mL
        Solute molar mass = $Mp g/mol
        Cossolvent molar mass = $Mc g/mol
        Solvent molar mass = $Mw g/mol

        Number of cossolvent ($(basename(cossolvent_pdbfile))) molecules = $nc molecules
        Number of solvent ($(basename(solvent_pdbfile))) molecules = $nw molecules

        Final cossolvent concentration = $cc_f mol/L
        Final solvent concentration = $cw_f mol/L
                                    = $(CMC*cw_f) molecules/Å³
                                    
        Final solvent density = $ρ g/mL
        Final volume fraction = $vv 
        Final molar fraction = $(nc/(nc+nw))

        ==================================================================
        """)

    open(input, "w") do io
        println(io,
            """
            # 
            # Packmol input file
            # 
            # Generated by MolSimToolkit.jl
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

