using PDBTools

export mvalue, parse_mvalue_server_sasa, run_gmx_sasa

abstract type CosolventType end
struct Urea <: CosolventType end

"""
    mvalue(; pdbname, sasas, type=1)

Calculates the m-value transfer free energy of a protein in 1M urea solution using the Tanford transfer model,
as implemented by Moeser and Horinek (https://doi.org/10.1021/acs.jpcb.7b02138).

# Arguments

- `pdbname::AbstractString`: Path to the PDB file of the protein structure.
- `sasas::Dict{String, Dict{Symbol, Float64}}`: A dictionary containing the change in solvent accessible surface area (SASA)
  upon denaturation for each amino acid type. This data can be obtained from the m-value server or calculated using GROMACS:
    - The output of the server can be parsed using the `parse_mvalue_server_sasa` function defined in this module.
    - Alternatively, the SASA values can be calculated using GROMACS with the `run_gmx_sasa` function defined in this module.
       This function requires GROMACS to be installed and accessible from the command line. 

- `type::Int`: Specifies which SASA value to use from the provided data, because the server provides minimum, average, and maximum values,
    according to different denatured models for the protein. The recommended value is `2` for comparison with experimental data.
    Normally, GROMACS calculations will provide a single value, so `type=1` should be used in that case.

# Returns

A named tuple with the following fields:
- `tot`: Total transfer free energy (kcal/mol).
- `bb`: Contribution from the backbone (kcal/mol).
- `sc`: Contribution from the side chains (kcal/mol).
- `restype`: A dictionary with the transfer free energy contributions per residue type.

Each entry in the dictionary is a named tuple with `bb` and `sc` fields representing the backbone and side chain contributions, respectively.

# Example calls

```julia
# Using SASA values from the m-value server
sasas_from_server=parse_mvalue_server_sasa(server_output)
mvalue(; pdbname="protein.pdb", sasas=sasas_from_server, type=2)

# Using SASA values calculated with GROMACS
sasas_gmx=run_gmx_sasa(native_pdb="native.pdb", desnat_pdb="desnat.pdb")
mvalue(; pdbname="protein.pdb", sasas=sasas_gmx, type=1)
```

"""
function mvalue(::Type{<:CosolventType}=Urea; pdbname, sasas, type=1)
    protein = read_pdb(pdbname, "protein")
    residue_types = unique(resname.(protein)) 
    DeltaG_per_residue = Dict{String, @NamedTuple{bb::Float64, sc::Float64}}()
    for rname in residue_types
        DeltaG_per_residue[rname] = (bb = (last(tfe_asa(Urea, rname))) * (sasas[rname][:bb][type]/100), 
                                     sc = (first(tfe_asa(Urea, rname))) * (sasas[rname][:sc][type]/100))
    end
    DeltaG_BB = sum(getfield(DeltaG_per_residue[key], :bb) for key in keys(DeltaG_per_residue))
    DeltaG_SC = sum(getfield(DeltaG_per_residue[key], :sc) for key in keys(DeltaG_per_residue))
    DeltaG = DeltaG_BB + DeltaG_SC
    return (tot=DeltaG, bb=DeltaG_BB, sc=DeltaG_SC, restype=DeltaG_per_residue)
end

#=
    tfe_asa(::Type{Urea}, restype::AbstractString)

Returns the transfer free energy per unit area (kcal/nm^2) for side chain and backbone
for a given amino acid type in urea solution according to the Tanford transfer model.

Isolated ASA values are from the Supporting Table 2 of https://doi.org/10.1073/pnas.0507053102
(https://www.pnas.org/doi/suppl/10.1073/pnas.0507053102/suppl_file/07053table2.pdf)

The dg_aa values are GTFE+ values from Table S1 of https://doi.org/10.1021/jp409934q
(https://pubs.acs.org/doi/10.1021/jp409934q#_i14)

The transfer free energy per unit area is calculated as:

    Δg = (dg_aa / 1000) / (ASA / 100)

where dg_aa is in cal/mol and ASA in Å², resulting in kcal/nm².

The backbone contribution is assumed to be the same for all amino acids and is calculated
using the backbone ASA of glycine, according to the universal model of Moeser and Horinek.

=#
function tfe_asa(::Type{Urea}, restype::AbstractString)
    isolated_ASA = Dict{String,Tuple{Float64,Float64}}( 
                   # BB      SC   (Å^2)
        "ALA"	=> (46.2,	71.9),
        "PHE"	=> (38.4,	184.4),
        "LEU"	=> (35.3,	157.8),
        "ILE"	=> (30.9,	150.1),
        "VAL"	=> (36.1,	128.4),
        "PRO"	=> (35.6,	111.0),
        "MET"	=> (38.6,	164.8),
        "TRP"	=> (37.4,	228.9),
        "GLY"	=> (88.1,	0),
        "SER"	=> (44.0,	85.8),
        "THR"	=> (37.9,	114.6),
        "TYR"	=> (38.7,	198.1),
        "GLN"	=> (37.8,	155.4),
        "ASN"	=> (40.2,	125.3),
        "ASP"	=> (40.5,	118.2),
        "GLU"	=> (37.8,	148.4),
        "HIS"	=> (40.4,	162.1),
        "LYS"	=> (38.7,	187.1),
        "ARG"	=> (39.1,	216.9),
        "CYS"	=> (42.6,	103.5),
    ) 

    dg_aa = Dict( # cal/mol (1M)
        "BB" => -39,
        "ALA" => 1.01,
        "PHE" => -68.64,
        "LEU" => -40.10,
        "ILE" => -23.96,
        "VAL" => -7.18,
        "PRO" => -3.18,
        "MET" => -33.87,
        "TRP" => -126.99,
        "GLY" => 0.00,
        "SER" => -6.09,
        "THR" => -7.62,
        "TYR" => -30.61,
        "GLN" => -40.34,
        "ASN" => -24.32,
        "ASP" => 18.02,
        "GLU" => 15.09,
        "HIS" => -36.04,
        "HSD" => -36.04,
        "HSE" => -36.04,
        "LYS" => -8.29,
        "ARG" => -6.70,
        "CYS" => 0.00,
    )

    # converted to kcal/nm^2 here
    bb_contribution = (dg_aa["BB"]/1000) / (first(isolated_ASA["GLY"])/100) # universal model: equal for all AAs
    if restype == "GLY"
        return 0.0, bb_contribution
    end
    return (dg_aa[restype]/1000) / (last(isolated_ASA[restype])/100), bb_contribution # side chain
end

#=
    read_gmx_sasa_values(filename::String, n)

Reads the output of `gmx sasa` and returns the SASA values.
`n` is the number of surfaces calculated (1 for BB only, 2 for SC and BB, for example).

=#
function read_gmx_sasa_values(filename::String, n)
    local sasa_values
    open(filename, "r") do io
        for line in eachline(io)
            if startswith(line, r"@|#")
                continue
            end
            sasa_values = parse.(Float64, split(line)[3:2+n])
        end
    end
    return sasa_values
end

"""
    parse_mvalue_server_sasa(string::AbstractString)

Parses the SASA output from the m-value calculator server (http://best.bio.jhu.edu/mvalue/), into a dictionary
that can be directly used as input to the `mvalue` function.

The input string should contain lines formatted as follows, and correspond to the SASA values for each amino acid type:

```julia
ALA 		    8 	 (    11.1)     79.1 [   147.1] 	 | 	 (   -13.0)     51.4 [   115.8] 
PHE 		    3 	 (   166.9)    197.1 [   230.2] 	 | 	 (    29.4)     56.4 [    83.4] 
LEU 		    7 	 (   475.2)    532.2 [   589.3] 	 | 	 (    89.3)    145.3 [   201.3] 
...
LYS 		    6 	 (   171.5)    220.4 [   269.3] 	 | 	 (    -4.5)     42.0 [    88.5] 
ARG 		    1 	 (   110.2)    124.4 [   138.6] 	 | 	 (    17.1)     25.0 [    33.0] 
CYS 		    0 	 (     0.0)      0.0 [     0.0] 	 | 	 (     0.0)      0.0 [     0.0] 
```

This data can be found in the output of the server, under the title "Sidechain and Backbone changes in Accessible Surface Area."

The function returns a dictionary where each key is an amino acid three-letter code (e.g., "ALA", "PHE"), and the value 
is another dictionary with two keys: `:sc` for side chain SASA values and `:bb` for backbone SASA values. 
Each of these keys maps to a tuple containing three Float64 values representing the minimum, average, and maximum SASA values in Å².

"""
function parse_mvalue_server_sasa(string::AbstractString)
    sasa = Dict{String, Dict{Symbol, Tuple{Float64,Float64,Float64}}}()
    for line in split(string, "\n")
        # Replace (, ) and [, ], | with spaces
        line = replace(line, '(' => ' ', ')' => ' ', '[' => ' ', ']' => ' ', '|' => ' ')
        data = split(line)
        if length(data) < 8
            continue
        end
        resname = data[1]
        sc1 = parse(Float64, data[3])
        sc2 = parse(Float64, data[4])
        sc3 = parse(Float64, data[5])
        bb1 = parse(Float64, data[6]) 
        bb2 = parse(Float64, data[7])
        bb3 = parse(Float64, data[8])
        sasa[resname] = Dict(:sc => (sc1, sc2, sc3), :bb => (bb1, bb2, bb3))
    end
    return sasa
end

"""
    run_gmx_sasa(; native_pdb::AbstractString, desnat_pdb::AbstractString)

Calculates the change in solvent accessible surface area (SASA) upon denaturation for each amino acid type
using GROMACS. Returns a dictionary that can be directly used as input to the `mvalue` function.

# Arguments

- `native_pdb::AbstractString`: Path to the PDB file of the native protein structure.
- `desnat_pdb::AbstractString`: Path to the PDB file of the denatured protein structure.

# Returns

A dictionary where each key is an amino acid three-letter code (e.g., "ALA", "PHE"), and the value
is another dictionary with two keys: `:sc` for side chain SASA values and `:bb` for backbone SASA values.
Each of these keys maps to a tuple containing a single Float64 value representing the change in SASA upon denaturation in Å².

"""
function run_gmx_sasa(;
    native_pdb::AbstractString,
    desnat_pdb::AbstractString,
)
    sasas = Dict{String, Dict{Symbol, Float64}}()
    p = read_pdb(native_pdb, "protein")
    for rname in unique(resname.(eachresidue(p)))
        sasa_bb_native, sasa_sc_native = 100 .* run_gmx_sasa_single(native_pdb, rname) # returns in nm^2
        sasa_bb_desnat, sasa_sc_desnat = 100 .* run_gmx_sasa_single(desnat_pdb, rname)
        sasas[rname] = Dict(:sc => sasa_sc_desnat - sasa_sc_native, :bb => sasa_bb_desnat - sasa_bb_native)
    end
    return sasas # returns in Å^2
end

#
# Runs gmx sasa for a single residue type in a given PDB file.
#
function run_gmx_sasa_single(pdbname, resname)
    index_file = tempname()*".ndx"
    sasa_file = tempname()*".xvg"
    p = read_pdb(pdbname, "protein and not element H")
    inds_protein = index.(p)
    inds_sidechain = index.(select(p, "resname $resname and not (name N CA C O)"))
    inds_backbone = index.(select(p, "resname $resname and name N CA C O"))
    open(index_file, "w") do io
        write(io, "[ SC ]\n")
        for i in inds_sidechain
            write(io, "$i\n")
        end
        write(io, "[ BB ]\n")
        for i in inds_backbone
            write(io, "$i\n")
        end
        write(io, "[ PROT ]\n")
        for i in inds_protein
            write(io, "$i\n")
        end
    end
    if length(inds_sidechain) == 0 # special case for GLY
        sasa_sc = 0.0
        try
            run(pipeline(`gmx sasa -s $pdbname -probe 0.14 -ndots 500 -surface PROT -output BB -n $index_file -o $sasa_file`, stdout=devnull, stderr=devnull))
        catch
            error("Error running gmx sasa for of $resname in $pdbname")
        end
        sasa_bb = read_gmx_sasa_values(sasa_file, 1)[1]
    else
        try
            run(pipeline(`gmx sasa -s $pdbname -probe 0.14 -ndots 500 -surface PROT -output SC BB -n $index_file -o $sasa_file`, stdout=devnull, stderr=devnull))
        catch
            error("Error running gmx sasa for of $resname in $pdbname")
        end
        sasa_temp = read_gmx_sasa_values(sasa_file, 2)
        sasa_sc, sasa_bb = sasa_temp[1], sasa_temp[2]
    end
    rm(sasa_file, force=true)
    rm(index_file, force=true)
    return sasa_bb, sasa_sc
end
