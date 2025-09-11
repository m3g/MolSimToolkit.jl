using PDBTools

export mvalue, parse_mvalue_server_sasa, run_gmx_sasa
#export sasa_desnat_average

abstract type CosolventType end
struct Urea <: CosolventType end
struct TMAO <: CosolventType end

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
function mvalue(Cosolvent::Type{<:CosolventType}=Urea; pdbname, sasas, type=1)
    protein = read_pdb(pdbname, "protein")
    residue_types = unique(resname.(protein)) 
    DeltaG_per_residue = Dict{String, @NamedTuple{bb::Float64, sc::Float64}}()
    for rname in residue_types
        DeltaG_per_residue[rname] = (bb = (last(tfe_asa(Cosolvent, rname))) * (sasas[rname][:bb][type]/100), 
                                     sc = (first(tfe_asa(Cosolvent, rname))) * (sasas[rname][:sc][type]/100))
    end
    DeltaG_BB = sum(getfield(DeltaG_per_residue[key], :bb) for key in keys(DeltaG_per_residue))
    DeltaG_SC = sum(getfield(DeltaG_per_residue[key], :sc) for key in keys(DeltaG_per_residue))
    DeltaG = DeltaG_BB + DeltaG_SC
    return (tot=DeltaG, bb=DeltaG_BB, sc=DeltaG_SC, restype=DeltaG_per_residue)
end

function _tfe_sasa(restype, tfe_bb, tfe_sc_restype, isolated_ASA_bb, isolated_ASA_sc)
    # converted to kcal/nm^2 here
    bb_contribution = (tfe_bb/1000) / (first(isolated_ASA_bb)/100) 
    if restype == "GLY"
        return 0.0, bb_contribution
    end
    return (tfe_sc_restype/1000) / (last(isolated_ASA_sc)/100), bb_contribution # side chain, backbone
end

#=
    tfe_asa(::Type{Urea}, restype::AbstractString)

Returns the transfer free energy per unit area (kcal/nm^2) for side chain and backbone
for a given amino acid type in urea solution according to the Tanford transfer model,
as implemente by Moeser and Horinek (https://pubs.acs.org/doi/10.1021/jp409934q#_i14).

=#
function tfe_asa(::Type{Urea}, restype::AbstractString;
    tfe_sc_bb = tfe_sc_bb_moeser_and_horinek,
    isolated_ASA = isolated_ASA,
)
    col = 7 # for urea in data tables 
    return _tfe_sasa(
        restype,
        tfe_sc_bb["BB"][col],
        tfe_sc_bb[restype][col],
        isolated_ASA["GLY"], # universal model: same for all residues
        isolated_ASA[restype],
    )
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

#
# Data section
#

#=

Amino acid side-chain and peptide backbone unit transfer free energies (cal/mol) from water to 1M osmolyte
Values for Urea GTFE+ values from Table S1 of https://doi.org/10.1021/jp409934q (https://pubs.acs.org/doi/10.1021/jp409934q#_i14)

=#
const tfe_sc_bb_moeser_and_horinek = Dict(
#                TMAO   Sarcosine     Betaine     Proline    Sorbitol    Sucrose         Urea
    "ALA" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,       1.01),
    "PHE" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -68.64),
    "LEU" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -40.10),
    "ILE" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -23.96),
    "VAL" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -7.18),
    "PRO" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -3.18),
    "MET" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -33.87),
    "TRP" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,    -126.90),
    "GLY" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,       0.00),
    "SER" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -6.09),
    "THR" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -7.62),
    "TYR" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -30.61),
    "GLN" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -40.34),
    "ASN" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -24.32),
    "ASP" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      18.02),
    "GLU" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      15.09),
    "HIS" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -36.04),
    "HSD" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -36.04),
    "HSE" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,     -36.04),
    "LYS" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -8.29),
    "ARG" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,      -6.70),
    "CYS" => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,       0.00),
    "BB"  => (   0.00,       0.00,       0.00,       0.00,       0.00,       0.00,        -39),
)

#=

Amino acid side-chain and peptide backbone unit transfer free energies (cal/mol) from water to 1M osmolyte
Supplementary Table 1 of https://doi.org/10.1073/pnas.0507053102

=#
const tfe_sc_bb_auton_and_bolen = Dict(
#                TMAO   Sarcosine     Betaine     Proline    Sorbitol    Sucrose         Urea
    "ALA" => ( -14.64,      10.91,       4.77,      -0.07,      16.57,      22.05,      -4.69),
    "PHE" => (  -9.32,     -12.64,    -112.93,     -71.26,      26.38,     -96.35,     -83.11),
    "LEU" => (  11.62,      38.33,     -17.73,       4.77,      39.07,      37.11,     -54.57),
    "ILE" => ( -25.43,      39.98,      -1.27,      -2.72,      36.90,      28.12,     -38.43),
    "VAL" => (  -1.02,      29.32,     -19.63,       7.96,      24.65,      33.92,     -21.65),
    "PRO" => (-137.73,     -34.23,    -125.16,     -63.96,      -4.48,     -73.02,     -17.65),
    "MET" => (  -7.65,       8.18,     -14.16,     -35.12,      20.97,      -6.66,     -48.34),
    "TRP" => (-152.87,    -113.03,    -369.93,    -198.37,     -67.23,    -215.27,    -141.46),
    "GLY" => (      0,          0,          0,          0,          0,          0,          0),
    "SER" => ( -39.04,     -27.98,     -41.85,     -33.49,      -1.58,      -2.79,     -20.56),
    "THR" => (   3.57,      -7.54,       0.33,     -18.33,      13.20,      20.82,     -22.09),
    "TYR" => (-114.32,     -26.37,    -213.09,    -138.41,     -53.50,     -78.41,     -45.08),
    "GLN" => (  41.41,     -10.19,       7.57,     -32.26,     -23.98,     -40.87,     -54.81),
    "ASN" => (  55.69,     -40.93,      33.17,     -17.71,     -21.21,     -28.28,     -38.79),
    "ASP" => ( -66.67,     -14.20,    -116.56,     -90.51,     -83.88,     -37.17,       3.55),
    "GLU" => ( -83.25,     -12.61,    -112.08,     -89.17,     -70.05,     -41.65,       0.62),
    "HIS" => (  42.07,     -20.80,     -35.97,     -45.10,     -42.45,    -118.66,     -50.51),
    "LYS" => (-110.23,     -27.42,    -171.99,     -59.87,     -32.47,     -39.60,     -22.76),
    "ARG" => (-109.27,     -32.24,    -109.45,     -60.18,     -24.65,     -79.32,     -21.17),
    "CYS" => (      0,          0,          0,          0,          0,          0,          0), # not reported
    "BB"  => (     90,         52,         67,         48,         35,         62,        -39),
)

#= 

Isolated ASA values are from the Supporting Table 2 of https://doi.org/10.1073/pnas.0507053102
(https://www.pnas.org/doi/suppl/10.1073/pnas.0507053102/suppl_file/07053table2.pdf)

=#
const isolated_ASA = Dict{String,Tuple{Float64,Float64}}( 
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

#=

Average SASA values from lower and upper bounds of the denatured state ensemble,
reported in Supplementary Table 2 of https://doi.org/10.1073/pnas.0507053102

These are the average values of Table 1 of https://doi.org/10.1021/bi962819o

Can be used for testing, but **do not** provide the same accuracy as the values
calculated with GROMACS or obtained from the m-value server.

=# 
const sasa_desnat_average = Dict(
    "ALA" => Dict(:bb => 27.9,	:sc => 55.1),
    "PHE" => Dict(:bb => 24.3,	:sc => 128.8),
    "LEU" => Dict(:bb => 22.7,	:sc => 109.6),
    "ILE" => Dict(:bb => 20.0,	:sc => 117.1),
    "VAL" => Dict(:bb => 20.4,	:sc => 96.4),
    "PRO" => Dict(:bb => 22.5,	:sc => 87.0),
    "MET" => Dict(:bb => 25.3,	:sc => 122.4),
    "TRP" => Dict(:bb => 23.6,	:sc => 156.6),
    "GLY" => Dict(:bb => 65.2,	:sc => 0.0),
    "SER" => Dict(:bb => 29.4,	:sc => 66.5),
    "THR" => Dict(:bb => 24.1,	:sc => 84.3),
    "TYR" => Dict(:bb => 25.6,	:sc => 141.7),
    "GLN" => Dict(:bb => 25.3,	:sc => 116.9),
    "ASN" => Dict(:bb => 25.2,	:sc => 90.1),
    "ASP" => Dict(:bb => 26.0,	:sc => 87.0),
    "GLU" => Dict(:bb => 25.7,	:sc => 113.4),
    "HIS" => Dict(:bb => 24.2,	:sc => 111.5),
    "LYS" => Dict(:bb => 26.1,	:sc => 150.7),
    "ARG" => Dict(:bb => 25.1,	:sc => 171.1),
    "CYS" => Dict(:bb => 26.4,	:sc => 73.0),
)
