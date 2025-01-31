
# Molecular Minimum Distances

Computes the minimum distance between *molecules*, which are represented as arrays of coordinates in two or three dimensions. 

To understand the utility and purpose of this package, consider the image below:

![nearest.png](./images/molecular_minimum_distances/nearest.png)

Here, there is one *blue* molecule, with 6 atoms, and several *red* molecules, with 2 atoms each. The package has identified which are the molecules of the *red* set that have at leat one atom within a cutoff from the atoms of the *blue* molecule, and annotated the corresponding atoms and the distances.

## Features

- Fast [cell-list approach](https://github.com/m3g/CellListMap.jl), to compute minimum-distance for thousands, or millions of atoms. 
- General periodic boundary conditions supported. 
- Advanced mode for in-place calculations, for non-allocating iterative calls (for analysis of MD trajectories, for example).
- Modes for the calculation of minimum-distances in sets of molecules.

## Most typical use: Understanding solvation

The most typical scenario is that of a protein, or another macromolecule, in a box of solvent. For example, here we download a frame of a protein which was simulated in a mixture of water and TMAO: 

```julia-repl
julia> using MolSimToolkit, PDBTools

julia> atoms = MolSimToolkit.MolecularMinimumDistances.download_example()
   Array{Atoms,1} with 62026 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
                                                       ⋮ 
   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  0.00  1.00     1    WAT2     62024
   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025
   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026
```

Next, we extract the protein coordinates, and the TMAO coordinates:

```julia-repl
julia> protein = coor(atoms,"protein")
1463-element Vector{SVector{3, Float64}}:
 [-9.229, -14.861, -5.481]
 [-10.048, -15.427, -5.569]
 [-9.488, -13.913, -5.295]
 ⋮
 [6.408, -12.034, -8.343]
 [6.017, -10.967, -9.713]

julia> tmao = coor(atoms,"resname TMAO")
2534-element Vector{SVector{3, Float64}}:
 [-23.532, -9.347, 19.545]
 [-23.567, -7.907, 19.381]
 [-22.498, -9.702, 20.497]
 ⋮
 [13.564, -16.517, 12.419]
 [12.4, -17.811, 12.052]
```

The system was simulated with periodic boundary conditions, with sides in this frame of `[83.115, 83.044, 83.063]`, and this information will be provided to the minimum-distance computation.

Finally, we find all the TMAO molecules having at least one atom closer than 12 Angstroms to the protein, using the current package (TMAO has 14 atoms):

```julia-repl
julia> list = minimum_distances(
           xpositions=tmao, # solvent
           ypositions=protein, # solute
           xn_atoms_per_molecule=14,
           cutoff=12.0,
           unitcell=[83.115, 83.044, 83.063]
       )
181-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(true, 2526, 97, 9.652277658666891)

julia> count(x -> x.within_cutoff, list)
33
```

Thus, 33 TMAO molecules are within the cutoff distance from the protein, and the distances can be used to study the solvation of the protein.

!!! tip
    The `coordination_number` function of this package essentially performs the above calculation iteratively 
    along a trajectory. The source of of such function is simple and can be used to further understand the utility
    and usage of the minimum-distance calculations.

## Performance

This package exists because this computation is fast. For example, let us choose the water molecules instead, and benchmark the time required to compute this set of distances:
```julia-repl
julia> water = coor(atoms,"resname TIP3")
58014-element Vector{SVector{3, Float64}}:
 [-28.223, 19.92, -27.748]
 [-27.453, 20.358, -27.476]
 [-27.834, 19.111, -28.148]
 ⋮
 [13.218, -3.647, -34.453]
 [12.618, -4.977, -34.303]

julia> using BenchmarkTools

julia> @btime minimum_distances(
           xpositions=$water, # solvent
           ypositions=$protein, # solute
           xn_atoms_per_molecule=3,
           cutoff=12.0,
           unitcell=[83.115, 83.044, 83.063]
       );
  6.288 ms (3856 allocations: 13.03 MiB)
```

To compare, a naive algorithm to compute the same thing takes roughly 400x more for this
system size:

```julia-repl
julia> @btime MolSimToolkit.MolecularMinimumDistances.naive_md($water, $protein, 3, [83.115, 83.044, 83.063], 12.0);
  2.488 s (97 allocations: 609.16 KiB)
```

And the computation can be made faster and in-place using the more advanced interface that allows preallocation of main necessary arrays:

```julia-repl
julia> sys = CrossPairs(
           xpositions=water, # solvent
           ypositions=protein, # solute
           xn_atoms_per_molecule=3,
           cutoff=12.0,
           unitcell=[83.115, 83.044, 83.063]
       )
CrossPairs system with:

Number of atoms of set: 58014
Number of atoms of target structure: 1463
Cutoff: 12.0
unitcell: [83.12, 0.0, 0.0, 0.0, 83.04, 0.0, 0.0, 0.0, 83.06]
Number of molecules in set: 4144

julia> @btime minimum_distances!($sys);
  2.969 ms (196 allocations: 22.80 KiB)
```

The remaining allocations occur only for the launching of multiple threads:

```julia-repl
julia> sys = CrossPairs(
           xpositions=water, # solvent
           ypositions=protein, # solute
           xn_atoms_per_molecule=14,
           cutoff=12.0,
           unitcell=[83.115, 83.044, 83.063],
           parallel=false # default is true
       );

julia> @btime minimum_distances!($sys);
  15.249 ms (0 allocations: 0 bytes)
```

## Example input files

The examples here use a molecular system, but the package actually only considers the coordinates of the atoms and the number of atoms of each molecule. Thus, more general distance problems can be tackled.

The input atomic positions used in the following examples can be obtained with:

```julia-repl
julia> using MolSimToolkit, PDBTools

julia> system = MolSimToolkit.MolecularMinimumDistances.download_example() 
   Array{Atoms,1} with 62026 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     ALA     A        1        1   -9.229  -14.861   -5.481  0.00  1.00     1    PROT         1
       2  HT1     ALA     A        1        1  -10.048  -15.427   -5.569  0.00  0.00     1    PROT         2
       3  HT2     ALA     A        1        1   -9.488  -13.913   -5.295  0.00  0.00     1    PROT         3
                                                       ⋮ 
   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  0.00  1.00     1    WAT2     62024
   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025
   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026

```
The system consists of a protein (with 1463 atoms), solvated by 181 TMAO molecules (with 14 atoms each), 19338 water molecules, and some ions. 

These coordinates belong to a snapshot of a simulation which was performed with cubic periodic boundary conditions, with a box side of `84.48` Angstrom. 

The coordinates of each of the types of molecules can be extracted from the `system` array of atoms with (using `PDBTools` - `v0.13` or greater):

```julia-repl
julia> protein = coor(system,"protein")
1463-element Vector{StaticArrays.SVector{3, Float64}}:
 [-9.229, -14.861, -5.481]
 [-10.048, -15.427, -5.569]
 [-9.488, -13.913, -5.295]
 ⋮
 [6.408, -12.034, -8.343]
 [6.017, -10.967, -9.713]

julia> tmao = coor(system,"resname TMAO")
2534-element Vector{StaticArrays.SVector{3, Float64}}:
 [-23.532, -9.347, 19.545]
 [-23.567, -7.907, 19.381]
 [-22.498, -9.702, 20.497]
 ⋮
 [13.564, -16.517, 12.419]
 [12.4, -17.811, 12.052]

julia> water = coor(system,"water")
58014-element Vector{StaticArrays.SVector{3, Float64}}:
 [-28.223, 19.92, -27.748]
 [-27.453, 20.358, -27.476]
 [-27.834, 19.111, -28.148]
 ⋮
 [13.218, -3.647, -34.453]
 [12.618, -4.977, -34.303]
```

Using these vectors of coordinates, we will illustrate the use of the current package.

##  Shortest distances from a solute

The simplest usage consists of finding for each molecule of one set the atoms of the other set which are closer to them. For example, here we want the atoms of the proteins which are closer to each TMAO molecule (14 atoms), within a cutoff of `12.0` Angstroms.

The simulations was performed with periodic boundary conditions, in a cubic box of sides `[84.48, 84.48, 84.48]`. We compute the minimum distances with:

```julia-repl
julia> list = minimum_distances(
           xpositions=tmao, # solvent
           ypositions=protein, # solute
           xn_atoms_per_molecule=14,
           cutoff=12.0,
           unitcell=[84.48, 84.48, 84.48]
       )
181-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(true, 2526, 97, 9.652277658666891)
```

The `list` contains, for each *molecule* of TMAO, a `MinimumDistance` object, containing the following fields, 
exemplified by printing the last entry of the list:
```julia-repl
julia> list[end]
MinimumDistance{Float64}(true, 2526, 97, 9.652277658666891)

Distance within cutoff, within_cutoff = true
x atom of pair, i = 2526
y atom of pair, j = 97
Distance found, d = 9.652277658666891
```

The fields `within_cutoff`, `i`, `j`, and `d` show if a distance was found within the cutoff,
the indices of the atoms involved in the contact, and their distance.

Getter functions are available to extract eac hof these fields, to add some convenience:
`within_cutoff`, `iatom`, `jatom`, and `distance`.

!!! note
    If the solute has more than one molecule, this will not be taken into 
    consideration in this mode. All molecules will be considered as part
    of the same structure (the number of atoms per molecule of the `protein` is not a parameter here).

## All shortest distances

A similar call of the previous section can be used to compute, for each molecule of a set of molecules, which is the closest atom
of every other molecule of another set. 

In the example, we can compute for each TMAO molecule, which is the closest atom of water, and vice-versa. The difference from the previous call
is that now wee need to provide the number of atoms of both TMAO and water:

```julia-repl
julia> water_list, tmao_list = minimum_distances(
           xpositions=water,
           ypositions=tmao,
           xn_atoms_per_molecule=3,
           yn_atoms_per_molecule=14,
           unitcell=[84.48, 84.48, 84.48],
           cutoff=12.0
       );

julia> water_list
19338-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 2, 1512, 4.779476331147592)
 MinimumDistance{Float64}(true, 6, 734, 2.9413928673334357)
 MinimumDistance{Float64}(true, 8, 859, 5.701548824661595)
 ⋮
 MinimumDistance{Float64}(true, 58010, 1728, 3.942870781549911)
 MinimumDistance{Float64}(true, 58014, 2058, 2.2003220218867936)

julia> tmao_list
181-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 12, 22520, 2.1985345118965056)
 MinimumDistance{Float64}(true, 20, 33586, 2.1942841657360606)
 MinimumDistance{Float64}(true, 37, 26415, 2.1992319113726926)
 ⋮
 MinimumDistance{Float64}(true, 2512, 37323, 2.198738501959709)
 MinimumDistance{Float64}(true, 2527, 33664, 2.1985044916943015)
```

Two lists were returned, the first containing, for each water molecule, `MinimumDistance` data associated to the closest TMAO molecule
(meaning the atoms involved in the contact and their distance). Similarly, the second list contains, for each TMAO molecule, the `MinimumDistance` data associated to each TMAO molecule. 

## Shortest distances within molecules

There is an interface to compute the shortest distances of molecules within a set of molecules. That is, given one group of molecules, compute for each molecule which is the shortest distance among the other molecules of the same type. 

A typical call would be:

```julia-repl
julia> water_list = minimum_distances(
           xpositions=water,
           xn_atoms_per_molecule=3,
           unitcell=[84.48, 84.48, 84.48],
           cutoff=12.0
       )
19338-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 2, 33977, 2.1997806708851724)
 MinimumDistance{Float64}(true, 4, 43684, 2.1994928961012814)
 MinimumDistance{Float64}(true, 9, 28030, 2.1997583958244142)
 ⋮
 MinimumDistance{Float64}(true, 58010, 22235, 2.1992096307537414)
 MinimumDistance{Float64}(true, 58012, 9318, 2.20003227249056)
```

Which contains for each water molecule the atoms involved in the closest contact to any other water molecule, and the distances (within the cutoff).
A pictorial representation of a result of this type is, for a simpler system:

![self pairs](./images/molecular_minimum_distances/self_pair.png)

This can be used for the identification of connectivity networks, for example, or for some types of clustering.

## Advanced usage

### System build and update

If the molecular minimum distances will be computed many times for similar systems, it is possible
to construct the system and update its properties. The use of the interface of `CellListMap`
is required (requires `CellListMap` version `0.7.24` or greater). 

For example, let us build one system with a protein and water:

```julia-repl
julia> using MolSimToolkit, PDBTools

julia> system = MolecularMinimumDistances.download_example();

julia> protein = coor(system, "protein");

julia> water = coor(system, "water");
```

We now build the `CrossPairs`  type of system, instead of calling the `minimum_distances` function directly:

```julia-repl
julia> sys = CrossPairs(
           xpositions=water, # solvent
           ypositions=protein, # solute
           xn_atoms_per_molecule=3,
           cutoff=12.0,
           unitcell=[84.48, 84.48, 84.48]
       )
CrossPairs system with:

Number of atoms of set x: 58014
Number of molecules in set x: 19338
Number of atoms of target structure y: 1463
Cutoff: 12.0
unitcell: [84.48, 0.0, 0.0, 0.0, 84.48, 0.0, 0.0, 0.0, 84.48]
```

Now `sys`  contains the necessary arrays for computing the list of minimum distances. We use now the
`minimum_distances!`  function (with the `!`), to update that list:

```julia-repl
julia> minimum_distances!(sys)
19338-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(true, 58011, 383, 10.24673074692606)
 MinimumDistance{Float64}(false, 0, 0, Inf)
```

The system can be now updated: the positions, cutoff, or unitcell can be modified, with the 
following interfaces:

### Updating positions

To update the positions, modify the `sys.xpositions` (or `ypositions`)  array. We will
boldy demonstrate this by making the first atom of the `x` set to be close to the first
atom of the protein, and recomputing the distances:

```julia-repl
julia> using StaticArrays

julia> sys.xpositions[2] = sys.ypositions[1] + SVector(1.0,0.0,0.0);

julia> minimum_distances!(sys)
19338-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 2, 4, 0.9202923448556931)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(true, 58011, 383, 10.24673074692606)
 MinimumDistance{Float64}(false, 0, 0, Inf)
```

### Updating the cutoff, unitcell and parallel flag

The `cutoff`, `unitcell` and `parallel` data of the `sys` objects can be modified 
directly. For example:
```julia-repl
julia> sys
CrossPairs system with:

Number of atoms of set x: 58014
Number of molecules in set x: 19338
Number of atoms of target structure y: 1463
Cutoff: 15.0
unitcell: [100.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 100.0]

julia> sys.cutoff = 10.0
10.0

julia> sys.unitcell = [84.4, 84.4, 84.4]
3-element Vector{Float64}:
 84.4
 84.4
 84.4

julia> sys.parallel = false
false

julia> sys
CrossPairs system with:

Number of atoms of set x: 58014
Number of molecules in set x: 19338
Number of atoms of target structure y: 1463
Cutoff: 10.0
unitcell: [84.4, 0.0, 0.0, 0.0, 84.4, 0.0, 0.0, 0.0, 84.4]
```

!!! note
    It is not possible to update the `unitcell` from a Orthorhombic to a general Triclinic cell. If the system
    will be Triclinic at any moment, the `unitcell` must be initialized with the full matrix instead of a 
    vector of sides.

### Index of molecules

Additionally, the low level interface allows the definition of more general groups of particles, in the sense that "molecule" can have different number of atoms in the same set. Therefore, one needs to provide *a function* that returns the index of the molecule of each atom, given the index of the atom. 

Briefly, if a set of atoms belong to molecules of the same number of atoms, one can compute the index of each molecule using
```julia
mol_indices(i,n) = div((i - 1), n) + 1
```
where `i` is the atom index in the array of coordinates, and `n` is the number of atoms per molecule. This is the default assumed in the basic interface, and can be called with:
```julia-repl
julia> using StaticArrays

julia> x = rand(SVector{3,Float64},9); # 3 water molecules

julia> mol_indices(2,3) # second atom belongs to first molecule
1

julia> mol_indices(4,3) # fourth atom belongs to second molecule
2
```

Typically, as we will show, this function will be used for setting up molecule indices.

However, more general indexing can be used. For instance, let us suppose that the 9 atoms of the `x` array of coordinates above belong to `2` molecules, with `4` and `5` atoms each. Then, we could define, for example:

```julia-repl
julia> my_mol_indices(i) = i <= 4 ? 1 : 2
my_mol_indices (generic function with 1 method)

julia> my_mol_indices(4)
1

julia> my_mol_indices(5)
2
```

Since the function can close-over an array of molecular indices, the definition can be completely general, that is:

```julia-repl
julia> molecular_indices = [ 1, 3, 3, 2, 2, 1, 3, 1, 2 ];

julia> my_mol_indices(i) = molecular_indices[i]
my_mol_indices (generic function with 1 method)

julia> my_mol_indices(1)
1

julia> my_mol_indices(5)
2
```

In summary, this function that given the index of the atom returns the index of the corresponding molecule must be provided in the advanced interface, and typically will be just a closure around the number of atoms per molecule, using the already available `mol_indices` function. 

#### Example

Let us mix water and TMAO molecules in the same set, and use a general function to compute the indices of the molecules of each atom: 

```julia-repl
julia> system = MolecularMinimumDistances.download_example();

julia> protein = coor(system, "protein");

julia> tmao_and_water = select(system, "resname TMAO or resname TIP3")
   Array{Atoms,1} with 60548 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
    1479    N    TMAO     A        1      120  -23.532   -9.347   19.545  0.00  1.00     1    TMAO      1479
    1480   C1    TMAO     A        1      120  -23.567   -7.907   19.381  0.00  1.00     1    TMAO      1480
    1481   C2    TMAO     A        1      120  -22.498   -9.702   20.497  0.00  1.00     1    TMAO      1481
                                                       ⋮ 
   62024  OH2    TIP3     C     9339    19638   13.485   -4.534  -34.438  0.00  1.00     1    WAT2     62024
   62025   H1    TIP3     C     9339    19638   13.218   -3.647  -34.453  0.00  1.00     1    WAT2     62025
   62026   H2    TIP3     C     9339    19638   12.618   -4.977  -34.303  0.00  1.00     1    WAT2     62026

julia> findfirst(at -> at.resname == "TIP3", tmao_and_water)
2535
```

Thus, the `tmao_and_water` atom array has two different types of molecules, TMAO with 14 atoms, and water with 3 atoms. 
The first atom of a water molecule is atom `2535` of the array. We extract the coordinates of the atoms with:
```julia-repl
julia> solvent = coor(tmao_and_water)
60548-element Vector{SVector{3, Float64}}:
 [-23.532, -9.347, 19.545]
 [-23.567, -7.907, 19.381]
 [-22.498, -9.702, 20.497]
 ⋮
 [13.218, -3.647, -34.453]
 [12.618, -4.977, -34.303]
```

 And now we define a function that, given the index of the atom, returns the molecule to which it belongs:

```julia-repl
julia> function mol_indices(i) 
           if i < 2535 # TMAO (14 atoms per molecule) 
               div(i-1,14) + 1 
           else # water (3 atoms per molecule)
               mol_indices(2534) + div(i-2534-1,3) + 1
           end
       end
mol_indices (generic function with 3 method)
```

The function above computes the molecular indices for TMAO in the standard way, and computes the water 
molecular indices by first summing the molecule index of the last TMAO molecule, and subtracting from the
atomic index of water the last index of the last TMAO atom. We can test this: 

```julia-repl
julia> mol_indices(14) # last atom of first TMAO
1

julia> mol_indices(15) # first atom of second TMAO
2

julia> mol_indices(2534) # last atom of last TMAO
181

julia> mol_indices(2535) # first atom of first water
182

julia> mol_indices(2537) # last atom of first water
182

julia> mol_indices(2538) # first atom of second water
183
```

With this function, we can construct the system using it instead of the `xn_atoms_per_molecule` integer variable,
to obtain the solvation of the protein by both TMAO and water in a single run:

```julia-repl
julia> sys = CrossPairs(
           xpositions=solvent, # solvent = coor(tmao_and_water)
           ypositions=protein, # solute
           xmol_indices = mol_indices,
           cutoff=12.0,
           unitcell=[84.48, 84.48, 84.48]
       )
CrossPairs system with:

Number of atoms of set x: 60548
Number of molecules in set x: 19519
Number of atoms of target structure y: 1463
Cutoff: 12.0
unitcell: [84.48, 0.0, 0.0, 0.0, 84.48, 0.0, 0.0, 0.0, 84.48]
```

As we can see, the number of molecules is correct (the sum of the number of water and tmao
molecules). And the list of minimum distances will retrive the information of the closest
protein atom to all solvent molecules of the set:

```julia-repl
julia> minimum_distances!(sys)
19519-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(true, 60545, 383, 10.24673074692606)
 MinimumDistance{Float64}(false, 0, 0, Inf)
```

## Citation

If this package was useful, please cite the article describing the main algorithms on which it is based:

L. Martínez, **CellListMap.jl: Efficient and customizable cell list implementation for calculation of pairwise particle properties within a cutoff.**
*Computer Physics Communications* **279**, 108452 (2022). 

DOI: [10.1016/j.cpc.2022.108452](https://doi.org/10.1016/j.cpc.2022.108452)

## Help entries

```@autodocs
Modules = [MolSimToolkit.MolecularMinimumDistances]
Order = [:function, :type]
```























## Details of the illustration

The initial illustration here consists of a toy solute-solvent example, where the solute is an approximately hexagonal molecule, and the solvent is composed by 40 diatomic molecules. The toy system is built as follows:

```julia
using MolecularMinimumDistances, StaticArrays
# x will contain the "solvent", composed by 40 diatomic molecules
T = SVector{2,Float64}
x = T[]
cmin = T(-20,-20)
for i in 1:40
    v = cmin .+ 40*rand(T)
    push!(x, v)
    theta = 2pi*rand()
    push!(x, v .+ T(sin(theta),cos(theta)))
end
# y will contain the "solute", composed by an approximate hexagonal molecule
y = [ T(1,1), T(1,-1), T(0,-1.5), T(-1,-1), T(-1,1), T(0,1.5) ]
```

Next, we compute the minimum distances between each molecule of `x` (the solvent)
and the solute. In the input we need to specify the number of atoms of each molecule
in `x`, and the cutoff up to which we want the distances to be computed:

```julia-repl
julia> list = minimum_distances(
           xpositions=x,
           ypositions=y,
           xn_atoms_per_molecule=2,
           unitcell=[40.0, 40.0],
           cutoff=10.0
       )
40-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 2, 3, 1.0764931248364737)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 ⋮
 MinimumDistance{Float64}(true, 74, 5, 7.899981412729262)
 MinimumDistance{Float64}(false, 0, 0, Inf)
```

The output is a list of `MinimumDistance` data structures, one for each molecule in `x`. The `true` indicates that a distance smaller than the cutoff was found, and for these the indices of the atoms in `x` and `y` associated are reported, along with the distance between them.

In this example, from the 40 molecules of `x`, eleven had atoms closer than the cutoff to some
atom of `y`:
```julia-repl
julia> count(x -> x.within_cutoff, list)
11
```

We have an auxiliary function to plot the result, in this case where the "atoms" are bi-dimensional:

```julia
using Plots
import MolecularMinimumDistances: plot_md!
p = plot(lims=(-20,20),framestyle=:box,grid=false,aspect_ratio=1)
plot_md!(p, x, 2, y, 6, list, y_cycle=true)
```
will produce the illustration plot above, in which the nearest point between the two sets is identified.

