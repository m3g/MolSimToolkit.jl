```@meta
CollapsedDocStrings = true
```
# Contact analysis

Contact analysis routines allow the computation and tracking of contacts along a simulation.

## Distance between two residues

This function computes the minimum distance between the atoms of two residues. It can be 
used standalone, but it is also the basis for the computation of contact maps. 

```@docs
residue_residue_distance
```

## Contact map

The `contact_map` function is used to evaluate contacts within a structure, or among
two structures.

```@docs
contact_map
```

The output of `contact_map` is a `ContactMap` data structure:

```@docs
ContactMap
```



