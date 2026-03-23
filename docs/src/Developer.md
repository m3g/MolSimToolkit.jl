# Developer zone

```@meta
CollapsedDocStrings = true
```

## Simulation

```@autodocs
Modules = [ MolSimToolkit ]
Pages = [ "datastructures/Simulation.jl", "datastructures/Trajectory.jl", "datastructures/Frame.jl" ]
Filter = (f) -> !(nameof(f) === :unitcell)
```

## Positions

```@autodocs
Modules = [ MolSimToolkit ]
Pages = [ "datastructures/Positions.jl" ]
Order = [ :function, :type ]
```

## Unit cell

```@autodocs
Modules = [ MolSimToolkit ]
Pages = [ "datastructures/Frame.jl" ]
Filter = (f) -> (nameof(f) === :unitcell)
```

## Wrap coordinates

```@autodocs
Modules = [ MolSimToolkitShared ]
Pages = [ "wrap.jl" ]
```




