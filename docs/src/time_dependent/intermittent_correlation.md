# Correlation functions

```@meta
CollapsedDocStrings = true
```

!!! tip
    `intermittent_correlation` also applies to the identities of hydrogen
    bonds, as computed by [`hydrogen_bond_occupancy`](@ref) — see
    [Persistence of hydrogen bonds](@ref hbond-persistence).

!!! tip
    `intermittent_correlation_profile` resolves the correlation function of a
    site occupancy by the distance between the site and the solvent, and
    `residence_time` converts each of those correlation functions into a
    characteristic time — see
    [Residence time as a function of the distance](@ref occupancy-residence-profile).

```@autodocs
Modules = [ MolSimToolkit ]
Pages = [ 
    "intermittent_correlation.jl",
]
```


