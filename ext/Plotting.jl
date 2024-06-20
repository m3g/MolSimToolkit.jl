module Plotting
    using TestItems: @testitem
    using MolSimToolkit
    import Plots
    include("./MolSimStyle.jl")
    include("./BlockAverages.jl")
    include("./REMD.jl")
    include("./SecondaryStructure.jl")
end