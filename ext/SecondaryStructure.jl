function MolSimToolkit.ss_heatmap(
    ssmap::AbstractMatrix{<:Integer}; 
    scalex=1.0, 
    kargs...
)
    default = Dict(
        :framestyle => :box,
        :color => Plots.palette(:tab20c,10),
        :clims => (0.5,10.5),
        :xlabel => "frame",
        :ylabel => "residue",
    )
    plt = Plots.heatmap(MolSimStyle, scalex*(1:size(ssmap,2)), 1:size(ssmap,1), ssmap; _kargs(default; kargs)...)
    return plt
end

