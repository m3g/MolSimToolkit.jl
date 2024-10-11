# Interface of the SLURM cluster manager
struct CLUSTER_INTERFACE_TEST <: ClusterType end

# Add this cluster interface to the list of cluster types
push!(cluster_types, CLUSTER_INTERFACE_TEST)

# Function that identifies that this is an interface test
is_cluster_type(::Type{CLUSTER_INTERFACE_TEST}) = haskey(ENV, "CLUSTER_INTERFACE_TEST")

# Get the job ID
get_jobid(::Type{CLUSTER_INTERFACE_TEST}) = "1.test"

# Get the list of nodes reserved for the job: returns a list of strings with node names
function get_nodelist(::Type{CLUSTER_INTERFACE_TEST})
    nodelist = [ "node$i" for i in 1:5 ]
    return nodelist
end

# Update the job with the list of nodes currently in use. Release the nodes that are not in use anymore.
# The nodelist is a list of strings with node names of nodes that ir in use.
function update_job!(::Type{CLUSTER_INTERFACE_TEST}, nodelist::AbstractVector{String})
    return nothing
end

# Command to execute a code in a node in this cluster interface
run_task_on_node(::Type{CLUSTER_INTERFACE_TEST}, task_run::Function, node::String, task) = task_run(task)
