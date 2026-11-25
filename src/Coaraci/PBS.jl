# Interface of the PBS cluster manager
struct PBS <: ClusterType end

# Get the job ID
get_jobid(::Type{PBS}) = ENV["PBS_JOBID"]

# Add this cluster interface to the list of cluster types
push!(cluster_types, PBS)

# Function that identifies the type of cluster manager
is_cluster_type(::Type{PBS}) = haskey(ENV, "PBS_JOBID")

# Get the list of nodes reserved for the job: returns a list of strings with node names
# (the node file has one line per reserved cpu, thus the node names are repeated)
get_nodelist(::Type{PBS}) = unique!(String.(filter(!isempty, strip.(readlines(ENV["PBS_NODEFILE"])))))

# Node that runs the job script: the primary execution host, which is the first node of the node file.
batch_host(::Type{PBS}) = first(get_nodelist(PBS))

# Release the nodes that are not in use anymore.
function update_job!(::Type{PBS}, ::AbstractVector{String}, nodes_to_release::AbstractVector{String})
    # PBS does not allow releasing the primary node
    nodes_to_release = filter(!=(batch_host(PBS)), nodes_to_release)
    if length(nodes_to_release) > 0
        run(`pbs_release_nodes -j $(get_jobid(PBS)) $nodes_to_release`)
    end
    return nothing
end

# Command to execute a code in a node in this cluster interface
run_task_on_node(::Type{PBS}, task_run::Function, node::String, task) = run(`ssh $node $(task_run(task))`)

