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

# Node that runs the job script, which cannot be released
batch_host(::Type{CLUSTER_INTERFACE_TEST}) = get(ENV, "CLUSTER_INTERFACE_TEST_BATCH_HOST", "node1")

# Update the job with the list of nodes currently in use. Release the nodes that are not in use anymore.
# `nodes_to_keep` are the names of the nodes in use, and `nodes_to_release` the names of the idle nodes.
# Here, an error is thrown if the batch host would be released (which terminates the job in a real
# cluster), and the refusal of the cluster to release nodes can be simulated by setting
# the CLUSTER_INTERFACE_TEST_FAIL_UPDATE environment variable.
function update_job!(::Type{CLUSTER_INTERFACE_TEST}, nodes_to_keep::AbstractVector{String}, nodes_to_release::AbstractVector{String})
    if batch_host(CLUSTER_INTERFACE_TEST) in nodes_to_release || !(batch_host(CLUSTER_INTERFACE_TEST) in nodes_to_keep)
        error("Batch host released: job terminated.")
    end
    haskey(ENV, "CLUSTER_INTERFACE_TEST_FAIL_UPDATE") && error("Cluster refused to release nodes.")
    return nothing
end

# Command to execute a code in a node in this cluster interface
run_task_on_node(::Type{CLUSTER_INTERFACE_TEST}, task_run::Function, node::String, task) = task_run(task)
