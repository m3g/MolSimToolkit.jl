# Interface of the SLURM cluster manager
struct SLURM <: ClusterType end

# Add this cluster interface to the list of cluster types
push!(cluster_types, SLURM)

# Function that identifies the type of cluster manager
is_cluster_type(::Type{SLURM}) = haskey(ENV, "SLURM_JOB_ID")

# Get the job ID
get_jobid(::Type{SLURM}) = ENV["SLURM_JOB_ID"]

# Get the list of nodes reserved for the job: returns a list of strings with node names
function get_nodelist(::Type{SLURM})
    nodelist = filter(s -> length(s) > 0, split(
        read(pipeline(`scontrol show hostnames $(ENV["SLURM_JOB_NODELIST"])`), String),
        "\n"
    ))
    return nodelist
end

# Update the job with the list of nodes currently in use. Release the nodes that are not in use anymore.
function update_job!(::Type{SLURM}, nodelist::AbstractVector{String})
    hostlist = read(pipeline(`scontrol show hostlist "$(join(nodelist,","))"`), String)
    run(`scontrol update JobId=$(ENV["SLURM_JOB_ID"]) NodeList=$hostlist`)
    ENV["SLURM_JOB_NODELIST"] = hostlist
    return nothing
end

# Command to execute a code in a node in this cluster interface
run_task_on_node(::Type{SLURM}, task_run::Function, node::String, task) = run(`ssh $node $(task_run(task))`)


