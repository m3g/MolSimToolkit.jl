# Interface of the PBS cluster manager
struct PBS <: ClusterType end

# Get the job ID
get_jobid(::Type{PBS}) = ENV["PBS_JOBID"]

# Add this cluster interface to the list of cluster types
push!(cluster_types, PBS)

# Function that identifies the type of cluster manager
is_cluster_type(::Type{PBS}) = haskey(ENV, "PBS_JOBID")

# Get the list of nodes reserved for the job: returns a list of strings with node names
function get_nodelist(::Type{PBS})
    nodelist = filter(s -> length(s) > 0, split(
        read(pipeline(`uniq /var/spool/pbs/aux/$(ENV["PBS_JOBID"])`), String),
        "\n"
    ))
    return nodelist
end

# Update the job with the list of nodes currently in use. Release the nodes that are not in use anymore.
function update_job!(::Type{PBS}, nodelist)
    return nothing
    if length(nodes_currently_in_use) < length(nodelist)
        resize!(nodelist, length(nodes_currently_in_use)) 
        nodelist .= nodes_currently_in_use
        hostlist = read(pipeline(`scontrol show hostlist "$(join(nodelist,","))"`), String)
        run(`scontrol update JobId=$(ENV["SLURM_JOB_ID"]) NodeList=$hostlist`)
        ENV["SLURM_JOB_NODELIST"] = hostlist
        print_flush(log, "Nodes released. Keeping nodes: $(ENV["SLURM_JOB_NODELIST"])")
    end
    return nothing
end

# Command to execute a code in a node in this cluster interface
run_task_on_node(::Type{PBS}, task_run::Function, node::String, task) = run(`ssh $node $(task_run(task))`)

