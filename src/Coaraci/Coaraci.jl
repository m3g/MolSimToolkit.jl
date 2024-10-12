module Coaraci

import ..MolSimToolkit: version
import Dates: now, format, canonicalize
using Base.Threads: @spawn, nthreads

# Functions that must be defined by the user for the custom simulation type
function task_run end
function task_finished end
function task_title end

# List of supported cluster types
abstract type ClusterType end
const cluster_types = Type{<:ClusterType}[]
#
# Interfaces for different cluster managers
#
include("./SLURM.jl")
include("./PBS.jl")
include("./CLUSTER_INTERFACE_TEST.jl")

# node data
mutable struct Node
    const name::String
    running_tasks::Int
end
name(node::Node) = node.name
import Base: ==
==(a::Node, b::Node) = a.name == b.name

# Internal Options
@kwdef mutable struct Options
    # print, or not, to stdout
    stdout::Bool = false
    # sleep before checking again for finished tasks (s)
    sleep::Float64 = 0.5
end
const options = Options()

# Type to name tasks available
struct AvailableTask end

# Custom printing that does not buffer
function print_flush(io, val) 
    println(io, val); flush(io)
    options.stdout && println(val)
    return nothing
end

# Function that checks input arguments and running options
function check_input(task_list, nodelist; task_run, task_finished, task_title)
    if !hasmethod(task_run, (eltype(task_list),))
        throw(ArgumentError("Coaraci.task_run function not defined for type $(eltype(task_list))"))
    end
    if !hasmethod(task_finished, (eltype(task_list),))
        throw(ArgumentError("Coaraci.task_finished function not defined for type $(eltype(task_list))"))
    end
    if !hasmethod(task_title, (eltype(task_list),))
        throw(ArgumentError("Coaraci.task_title function not defined for type $(eltype(task_list))"))
    end
    if length(nodelist) == 0
        throw(ArgumentError("No nodes available for running simulations."))
    end
    return nothing
end

function find_cluster_type()
    i_cluster_type = findfirst(is_cluster_type, cluster_types)
    isnothing(i_cluster_type) && throw(ArgumentError("Cluster type not detected."))
    return cluster_types[i_cluster_type]
end

function simulate(
    task_list::AbstractVector; 
    logfile="coaraci.log",
    ntasks_per_node=1,
    cluster_type=find_cluster_type(),
    task_run=task_run,
    task_finished=task_finished,
    task_title=task_title,
)

    # Remove log file, if present:
    rm(logfile; force=true)

    # The following lines obtain a list of the
    # nodes that were reserved for this job. The list is put into channels, which
    # will be taken and freed as the task start and finish. 
    nodelist = Node[]
    for name in get_nodelist(cluster_type)
        push!(nodelist, Node(name, 0))
    end

    # Check if all the input arguments are correct
    check_input(task_list, nodelist; task_run, task_finished, task_title)

    # Lock to control access to the list of nodes in use
    lock_nodes = ReentrantLock()

    # List of (Julia) tasks
    julia_tasks = Task[]

    # Channels with available execution tasks (nodes * ntasks_per_node 
    available_tasks = Channel{AvailableTask}(length(nodelist) * ntasks_per_node)
    for _ in nodelist
        for _ in 1:ntasks_per_node
            put!(available_tasks, AvailableTask())
        end
    end

    # LOOP that runs the tasks
    open(logfile, "w") do log # open log file to write
        init_time = now()
        print_flush(log, "=============================================================")
        print_flush(log, "MolSimToolkit version: $version")
        print_flush(log, "Starting Coaraci managed submissions: $(format(init_time, "yyyy-mm-dd HH:MM:SS"))")
        print_flush(log, "=============================================================")
        print_flush(log, "Cluster type: $cluster_type")
        print_flush(log, "JOB ID: $(get_jobid(cluster_type))")
        print_flush(log, "Nodes: $(join(name.(nodelist), ", "))") # print nodelist to log
        print_flush(log, "Number of tasks per node: $(ntasks_per_node)")
        print_flush(log, "-------------------------------------------------------------")

        # Remove from task list the tasks that are already finished
        tasks_to_run = eltype(task_list)[]
        for task in task_list
            if task_finished(task)
                print_flush(log, "$(task_title(task)) is finished: won't run.")
            else
                push!(tasks_to_run, task)
            end
        end
        print_flush(log, "Number of tasks to run: $(length(tasks_to_run))")
        print_flush(log, "-------------------------------------------------------------")

        # Create a list that will contain the list of nodes currently in use,
        # to retain them if other nodes are not required anymore
        nodes_currently_in_use = Node[] 
        for (irun, task) in enumerate(tasks_to_run)
            local node
            local inode = nothing
            # Take an available channel task (will be blocked if no task is available)
            used_task = take!(available_tasks)
            # Update list of nodes currently in use
            @lock lock_nodes begin
                inode = findfirst(n -> n.running_tasks < ntasks_per_node, nodelist)
                if isnothing(inode)
                    error("Could not find a node with available tasks.")
                end
                node = nodelist[inode]
                if node.running_tasks == 0
                    push!(nodes_currently_in_use, node)
                end
                node.running_tasks += 1
            end
            print_flush(log, "$irun: running $(task_title(task)) in node $(name(node)) ($(node.running_tasks) in this node)")
            # Launch process that executes the task in the available node
            t = @spawn begin
                # Executes the task in the node.
                run_task_on_node(cluster_type, task_run, name(node), task)
                # When the task finishes, report to the log file and return the node to the channel list.
                if task_finished(task)
                    print_flush(log, "$irun: $(task_title(task)) in node $(name(node)) finished successfully after $(canonicalize(now() - init_time)).")
                else
                    print_flush(log, "$irun: $(task_title(task)) in node $(name(node)) finished with errors after $(canonicalize(now() - init_time)).")
                end
                # Remove node from list of currently used nodes
                @lock lock_nodes begin
                    node.running_tasks -= 1
                    if node.running_tasks == 0
                        filter!(!=(node), nodes_currently_in_use) 
                    end
                end
                # Release available task
                put!(available_tasks, used_task)
            end
            push!(julia_tasks, t)
        end
        # This loop runs until all tasks finish. If there are free nodes, reconfigure
        # slurm to release the free nodes and keep only the nodes currently in use.
        while any(!istaskdone, julia_tasks)
            @lock lock_nodes begin 
                if length(nodes_currently_in_use) < length(nodelist) 
                    n_released_nodes = length(nodelist) - length(nodes_currently_in_use)
                    resize!(nodelist, length(nodes_currently_in_use)) 
                    nodelist .= nodes_currently_in_use
                    nodelist_names = name.(nodelist)
                    update_job!(cluster_type, name.(nodelist))
                    print_flush(log, "> $(n_released_nodes) nodes released. Keeping nodes: $(join(nodelist_names, ", "))")
                end
            end
            sleep(options.sleep)
        end

        # Report final
        print_flush(log, "-------------------------------------------------------------")
        print_flush(log, "finished all tasks:")
        print_flush(log, " - tasks finished successfully: $(sum(task_finished, task_list))")
        print_flush(log, " - tasks finished with errors: $(sum(!task_finished, task_list))")
        print_flush(log, "=============================================================")
        print_flush(log, "End Coaraci managed submissions: $(format(now(), "yyyy-mm-dd HH:MM:SS"))")
        print_flush(log, "Total running time: $(canonicalize(now() - init_time))")
        print_flush(log, "=============================================================")
    end
end

end
