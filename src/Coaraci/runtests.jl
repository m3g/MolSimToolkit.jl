@testsnippet CoaraciTests begin
    using MolSimToolkit: Coaraci
    # Print to stdout
    Coaraci.options.stdout = true
    if haskey(ENV, "CLUSTER_INTERFACE_TEST")
        pop!(ENV, "CLUSTER_INTERFACE_TEST")
    end
    @kwdef struct TestTask 
        title::String
        delay::Int
        output_file::String
    end
    task_finished(t::TestTask) = begin
        if t.title == "ERR"
            return false
        end
        isfile(t.output_file)
    end
    task_title(t::TestTask) = t.title
    function task_run(t::TestTask) 
        sleep(t.delay)
        open(t.output_file, "w") do io
            println(io, "Task $(t.title) finished.")
        end
        return nothing
    end
end

@testitem "Argument errors" setup=[CoaraciTests] begin
    task_list = [ TestTask("A", 1, tempname()), TestTask("B", 2, tempname()) ]
    @test_throws ArgumentError Coaraci.simulate(task_list; task_finished, task_title, task_run)
    ENV["CLUSTER_INTERFACE_TEST"] = "true"
    @test_throws ArgumentError Coaraci.simulate(task_list; task_finished, task_title)
    @test_throws ArgumentError Coaraci.simulate(task_list; task_run, task_title)
    @test_throws ArgumentError Coaraci.simulate(task_list; task_run, task_finished)
end

@testitem "Cluster interface test" setup=[CoaraciTests] begin
    ENV["CLUSTER_INTERFACE_TEST"] = "true"
    task_list = [ TestTask("A", 1, tempname()), TestTask("B", 2, tempname()) ]
    rm.(getfield.(task_list,:output_file); force=true)
    Coaraci.simulate(task_list; task_run, task_finished, task_title)
    @test all(task_finished, task_list)
    @test occursin("finished all tasks:", read("coaraci.log", String))
    @test occursin("finished with errors: 0", read("coaraci.log", String))

    # One task fails:
    rm.(getfield.(task_list,:output_file); force=true)
    push!(task_list, TestTask("ERR", 3, tempname()))
    Coaraci.simulate(task_list; task_run, task_finished, task_title)
    @test sum(task_finished, task_list) == 2
    @test occursin("finished with errors: 1", read("coaraci.log", String))

    # Two tasks per node
    task_list = [ TestTask("A", 1, tempname()), TestTask("B", 2, tempname()) ]
    rm.(getfield.(task_list,:output_file); force=true)
    Coaraci.simulate(task_list; task_run, task_finished, task_title, ntasks_per_node=2)
    @test sum(task_finished, task_list) == 2
    @test occursin("2 in this node", read("coaraci.log", String))
    @test occursin("Keeping nodes: node1", read("coaraci.log", String))

    # Three tasks per node, 4 tasks
    task_list = [ 
        TestTask("A", 2, tempname()), 
        TestTask("B", 2, tempname()),
        TestTask("C", 2, tempname()), 
        TestTask("D", 2, tempname()), 
    ]
    rm.(getfield.(task_list,:output_file); force=true)
    Coaraci.simulate(task_list; task_run, task_finished, task_title, ntasks_per_node=3)
    @test sum(task_finished, task_list) == 4
    @test occursin("3 in this node", read("coaraci.log", String))
    @test occursin("D in node node2 finished", read("coaraci.log", String))
    # Test don´t run tasks finished
    Coaraci.simulate(task_list; task_run, task_finished, task_title, ntasks_per_node=3)
    @test occursin("A is finished: won't run", read("coaraci.log", String))

    # Test overloading the internal functions
    Coaraci.task_finished(t::TestTask) = isfile(t.output_file)
    Coaraci.task_title(t::TestTask) = t.title
    function Coaraci.task_run(t::TestTask) 
        sleep(t.delay)
        open(t.output_file, "w") do io
            println(io, "Task $(t.title) finished.")
        end
        return nothing
    end
    rm.(getfield.(task_list,:output_file); force=true)
    Coaraci.simulate(task_list; ntasks_per_node=3)
    @test sum(task_finished, task_list) == 4
end
