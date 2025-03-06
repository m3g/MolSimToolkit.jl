# [Cluster submission management](@id Coaraci)

!!! warning
    This is an experimental feature. Breaking changes may occur without 
    a breaking package release.

The `Coaraci` module of `MolSimToolkit` manages the submission of multiple independent tasks in a single job of a cluster, where the job requested multiple nodes.
For example, a job requests 20 nodes from a SLURM computer cluster, and the goal is to run independent simulations on each of the nodes.
This script will distribute the simulations in the nodes and, when more nodes than simulations are available, release the spare nodes
to the cluster.

Currently fully supports the SLURM cluster manager. Partial support for PBS is available.

!!! compat
    The `Coaraci` module was introduced in MolSimToolkit v1.20.0

## SLURM/Gromacs example

### The submission script

Consider the following script to submit a job in a SLURM cluster. Here, we aim to run many (more than 20)
independent simulations, and we request 20 nodes from the cluster. Each simulation will run in 48 cores
inside each node.

```bash
#!/bin/bash
#SBATCH -n 48
#SBATCH --job-name test
#SBATCH --nodes 20
#SBATCH --partition=paralela
#SBATCH -o /home/users/test/proteina/job_files/out_tmao_replica13.txt
#SBATCH -e /home/users/test/proteina/job_files/err_tmao_replica13.txt

cd $SLURM_SUBMIT_DIR

#### Run the Julia script
julia simulations.jl
```

### The Julia script

Let us call this script `simulations.jl`, which is invoked by the submission script, as shown above.

```julia
using MolSimToolkit: Coaraci

# Define a structure that contains all the data required to run your simulations.
# none of the fields is mandatory, but here they will be used.
@kwdef struct SimulationData
    title::String
    dir::String
    cosolvent::String
    finish_file::String
end

#
# MANDATORY: Define these 3 functions:
#
# Define `task_run`: which is the command that runs a single simulation, possibly using as arguments
# the simulation data
Coaraci.task_run(s::SimulationData) = "/home/ander/run_script/run_single_simulation.sh $(s.dir) $(s.cosolvent)"

# Define a function `task_finished` that recognizes if a simulation was already run.
Coaraci.task_finished(s::SimulationData) = isfile(joinpath(s.dir,s.finish_file))

# Define a function `title` that sets the title of the run:
Coaraci.task_title(s::SimulationData) = s.title

# Simulation list: the following will create a list of simulations
# that have to be executed. Each simulation is defined by a directory.
# Here, we construct the directory names with the specifiers of each
# simulation. Important: the result is a list of **full** directories
# paths for each simulation.
base_dir="/home/users/test/protein"
concentrations = [ "0.1" "0.2" ]
replicas = [ "replica1", "replica2" ]
# This is the crucial part: build a proper list of simulation objects:
simulation_list = SimulationData[]
for conc in concentrations, rep in replicas
    sim = SimulationData(
        title="$conc/$rep",
        dir=joinpath(base_dir,conc,rep),
        cosolvent="tmao",
        finish_file="production.gro",
    )
    push!(simulation_list, sim)
end

# Execute the script that runs the simulations in parallel, one per node.
# The name of the log file of the script is set here.
Coaraci.simulate(
    simulation_list; 
    logfile="mysimulations.log",
    #ntasks_per_node=1
)
```

The optional `ntasks_per_node` parameter can be set to define how many tasks are run per node. 
In this case, for example, if `ntasks_per_node=2`, the script that runs each task inside the 
node should request only half of the CPUs of the node, for maximal efficiency. By default,
`ntasks_per_node=1`.  

### The script that runs a simulation inside a node

The following script receives 2 arguments, which are provided by the definition of `task_run`,
in the script above. At the end of all simulations, the `production.gro` file is generated, which 
is then identified as the marker of a finished run. A finished run will not be run again.

And, finally, the script that submit the actual simulations, **from within the node**, and
using, here 48 cores (in a single node), is, for example:

```bash
#!/bin/bash
simulation_dir=$1 # first argument of script

# second argument of script: in this case, the cosolvent: can be used
# to label or run specific analyzes for the same simulation
cosolvent=$2

n=48 # processors in this node

#### modules necessary to run Gromacs #########
module load hwloc
module load gnu12
module load openmpi4
module load gcc-runtime/12.2.0-gcc-12.2.0-lbbfl34
module load gcc-runtime/12.2.0-gcc-12.2.0-sqqkkcb
module load gromacs/2023-gcc-12.2.0-zktp5lx

# Important: go into simulation dir
cd $simulation_dir

# Run simulations: gromacs here
/home/users/apereira/softwares/packmol-20.15.1/packmol < box.inp > box.log
gmx_mpi grompp -f mim.mdp -c solvated.pdb -r solvated.pdb -p processed.top -o em.tpr -maxwarn 3 > /dev/null
mpirun -np $n gmx_mpi mdrun -v -deffnm em -ntomp 1  > /dev/null #  -npme 12 -dd 4 3 3

gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p processed.top -o nvt.tpr -maxwarn 3 > /dev/null
mpirun -np $n gmx_mpi mdrun -v -deffnm nvt -ntomp 1 > /dev/null

gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -p processed.top -o npt.tpr -maxwarn 3 > /dev/null
mpirun -np $n gmx_mpi mdrun -v -deffnm npt -ntomp 1 > /dev/null

gmx_mpi grompp -f prod.mdp -c npt.gro -r npt.gro -p processed.top -o prod.tpr -maxwarn 3 > /dev/null
mpirun -np $n gmx_mpi mdrun -v -deffnm prod -ntomp 1 > /dev/null
```

## PBS example

The following script submits a job the `parexp` queue of a PBS cluster:

```bash
#PBS -N meuteste
#PBS -q parexp
#PBS -l nodes=2:ppn=48
#PBS -m abe
#PBS -e erros
#PBS -o saida

cd $PBS_O_WORKDIR
julia submission.jl
```

Where `submission.jl` is the following script, which sets up two very simple tasks that just wait some `delay`
seconds to finish:

```julia
using MolSimToolkit: Coaraci
@kwdef struct TestTask 
    title::String
    delay::Int
    output_file::String
end
dir="/home/lovelace/proj/proj864/user/test/"
Coaraci.task_finished(t::TestTask) = isfile(t.output_file)
Coaraci.task_title(t::TestTask) = t.title
Coaraci.task_run(t::TestTask) = "$dir/run.sh $dir $(t.title) $(t.delay) $(t.output_file)"
task_list = [ 
    TestTask(title="A", delay=30, output_file="end1.txt"), 
    TestTask(title="B", delay= 2, output_file="end2.txt"),
]
Coaraci.simulate(task_list; ntasks_per_node=1)
```

and the `run.sh` script being executed is:

```bash
#!/bin/bash
cd $1 # running dir
sleep $3 # delay
hname=`hostname`
echo "running! $2 in $hname" > $4 # title and output_file name
```

Note that the delay is passed as the third argument of the execution. The script will run a 30s 
long task on the first node, and a 2s task on the second node. When the second task finishes, the
second node is released. Note: PBS does not allow releasing the primary node, thus it will be kept
even if no task is running on it.

## Example log file

A typical log file will look like:

```bash
=============================================================
Starting Coaraci managed submissions: 2024-10-11 10:48:54
=============================================================
Cluster type: MolSimToolkit.Coaraci.PBS
JOB ID: 144830.ada
Nodes: adano62, adano72
Number of tasks per node: 1
-------------------------------------------------------------
Number of tasks to run: 2
-------------------------------------------------------------
1: running A in node adano62 (1 in this node)
2: running B in node adano72 (1 in this node)
2: B in node adano72 finished successfully after 2 seconds, 507 milliseconds.
> 1 nodes released. Keeping nodes: adano62
1: A in node adano62 finished successfully after 30 seconds, 480 milliseconds.
-------------------------------------------------------------
finished all tasks:
 - tasks finished successfully: 2
 - tasks finished with errors: 0
=============================================================
End Coaraci managed submissions: 2024-10-11 10:49:24
Total running time: 30 seconds, 762 milliseconds
=============================================================
```
