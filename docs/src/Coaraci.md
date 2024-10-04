# Cluster submission management

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

# Define a function `task_finished` that recognizes if a simulation is finished:
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
        finish_file="finished.txt"
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
in the script above. At the end of all simulations, a `finished.txt` file is generated, which 
is then identified as the marker of a finished run.  

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

# Go into simulation dir
cd $simulation_dir

host=`hostname`
echo $host > running_in_hostname.txt # Just write to a file the current hostname, for checking

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

# Run analyzes
julia -t $n mddf_water.jl
julia -t $n mddf_${cossolvent}.jl

\rm -f running_in_hostname.txt
# Create the file that indicates that the simulation is finished
echo "finished simulation in node $host" > finished.txt
```
