#!/bin/bash
#!
#! Example SLURM job script for Wilkes (K20c Tesla nodes)
#! Last updated: Wed Jan  8 20:05:59 GMT 2014
#!


help()
{
	echo "This program takes output smaug data and gathers into a single outputfile"
        echo "./smauggather.sh numproc start finish step"

	echo "numproc is the number of processors"
        echo "start is the starting step"
        echo "finish is the final step"
        echo "step is optional and is the interval"

	exit 1
}




#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J sgather_ngdr
#! Which project should be charged (NB Wilkes projects end in '-GPU'):
#SBATCH -A SHEFFIELD-GPU
#! How many whole nodes should be allocated?
#SBATCH --nodes=2
#! How many (MPI) tasks in total should be launched?
#SBATCH --ntasks=4
#! How much wallclock time will be required?
#SBATCH --time=02:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --requeue

#! Do not change:
#SBATCH -p tesla

#! sbatch directives end here (put any additional directives above this line)
#!
#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
ppn=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module load default-wilkes                 # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:
#module load cuda/5.5
#module load mvapich2/gcc/2.0b/gdr
module load cuda/6.5
module load mvapich2-GDR/gnu/2.0-8_cuda-6.5

export MV2_USE_CUDA=1
export MV2_USE_GPUDIRECT=1
export MV2_RAIL_SHARING_POLICY=FIXED_MAPPING
export MV2_PROCESS_TO_RAIL_MAPPING=mlx5_0:mlx5_1
export MV2_RAIL_SHARING_LARGE_MSG_THRESHOLD=1G
export MV2_CPU_BINDING_LEVEL=SOCKET
export MV2_CPU_BINDING_POLICY=SCATTER
export CUDA_VISIBLE_DEVICE=0



#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

#! Override the default value of ppn (= Processes Per Node) if necessary
#! (e.g. because the memory required per process will otherwise overload the node).
#! Note: since each tesla node has 2 GPUs, your code may require you to set ppn=2 here.
#! Do this by uncommenting and editing the following line (ensure ppn<=12):
ppn=2

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 12:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
mpi_tasks_per_node=$[$ppn/$OMP_NUM_THREADS]
np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.


#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
#! Delete gpu_run_task.sh if you do _not_ wish to bind even MPI ranks to GPU0 and IB card0,
#! and odd ranks to GPU1 and IB card1. This binding may improve performance and is
#! REQUIRED for Turbostream.
#CMD="mpirun -ppn $mpi_tasks_per_node -np $np gpu_run_task.sh $application $options"
CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"


###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, ppn=$ppn (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"




start=0
end=50 
step=1
current=$start



numproc=4



#! Run options for the application:
#options="-d cuda D D"
#options="a"
options=" -genvall "


while [ "$current" -le $end ]
do

        #eval $CMD
        #mpirun -np $numproc bin/smaug -o gather $current
        #mpirun -npernode $mpi_tasks_per_node -np $np $application -o gather $current

	#! Full path to application executable: 
	#application="osu_latency"
	application="bin/smaug -o gather $current"

	mpirun -ppn $mpi_tasks_per_node -np $np $options $application
        current=`expr $current + $step`

done


