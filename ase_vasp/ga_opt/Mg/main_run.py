from random import random
import time, os, socket, sys

from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from slurmqueuerun import SLURMQueueRun
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.standardmutations import (
    MirrorMutation,
    PermutationMutation,
    RattleMutation,
)
from ase.ga.utilities import closest_distances_generator, get_all_atom_types
from ase.io import write


def jtg(job_name, traj_file):
    s = '#!/bin/bash\n'
    s += '\n'
    s += f'#SBATCH --job-name {job_name}\n'
    s += '#SBATCH --nodes=1\n'  # per structure
    s += '#SBATCH --ntasks-per-node=128\n'
    s += '#SBATCH --partition=wholenode\n'
    s += '#SBATCH --time=96:00:00\n'
    s += f'#SBATCH --output={tmp_folder}/{job_name}_%j.log\n'
    #s += f'#SBATCH --output=/dev/null\n'
    s += '\n'
    s += 'date\n'
    s += 'ml reset\n'
    s += 'ml load intel-mkl\n'
    s += 'export PATH=/home/x-graeme/vasp/vasp.6.4.3/bin:$PATH\n'
    s += 'export VASP_PP_PATH="/anvil/scratch/x-sjung3/prussian_blue/vasp_potpaws"\n'
    s += 'export ASE_VASP_COMMAND="mpirun -np $SLURM_NTASKS vasp_std"\n'
    s += f'python calc.py {traj_file}\n'
    s += 'date\n'
    return s


population_size = 100
mutation_probability = 0.3
n_to_test = 120
slurm_sleep_interval = 300  # seconds
slurm_job_prefix = "GA_" + os.path.basename(os.getcwd())
max_n_jobs_relax = 207  # be mindful of SLURM restrictions on number of running/pending jobs
max_n_jobs_ga = 5  # be mindful of SLURM restrictions on number of running/pending jobs

initial_db_size = 20  # should match population_size from initialize_db.py

print(f"Hostname: {socket.gethostname()}", flush=True)
print(f"Process ID: {os.getpid()}", flush=True)

# Initialize the different components of the GA
da = DataConnection('gadb.db')
tmp_folder = 'tmp_ga/'
slurm_run = SLURMQueueRun(da,
                          tmp_folder=tmp_folder,
                          job_prefix=slurm_job_prefix,
                          n_relax=max_n_jobs_relax,
                          n_ga=max_n_jobs_ga,
                          job_template_generator=jtg,
                          )

atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_to_optimize = len(atom_numbers_to_optimize)
skeleton = da.get_slab()
all_atom_types = get_all_atom_types(skeleton, atom_numbers_to_optimize)
blmin = closest_distances_generator(all_atom_types,
                                    ratio_of_covalent_radii=0.7)

comp = InteratomicDistanceComparator(n_top=n_to_optimize,
                                     pair_cor_cum_diff=0.015,
                                     pair_cor_max=0.7,
                                     dE=0.02,
                                     mic=False)

pairing = CutAndSplicePairing(skeleton, n_to_optimize, blmin, number_of_variable_cell_vectors=3)
mutations = OperationSelector([1., 1., 1.],
                              [MirrorMutation(blmin, n_to_optimize),
                               RattleMutation(blmin, n_to_optimize),
                               PermutationMutation(n_to_optimize)])

# Relax all unrelaxed structures (e.g. the starting population)
slurm_run.__cleanup__()
print("# unrelaxed:", da.get_number_of_unrelaxed_candidates())
print("# relaxed:", len(da.get_all_relaxed_candidates()))
print("# previously queued:", len(da.get_all_candidates_in_queue()))
sys.stdout.flush()
while True:  # outer loop for case where SLURM job gets terminated while this script is running
    # the only way to exit this loop is if there are no new candidates to submit AND no jobs are running
    for qid in da.get_all_candidates_in_queue():  # reset terminated structures so that calc.py can resume them
        if not slurm_run.is_running(qid):
            da.remove_from_queue(qid)

    need_to_run_more = False
    while da.get_number_of_unrelaxed_candidates() > 0:
        need_to_run_more = True
        print(f'{da.get_number_of_unrelaxed_candidates()} more to relax', flush=True)
        print(f'{slurm_run.number_of_jobs_running()} jobs running', flush=True)
        while (da.get_number_of_unrelaxed_candidates() > 0) and (not slurm_run.enough_jobs_running_relax()):
            a = da.get_an_unrelaxed_candidate()
            slurm_exit_code = slurm_run.relax(a)
            if slurm_exit_code:
                raise ValueError("Failed to submit relaxation job\n"
                                 "I refuse to continue\n"
                                 "Cancel all jobs, resolve issue, and try again\n"
                                 )
        time.sleep(slurm_sleep_interval)

    while slurm_run.number_of_jobs_running() > 0:
        need_to_run_more = True
        print(f'Waiting for {slurm_run.number_of_jobs_running()} jobs to finish', flush=True)
        time.sleep(slurm_sleep_interval)

    if not need_to_run_more:
        break
slurm_run.__cleanup__()

# create the population
population = Population(data_connection=da,
                        population_size=population_size,
                        comparator=comp)

# Submit new candidates until enough are running
n_tested = len(da.get_all_relaxed_candidates()) - initial_db_size
while n_tested < n_to_test:
    while (not slurm_run.enough_jobs_running_ga() and
        len(population.get_current_population()) >= 2 and
        n_tested < n_to_test):
        a1, a2 = population.get_two_candidates()
        a3, desc = pairing.get_new_individual([a1, a2])
        if a3 is None:
            continue
        da.add_unrelaxed_candidate(a3, description=desc)

        if random() < mutation_probability:
            a3_mut, desc = mutations.get_new_individual([a3])
            if a3_mut is not None:
                da.add_unrelaxed_step(a3_mut, desc)
                a3 = a3_mut
        slurm_exit_code = slurm_run.relax(a3)
        if slurm_exit_code:
            raise ValueError("Failed to submit relaxation job\n"
                             "I refuse to continue\n"
                             "Cancel all jobs, resolve issue, and try again\n"
                             )

        population.update()
        n_tested += 1
    print(f'{n_tested} candidates tested')
    print(f'{slurm_run.number_of_jobs_running()} jobs running')
    print(f'Current population: {len(population.get_current_population())}')
    sys.stdout.flush()
    time.sleep(slurm_sleep_interval)

# Wait for all jobs to finish
while slurm_run.number_of_jobs_running() > 0:
    print(f'Waiting for {slurm_run.number_of_jobs_running()} jobs to finish', flush=True)
    time.sleep(slurm_sleep_interval)

write('all_candidates.traj', da.get_all_relaxed_candidates())
