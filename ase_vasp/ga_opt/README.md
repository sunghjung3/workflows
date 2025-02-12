## Optimization with ASE Genetic Algorithm

----------
### About

Taken from the [ASE example](https://wiki.fysik.dtu.dk/ase/tutorials/ga/ga_optimize.html), but modified to accomodate additional features:

- Supports SLURM (`slurmqueuerun.py`)
- Actually waits SLURM jobs to finish
- If either `main_run.py` or a SLURM job is prematurely terminated, running `main_run.py` will resume exactly where it left off. No need to manually reset or adjust code.


### Set up

In this example, we attempt to find the globally optimal positions of placing a given number of Mg atoms inside an initial skeleton structure. Both the Mg atoms and the positions of the skeleton are adjustable.

The directories and file are set up the following way (relative to `$root_dir`):

- `$root_dir`
    - This `README.md` is located here
    - Place the POSCAR of the empty skeleton structure here inside `POSCAR` (no Mg atoms yet)

- `$root_dir/Mg`
    - `initialize_db.py`
        - Run this script from `$root_dir/Mg`. It will create directories `Mg-$i` for `i` values specified at the bottom of this script, and each directory will contain `gadb.db` with various structures of skeleton containing `i` Mg atoms.
        - The chemical element, population size, and the valid placement volume can be easily modified.
    - `main_run.py`
        - Main script that should be run in the background of a cluster. Responsible for submiting jobs, making mutations, and waiting.
        - Run from `$root_dir/Mg/Mg-$i` as the working directory: `nohup python ../main_run.py &`
        - Adjust the GA settings inside this script
            - `initial_db_size` in this script should match `population_size` from the `initialize_db.py` script. Important for resuming termiated jobs.
        - Also don't forget to modify the `jtg()` function to customize the ASE VASP job submission script (SBATCH lines, environmental variables, etc.)
            - See `ase_vasp/single_pt_calc/slurm_job-arrays` of this repo for examples on setting the ASE VASP environmental variables.
    - `slurmqueuerun.py`
        - Modified `ase.ga.pbs_queue_run.PBSQueueRun` to work with SLURM.
        - No need to change anything in this file.

- `$root_dir/Mg/Mg-$i`
    - `calc.py`
        - Customize VASP settings here. Special care should be given for settings that are POSCAR dependent, like `magmom`.


### Output

- `$root_dir/Mg/Mg-$i/all_candidates.traj` will be written at the end, sorted by energy.
- All the intermediate ASE GA files will be located in `$root_dir/Mg/Mg-$i/tmp_ga`
- All the ASE VASP files will be located in `$root_dir/Mg/Mg-$i/tmp_ga/cand$j-vasp`
