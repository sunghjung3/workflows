## Parallelized VASP single-point calculations using ASE's API and SLURM job arrays

----------

### Set up

- POTCAR
    - Either:
        - Place the POTCAR file here (e.g. `potpaw_PBE/<element>/POTCAR` for PBE XC), or
        - Edit the `VASP_PP_PATH` variable in `run.sh`
            - NOTE: `..` points to this root directory for this variable

- `run.sh`
        
    - Number of calculations to run concurrently: edit the last integer in the line
        
            #SBATCH --array=0-19  # 20 concurrent calculations

    - Email to send job status notifications (edit or delete line)

            #SBATCH --mail-user=<email>

    - Parallelization over nodes and MPI

            #SBATCH --nodes=1  # Run on 1 node
            #SBATCH --ntasks=48  # Run on 48 MPI ranks
    
    - SLURM partition on current cluster

            #SBATCH --partition=<partition>  # SLURM partition

    - Location of the trajectory to calculate: `traj_path`
        - Recommended to put the **absolute path**

    - File name prefix for output trajectory files: `out_traj_prefix`
        - Trajectory of converged calculations: `<out_traj_prefix>_ef.traj`
        - Trajectory of unconverged calculations: `<out_traj_prefix>_unconv.traj`


- `ase_vasp_run.py`

    - INCAR settings: modify parameters to `Vasp(...)` function


### Output

- Each task in the job array will create its own directory after its array index.
    Each of these directories will contain its own output trajectory files.
- After the conclusion of all tasks in the array, task 0 will collect all trajectory files
    into a single converged file and unconverged file in the root directory.
- The job array directories and the log files can be kept or deleted.
