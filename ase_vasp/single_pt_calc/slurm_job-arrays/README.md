## Parallelized VASP single-point calculations using ASE's API and SLURM job arrays

----------

### Set up

- POTCAR
    - Either:
        - Place the POTCAR file here (e.g. `potpaw_PBE/<element>/POTCAR` for PBE XC), or
        - Edit the `VASP_PP_PATH` variable in `run.sh`
            - NOTE: `..` points to this root directory for this variable

- `run.sh`
    - Location of the trajectory to calculate: `traj_path`
        - Recommended to put the **absolute path**
    - Number of calculations to run concurrently: edit the last integer in the line
        
            #SBATCH --array=0-19  # 20 concurrent calculations

    - Email to send job status notifications (edit or delete line)

            #SBATCH --mail-user=<email>

    - Parallelization over nodes and MPI

            #SBATCH --nodes=1  # Run on 1 node
            #SBATCH --ntasks=48  # Run on 48 MPI ranks
    
    - SLURM partition on current cluster

            #SBATCH --partition=<partition>  # SLURM partition

- `ase_vasp_run.py`

    - Name of trajectory file to write converged structures with energy and forces: `B_data_ef`
    - Name of trajectory file to write unconverged structures: `B_data_unconv`
    - INCAR settings: modify parameters to `Vasp(...)` function


### Output

- Each task in the job array will create its own directory after its array index.
    Each of these directories will contain its own output trajectory files.
- After the conclusion of all tasks in the array, task 0 will collect all trajectory files
    into a single converged file and unconverged file in the root directory.
- The job array directories and the log files can be kept or deleted.