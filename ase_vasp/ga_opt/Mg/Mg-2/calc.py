from ase.calculators.vasp import Vasp
import os


def vasp_calc(vasp_dir):
    # Element order: Fe   Co    C    N    X
    calc = Vasp(
                # System Settings
                xc='PBE',
                directory=vasp_dir,

                # Electronic Structure
                encut=520,
                prec='Normal',
                algo='Normal',
                lreal='Auto',
                ediff=1e-6,
                ismear=0,
                sigma=0.05,
                nelm=100,
                isym=0,
                
                # Ionic Relaxation
                ediffg=-0.01,
                isif=3,
                ibrion=2,
                nsw=1000,
                potim=0.1,
                
                # Magnetism and Spin
                ispin=2,
                magmom= 16*[2.0] + 16*[2.0] + 96*[0.0] + 96*[0.0] + 2*[0.0],  # 16*2.0 16*2.0 96*0 96*0 2*0
                
                # van der Waals
                ivdw=11,
                
                # DFT+U settings
                ldau=True,
                ldautype=2,
                ldauu=[7.0, 3.0, 0.0, 0.0, 0.0],
                lmaxmix=4,
                
                # kpoints
                kpts = [6,6,3],

                # Parallelization
                ncore=16,
            )
    return calc

def existing_contcar(vasp_dir):
    contcar_path = os.path.join(vasp_dir, "CONTCAR")
    if os.path.exists(contcar_path) and os.path.getsize(contcar_path) > 0:
        return contcar_path
    else:
        return None

if __name__ == '__main__':
    import sys

    from ase.io import read, trajectory
    from pymatgen.io.vasp.outputs import Vasprun


    fname = sys.argv[1]
    vasp_dir = f'./{fname[:-5]}-vasp'

    # check if we can resume from a previously terminated job
    contcar_path = existing_contcar(vasp_dir)
    if contcar_path is not None:
        a = read(contcar)
        traj = trajectory.Trajectory(fname, 'w')
        traj.write(a)
        traj.close()

    print(f'Now relaxing {fname}')
    a = read(fname)

    a.calc = vasp_calc(vasp_dir)
    a.info['key_value_pairs']['raw_score'] = -a.get_potential_energy()

    if Vasprun(os.path.join(vasp_dir, "vasprun.xml")).converged_electronic:
        print("Converged!")
    else:
        print("Failed to converge.")

    write(fname[:-5] + '_done.traj', a)
    print(f'Done relaxing {fname}')
