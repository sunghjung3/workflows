from ase.calculators.vasp import Vasp
from ase.io import Trajectory
from pymatgen.io.vasp.outputs import Vasprun

import sys, os

traj_file = sys.argv[1]
job_array_id = int(sys.argv[2])  # 0 to job_array_len - 1
job_array_len = int(sys.argv[3])


B_data = Trajectory(traj_file,'r')
B_data_ef = Trajectory('B_ef.traj','w')
B_data_unconv = Trajectory('B_unconv.traj','w')

lower_image = round( len(B_data) / job_array_len * job_array_id )
upper_image = round( len(B_data) / job_array_len * (job_array_id+1) )  # exclusive

calc = Vasp(prec = 'Medium',
            xc = 'PBE',
            setups = {'B':''},
            ediff = 1e-6,
            #ediffg = -0.01,
            #kpts = 28,
            kspacing = 0.16,
            lcharg = False,
            lwave = False,
            isym = 0,
            ispin = 1,
            ncore = 8,
            algo = 'VeryFast',
            lreal = 'Auto',
            #encut = 250,
            ismear = 1,
            sigma = 0.2,
            nsw = 0,
            ibrion = -1,
            isif = 3,
            #icharg = 1,
            #lorbit = 11,
            lasph = True,
            nelm = 40,
            #amix = 0.02,
            #bmix = 0.2,
            #lmaxmix = 4
            )

for k, atoms in enumerate(B_data[lower_image:upper_image]):
    status_file_name = "." + str(lower_image) + "_" + str(k + lower_image) + "_" + str(upper_image)
    with open(status_file_name, 'w') as f:
        pass  # create an empty file
    atoms.set_calculator(calc)
    atoms.get_potential_energy(force_consistent=True)
    if Vasprun("vasprun.xml").converged_electronic:
        B_data_ef.write(atoms)
    else:
        B_data_unconv.write(atoms)
    os.remove(status_file_name)
B_data.close()
B_data_ef.close()
B_data_unconv.close()
