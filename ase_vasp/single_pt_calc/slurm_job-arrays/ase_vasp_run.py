from ase.calculators.vasp import Vasp
from ase.io import Trajectory
from pymatgen.io.vasp.outputs import Vasprun

import sys, os

traj_file = sys.argv[1]
job_array_id = int(sys.argv[2])  # 0 to job_array_len - 1
job_array_len = int(sys.argv[3])
out_traj_prefix = sys.argv[4]


input_data = Trajectory(traj_file,'r')
conv_traj = Trajectory(out_traj_prefix + '_ef.traj','w')
unconv_traj = Trajectory(out_traj_prefix + '_unconv.traj','w')

lower_image = round( len(input_data) / job_array_len * job_array_id )
upper_image = round( len(input_data) / job_array_len * (job_array_id+1) )  # exclusive


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

for k, atoms in enumerate(input_data[lower_image:upper_image]):
    status_file_name = "." + str(lower_image) + "_" + str(k + lower_image) + "_" + str(upper_image)
    with open(status_file_name, 'w') as f:
        pass  # create an empty file
    atoms.set_calculator(calc)
    atoms.get_potential_energy(force_consistent=True)
    if Vasprun("vasprun.xml").converged_electronic:
        conv_traj.write(atoms)
    else:
        unconv_traj.write(atoms)
    os.remove(status_file_name)
input_data.close()
conv_traj.close()
unconv_traj.close()
