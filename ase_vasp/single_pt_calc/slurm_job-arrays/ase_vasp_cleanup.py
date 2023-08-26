from ase.io.trajectory import Trajectory
from ase.io.ulm import InvalidULMFileError
import sys, os


job_array_len = int(sys.argv[1])
out_traj_prefix= sys.argv[2]

conv_traj_name = out_traj_prefix + "_ef.traj"
unconv_traj_name = out_traj_prefix + "_unconv.traj"

conv_traj = Trajectory(conv_traj_name, 'w')
unconv_traj = Trajectory(unconv_traj_name, 'w')


for i in range(job_array_len):
    try:
        conv_traj_i = Trajectory( os.path.join(str(i), conv_traj_name), 'r' )
        for image in conv_traj_i:
            conv_traj.write(image)
        conv_traj_i.close()
    except InvalidULMFileError:
        print(f"{os.path.join(str(i), 'Pt_ef.traj')} is not a valid ULM file. Skipping...")

    try:
        unconv_traj_i = Trajectory( os.path.join(str(i), unconv_traj_name), 'r' )
        for image in unconv_traj_i:
            unconv_traj.write(image)
        unconv_traj_i.close()
    except InvalidULMFileError:
        print(f"{os.path.join(str(i), 'Pt_unconv.traj')} is not a valid ULM file. Skipping...")

conv_traj.close()
unconv_traj.close()
