from ase.io.trajectory import Trajectory
import sys, os


job_array_len = int(sys.argv[1])

B_data_ef = Trajectory("B_ef.traj", 'w')
B_data_unconv = Trajectory("B_unconv.traj", 'w')


for i in range(job_array_len):
    B_data_ef_i = Trajectory( os.path.join(str(i), "B_ef.traj"), 'r' )
    for image in B_data_ef_i:
        B_data_ef.write(image)
    B_data_ef_i.close()

    B_data_unconv_i = Trajectory( os.path.join(str(i), "B_unconv.traj"), 'r' )
    for image in B_data_unconv_i:
        B_data_unconv.write(image)
    B_data_unconv_i.close()


B_data_ef.close()
B_data_unconv.close()
