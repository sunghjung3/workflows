from ase.ga.pbs_queue_run import PBSQueueRun
from ase.io import write
from subprocess import Popen, PIPE
import os


class SLURMQueueRun(PBSQueueRun):
    """ASE PBSQueueRun with SLURM sucks."""
    def __init__(self, data_connection, tmp_folder, job_prefix,
                 n_relax, n_ga, job_template_generator,
                 qsub_command='sbatch', qstat_command='squeue',
                 find_neighbors=None, perform_parametrization=None):
        super(SLURMQueueRun, self).__init__(data_connection, tmp_folder,
                                            job_prefix, n_ga,
                                            job_template_generator,
                                            qsub_command, qstat_command,
                                            find_neighbors, perform_parametrization)
        self.n_relax = n_relax
        self.qstat_flags = '-o "%.18i %.100j %.8u %.2t %.10M %.6D %R"'  # full job name (100 chars)
    
    def relax(self, a):
        """Copy from parent class, but returns SLURM submission exit code"""
        self.__cleanup__()
        self.dc.mark_as_queued(a)
        if not os.path.isdir(self.tmp_folder):
            os.mkdir(self.tmp_folder)
        fname = '{}/cand{}.traj'.format(self.tmp_folder,
                                        a.info['confid'])
        write(fname, a)
        job_name = '{}_{}'.format(self.job_prefix, a.info['confid'])
        with open('tmp_job_file.job', 'w') as fd:
            fd.write(self.job_template_generator(job_name, fname))
        c = os.system(f'{self.qsub_command} tmp_job_file.job')
        return c  # 0 if successful

    def enough_jobs_running_ga(self):
        return super().enough_jobs_running()
        
    def enough_jobs_running_relax(self):
        return self.number_of_jobs_running() >= self.n_relax
    
    def relevant_jobs(self):
        self.__cleanup__()
        p = Popen([f'`which {self.qstat_command}` -u `whoami` {self.qstat_flags}'],
                  shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                  close_fds=True, universal_newlines=True)
        fout = p.stdout
        lines = fout.readlines()
        myjobs = []
        for line in lines:
            if line.find(self.job_prefix + '_') != -1:
                myjobs.append(line)
        return myjobs

    def number_of_jobs_running(self):
        return len(self.relevant_jobs())

    def is_running(self, aid):
        job_name = '{}_{}'.format(self.job_prefix, aid)
        for job in self.relevant_jobs():
            if job_name in job:
                return True
        return False
