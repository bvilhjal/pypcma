"""
SLURM related code
"""

import os

def generate_slurm_script(file_name, run_id, command_str, queue_id=None, out_file=None,
                          err_file=None, walltime='24:00:00', num_cores=1, max_memory=4000,
                          work_dir=None, email=None):
    
    
    with open(file_name, 'w') as f:
        f.write('#!/bin/bash\n') 
        f.write('#SBATCH -c %d\n' % num_cores)
        if work_dir is not None:
            f.write('#SBATCH --workdir=%s\n' % work_dir)
        f.write('#SBATCH --job-name=%s\n' % run_id) 
        if out_file is not None:
            f.write('#SBATCH --output=%s\n' % out_file)
        if err_file is not None:
            f.write('#SBATCH --error=%s\n' % err_file)
        f.write('#SBATCH --time=%s\n' % walltime)
        f.write('#SBATCH --mem=%d\n' % max_memory)
        if queue_id is not None:
            f.write('#SBATCH -p %s\n' % queue_id)
        if email is not None:
            f.write('#SBATCH --mail-type=ALL\n')
            f.write('#SBATCH --mail-user=%s\n' % (email))
        f.write('\n')
        f.write(command_str)
        f.write('\n')
        
    return file_name


def submit_job(run_id, command, parameters=None, out_file=None, err_file=None,
               walltime='12:00:00', queue_id='normal', num_cores=4, max_memory=4000,
               script_dir='.', job_dir='.', email=None, delete_script_after_submit=False):
    """
    Submit a command to the cluster
    """
    script_file_name = '%s/%s.bash' % (script_dir, run_id)
    command_str = command 
    if parameters is not None:
        for k in parameters:
            command_str = '%s --%s=%s'(command_str, k, str(parameters[k]))
    
    file_name = generate_slurm_script(script_file_name, run_id, command_str, out_file=out_file,
                                      err_file=err_file, walltime=walltime, queue_id=queue_id,
                                      num_cores=num_cores, max_memory=max_memory, work_dir=job_dir, email=email)
    os.system('sbatch %s' % file_name)
    
    if delete_script_after_submit:
        os.remove(file_name)
