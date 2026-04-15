# script to submit each job to sumner2 #

import os

shell_scripts = [i for i in os.listdir() if i.endswith('.sh')]

for i in shell_scripts:
	os.system(f'sbatch {i}')
