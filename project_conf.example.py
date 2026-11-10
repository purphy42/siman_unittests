# -*- coding: utf-8 -*-
from __future__ import division, unicode_literals, absolute_import

NEW_BATCH = True
check_job = 0
PBS_PROCS = False
WALLTIME_LIMIT = True
path_to_phonopy = '/home/arseniy/.local/bin/phonopy'

CLUSTERS = {}
DEFAULT_CLUSTER = 'razor128'
AUTO_UPDATE_DB = True
PATH2PROJECT = '/siman_develop/'
PATH_TO_PROJECT_ON_COMP = '/home/arseniy/Desktop/work/siman_develop/'

CLUSTERS['razor128'] = {'address': 'razor128', 'vasp_com': 'srun vasp_std', 'homepath': '/home/a.burov', 'schedule': 'SLURM', 'corenum': 16, 'partition': 'zen4', 'scratch': '/scr/a.burov/', 'modules': 'source /etc/profile.d/modules.sh; module load zen4/vasp/5.4.4_openmpi_gcc_aocl; \nulimit -s unlimited\n'}

CLUSTERS['razor128_cp'] = {'address': 'razor128', 'vasp_com': 'srun vasp_std', 'homepath': '/home/a.burov/', 'schedule': 'SLURM', 'corenum': 16, 'partition': 'zen4', 'scratch': '/scr/a.burov/', 'modules': 'source /etc/profile.d/modules.sh; module load zen4/cp-vasp/6.3.2; \nulimit -s unlimited\n'}

CLUSTERS['razor128_res'] = {'address': 'razor128', 'vasp_com': 'srun vasp_std', 'homepath': '/home/a.burov/siman_develop/', 'schedule': 'SLURM', 'corenum': 16, 'partition': 'zen4', 'scratch': '/scr/a.burov/', 'modules': 'source /etc/profile.d/modules.sh; module load zen4/cp-vasp/6.3.2; \nulimit -s unlimited\n'}

PATH2POTENTIALS = '/home/arseniy/Desktop/vasp/potpaw_paw/potpaw_PBE_MPIE'
PATH2EDITOR = 'subl'
pmgkey = ''
mpkey = ''
PATH2DATABASE = '/home/arseniy/Desktop/work/siman_develop/'
PATH2JMOL = '/home/arseniy/Desktop/Jmol/Jmol-16.4.7-binary/jmol-16.4.7/jmol.sh'
PATH2VESTA = '/home/arseniy/Desktop/VESTA-gtk3/VESTA'
PATH2NEBMAKE = '/home/arseniy/Desktop/study/comp_chemistry/final_project/vts/nebmake.pl'
PATH2BADER = '/home/arseniy/Downloads/bader_lnx_64/bader'
geo_folder = './'
RAMDISK = None
EXCLUDE_NODES = 0
CIF2CELL = True
MANUALLY_ADDED = []
