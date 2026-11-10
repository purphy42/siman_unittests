# Copy this file to project_conf.py and fill in the paths on this computer.
# Do not commit project_conf.py.

NEW_BATCH = True
check_job = 0
PBS_PROCS = False
WALLTIME_LIMIT = True
path_to_phonopy = 'phonopy'

CLUSTERS = {}
DEFAULT_CLUSTER = 'razor128'
AUTO_UPDATE_DB = True
PATH2PROJECT = '/siman_unittests/'
PATH_TO_PROJECT_ON_COMP = '/absolute/path/to/this/clone'

CLUSTERS['razor128'] = {
    'address': 'cluster-hostname',
    'vasp_com': 'srun vasp_std',
    'homepath': '/home/USER',
    'schedule': 'SLURM',
    'corenum': 16,
    'partition': 'partition-name',
    'scratch': '/scr/USER/',
    'modules': 'module load vasp\nulimit -s unlimited\n',
}

CLUSTERS['razor128_cp'] = {
    'address': 'cluster-hostname',
    'vasp_com': 'srun vasp_std',
    'homepath': '/home/USER',
    'schedule': 'SLURM',
    'corenum': 16,
    'partition': 'partition-name',
    'scratch': '/scr/USER/',
    'modules': 'module load cp-vasp\nulimit -s unlimited\n',
}

CLUSTERS['razor128_res'] = {
    'address': 'cluster-hostname',
    'vasp_com': 'srun vasp_std',
    'homepath': '/home/USER/siman_unittests/',
    'schedule': 'SLURM',
    'corenum': 16,
    'partition': 'partition-name',
    'scratch': '/scr/USER/',
    'modules': 'module load vasp\nulimit -s unlimited\n',
}

PATH2POTENTIALS = '/path/to/potpaw_PBE'
PATH2EDITOR = 'subl'
pmgkey = ''
mpkey = ''
PATH2DATABASE = '/absolute/path/to/this/clone'
PATH2JMOL = '/path/to/jmol.sh'
PATH2VESTA = '/path/to/VESTA'
PATH2NEBMAKE = '/path/to/nebmake.pl'
PATH2BADER = '/path/to/bader'
geo_folder = './'
RAMDISK = None
EXCLUDE_NODES = 0
CIF2CELL = True
MANUALLY_ADDED = []
