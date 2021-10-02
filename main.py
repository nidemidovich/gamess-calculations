from job import Job
from utils import read_config, read_sys_input


# calc params
opt_params = {
    'contrl': {
        'scftyp': 'rhf',
        'runtyp': 'optimize',
        'units': 'angs',
        'maxit': 199
    },
    'system': {'mwords': 500},
    'basis': {
        'gbasis': 'n31',
        'ngauss': 6,
        'ndfunc': 1
    },
    'statpt': {
        'opttol': 0.0001,
        'nstep': 1000
    },
    'guess': {'guess': 'huckel'},
    'dft': {'dfttyp': 'b3lyp'}
}
nrg_params = {
    'contrl': {
        'scftyp': 'rhf',
        'runtyp': 'energy',
        'tddft': 'excite',
        'units': 'angs'
    },
    'system': {
        'mwords': 5000,
        'memddi': 2000
    },
    'basis': {
        'gbasis': 'n31',
        'ngauss': 6,
        'ndfunc': 1,
        'npfunc': 1
    },
    'guess': {'guess': 'huckel'},
    'dft': {'dfttyp': 'b3lyp'}
}

# preparation
config = read_config()
material, material_name = read_sys_input()

# main
job = Job(
    material=material,
    material_name=material_name,
    config=config,
    opt_params=opt_params,
    nrg_params=nrg_params
)
job.job()
