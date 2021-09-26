optimization_params = {
    'contrl': {
        'scftyp': 'rhf',
        'runtyp': 'optimize',
        'units': 'angs'
    },
    'system': {'mwords': 50},
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
    'dft': {'dfttyp': 'b31yp'}
}

energy_params = {
    'contrl': {
        'scftyp': 'rhf',
        'runtyp': 'energy',
        'units': 'angs'
    },
    'system': {'mwords': 50},
    'basis': {
        'gbasis': 'n31',
        'ngauss': 6,
        'ndfunc': 1
    },
    'guess': {'guess': 'huckel'},
    'dft': {'dfttyp': 'b31yp'}
}
