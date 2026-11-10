"""
SIMAN benchmark.

Set one tag below to True and run this file. dry_run writes the jobs locally
(run=0). run submits every job again (run=2), including calculations that
already finished or are still marked as running. read downloads and parses
the outputs. copy puts finished folders in siman_unittests/.

The calculation name is structure.purpose. purpose is sp, relax, dos, aimd,
eos, band, bader, phonon, polaron, occ, ncl, sol, or cp. The
constant-potential job is submitted by run only after the solvent job has a
Fermi level.
"""

import copy as pycopy
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from siman import header
from siman.calc_manage import add, add_des, choose_cluster, get_structure_from_matproj, res
from siman.core.cluster_run_script import complete_run, make_run, prepare_run
from siman.functions import push_to_server, run_on_server
from siman.database import read_database, write_database
from siman.geo import replic, scale_cell_uniformly
from siman.header import _update_configuration
from siman.set_functions import read_vasp_sets


VASP = 'razor128'       # zen4 partition, zen4/vasp
VASP_SOL = 'razor128_cp'  # zen4 partition, zen4/cp-vasp
FOLDER = 'siman_bench'
RESULT_DIR = 'siman_unittests'

dry_run = False  # write every job locally, run=0
run = False      # submit every job again, including finished ones, run=2
read = False     # download and parse the outputs
copy = False    # copy finished folders to siman_unittests/

LI_ID = 'mp-135'
MNO_ID = 'mp-19006'

# (structure, purpose, versions, cluster)
CALCS = [
    ('Li', 'sp', [1], VASP),
    ('Li', 'dos', [1], VASP),
    ('Li', 'eos', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], VASP),
    ('Li', 'band', [1], VASP),
    ('Li2', 'relax', [1], VASP),
    ('Li2', 'aimd', [1], VASP),
    ('Li2', 'bader', [1], VASP),
    ('Li2', 'phonon', [1], VASP),
    ('MnO', 'polaron', [1, 2], VASP),
    ('MnO', 'occ', [1], VASP),
    ('MnO', 'ncl', [1], VASP),
    ('Li_slab', 'sol', [1], VASP_SOL),
    ('Li_slab', 'cp', [1], VASP_SOL),
]


def calc_name(it, purpose):
    return f'{it}.{purpose}'


bench_sets = [
    ('sp', 'static', {
        'NSW': 0, 'IBRION': -1, 'ISTART': 0, 'ICHARG': 2,
        'ISPIN': 1, 'ENCUT': 400, 'ENAUG': 600, 'KSPACING': 0.5,
        'ISMEAR': 1, 'SIGMA': 0.1, 'PREC': 'Normal', 'ALGO': 'Normal',
        'NELM': 60, 'LREAL': False, 'LWAVE': False, 'LCHARG': True,
        'savefile': 'ocz',
    }, 'over'),
    ('relax', 'sp', {
        'NSW': 5, 'IBRION': 2, 'POTIM': 0.2, 'ISIF': 2, 'EDIFFG': -0.05,
        'LCHARG': False,
    }, 'over'),
    ('dos', 'sp', {
        'ISMEAR': -5, 'SIGMA': 0.05, 'LORBIT': 11, 'NEDOS': 2001,
        'NSW': 0, 'IBRION': -1, 'LCHARG': False, 'savefile': 'odz',
    }, 'over'),
    ('aimd', 'sp', {
        'IBRION': 0, 'MDALGO': 2, 'SMASS': 0, 'POTIM': 2, 'NSW': 10,
        'TEBEG': 300, 'TEEND': 300, 'ISIF': 2, 'ISMEAR': 1,
        'LWAVE': False, 'LCHARG': False, 'LREAL': False, 'savefile': 'ozt',
    }, 'over'),
    ('band', 'sp', {
        'ICHARG': 2, 'LORBIT': 11, 'NSW': 0, 'IBRION': -1,
        'ISMEAR': 1, 'SIGMA': 0.1, 'LWAVE': False, 'LCHARG': False,
        'savefile': 'ox',
        'k_band_structure': [
            5, ('G', 0, 0, 0), ('H', 0.5, -0.5, 0.5),
            ('N', 0, 0, 0.5), ('G', 0, 0, 0),
        ],
    }, 'over'),
    ('bader', 'sp', {
        'NSW': 0, 'IBRION': -1, 'LAECHG': '.TRUE.', 'LCHARG': '.TRUE.',
        'LWAVE': False, 'savefile': 'oca',
    }, 'over'),
    ('phonon', 'sp', {
        'IBRION': 8, 'NFREE': 2, 'POTIM': 0.015, 'NSW': 20,
        'EDIFF': 1e-6, 'ISIF': 2, 'LREAL': False, 'PREC': 'Accurate',
        'ISMEAR': 1, 'SIGMA': 0.1, 'LWAVE': False, 'LCHARG': False,
        'NPAR': None, 'savefile': 'ox',
    }, 'over'),
    ('eos', 'sp', {
        'NSW': 0, 'IBRION': -1, 'ISIF': 2, 'PREC': 'Accurate',
    }, 'over'),
    ('polaron', 'sp', {
        'ISPIN': 2, 'ISMEAR': 0, 'SIGMA': 0.05, 'LORBIT': 11,
        'LMAXMIX': 4, 'MAGMOM': {'Mn': 5.0, 'O': 0.6},
        'LDAU': True, 'LDAUTYPE': 2, 'LDAUL': {'Mn': 2, 'O': -1},
        'LDAUU': {'Mn': 3.9, 'O': 0.0}, 'LDAUJ': {'Mn': 0.0, 'O': 0.0},
        'LASPH': True, 'LDAUPRINT': 2,
    }, 'over'),
    ('occ', 'polaron', {
        'NSW': 0, 'IBRION': -1, 'LDAUPRINT': 2, 'savefile': 'oxz',
    }, 'over'),
    ('ncl', 'sp', {
        'ISPIN': 2, 'LNONCOLLINEAR': '.TRUE.', 'LSORBIT': '.FALSE.',
        'SAXIS': [0, 0, 1], 'LORBIT': 11, 'ISMEAR': 0, 'SIGMA': 0.05,
        'LMAXMIX': 4, 'MAGMOM': {'Mn': 5.0, 'O': 0.6},
        'NSW': 0, 'IBRION': -1, 'LWAVE': False, 'LCHARG': False,
        'savefile': 'oz',
    }, 'over'),
    ('sol', 'sp', {
        'ISMEAR': 0, 'SIGMA': 0.1, 'LSOL': True, 'ISOL': 2, 'EB_K': 89.78,
        'LDIPOL': True, 'IDIPOL': 3, 'DIPOL': [0.5, 0.5, 0.5],
        'LVTOT': True,
        'LCEP': False, 'ALGO': 'Normal', 'NELM': 120,
        'LWAVE': False, 'LCHARG': False,
        'savefile': 'oczl', 'savesol': 'phi,ss,ri',
    }, 'over'),
    ('cp', 'sol', {
        'LCEP': True, 'NESCHEME': 1, 'CAP_MAX': 7, 'FERMICONVERGE': 0.001,
        'NELM': 200,
    }, 'over'),
]


def load_project():
    here = os.path.dirname(os.path.abspath(__file__))
    os.chdir(here)
    _update_configuration(os.path.join(here, 'project_conf.py'))
    header.override_cluster_address = 1
    header.conv, header.varset, _ = read_database()
    read_vasp_sets(bench_sets, override_global=1)
    write_database()


def _li():
    return get_structure_from_matproj(
        mat_proj_id=LI_ID, it='Li', it_folder=FOLDER, ver=1)


def _mno():
    primitive = get_structure_from_matproj(
        mat_proj_id=MNO_ID, it='MnO', it_folder=FOLDER, ver=1)
    return replic(primitive, mul=(2, 1, 1))


def li_pair(li):
    """Two Li atoms, so relaxation, AIMD, Bader, and phonons have internal motion."""
    pair = replic(li, mul=(2, 1, 1))
    pair.name = 'Li2'
    return pair


def li_displaced(pair):
    """Shift one atom so the ionic relaxation has a nonzero force."""
    st = pycopy.deepcopy(pair)
    x = np.array(st.xcart, dtype=float)
    x[1, 0] += 0.15
    st.xcart = x
    st.xcart2xred()
    st.name = 'Li2_relax'
    return st


def with_moments(st):
    st = pycopy.deepcopy(st)
    st.magmom = [5.0 if el == 'Mn' else 0.6 for el in st.get_elements()]
    return st


def li_slab(li):
    """2x2x4 primitive Li with 12 A of vacuum along c."""
    slab = replic(li, mul=(2, 2, 4))
    slab = slab.add_vacuum(vector_i=2, thickness=12.0)
    slab.name = 'Li_slab'
    os.makedirs('structures/siman_bench', exist_ok=True)
    slab.write_poscar('structures/siman_bench/Li_slab.sol.POSCAR')
    return slab


def polaron_pair(mno):
    """Two Mn sites, each with one Mn-O axis expanded by 0.1 A."""
    mn = [i for i, el in enumerate(mno.get_elements()) if el == 'Mn']
    start = mno.localize_polaron_dist(mn[0], d=0.1, nn=6, axis=2)
    end = mno.localize_polaron_dist(mn[1], d=0.1, nn=6, axis=2)
    for st, center in ((start, mn[0]), (end, mn[1])):
        moments = [0.6 if el == 'O' else 5.0 for el in st.get_elements()]
        moments[center] = 4.0
        st.magmom = moments
    os.makedirs('structures/siman_bench', exist_ok=True)
    start.write_poscar('structures/siman_bench/MnO.polaron.1.POSCAR')
    end.write_poscar('structures/siman_bench/MnO.polaron.2.POSCAR')
    return start, end


def structures():
    li = _li()
    pair = li_pair(li)
    displaced = li_displaced(pair)
    slab = li_slab(li)
    mno = _mno()
    magnetic = with_moments(mno)
    start, end = polaron_pair(mno)
    eos = scale_cell_uniformly(li, scale_region=(-15, 15), n_scale_images=10)
    os.makedirs('structures/siman_bench', exist_ok=True)
    for version, st in enumerate(eos, start=1):
        st.write_poscar(f'structures/siman_bench/Li.eos.{version}.POSCAR')
    pair.write_poscar('structures/siman_bench/Li2.POSCAR')
    displaced.write_poscar('structures/siman_bench/Li2.relax.POSCAR')
    return {
        ('Li', 'sp'): {1: li},
        ('Li', 'dos'): {1: li},
        ('Li', 'eos'): {version: st for version, st in enumerate(eos, start=1)},
        ('Li', 'band'): {1: li},
        ('Li2', 'relax'): {1: displaced},
        ('Li2', 'aimd'): {1: pair},
        ('Li2', 'bader'): {1: pair},
        ('Li2', 'phonon'): {1: pair},
        ('MnO', 'polaron'): {1: start, 2: end},
        ('MnO', 'occ'): {1: magnetic},
        ('MnO', 'ncl'): {1: magnetic},
        ('Li_slab', 'sol'): {1: slab},
        ('Li_slab', 'cp'): {1: slab},
    }


def _db_calc(it, purpose, version):
    key = (it, purpose, version)
    if key not in header.db:
        return None
    return header.db[key]


def solvent_fermi():
    cl = _db_calc('Li_slab', 'sol', 1)
    value = getattr(cl, 'efermi', None) if cl is not None else None
    if value is None:
        return None
    return float(value)


def _finished(cl):
    return cl is not None and '4' in str(getattr(cl, 'state', ''))


def _running(cl):
    return cl is not None and '3' in str(getattr(cl, 'state', ''))


def _in_bench(cl):
    return cl is not None and FOLDER in str(getattr(cl, 'dir', '')).replace('\\', '/')


def _add_one(it, purpose, version, cluster, input_st, run_mode, update_set=None, force=False, overwrite=False):
    add_des(header.struct_des, it, FOLDER, des='SIMAN benchmark', override=True)
    existing = _db_calc(it, purpose, version)
    if not overwrite and (run_mode == 2 or force) and (_finished(existing) or _running(existing)):
        print(f'{calc_name(it, purpose)}.{version} is already {existing.state}')
        return
    if run_mode == 0 and _in_bench(existing) and not force:
        print(f'{calc_name(it, purpose)}.{version} is already written')
        return

    # A previous submission sets run_flag. With run=0, SIMAN then returns
    # without rewriting. Clear it so a multi-version batch can be rebuilt
    # and submitted once by the combined script.
    if overwrite and run_mode == 0 and existing is not None:
        existing.run_flag = False

    saved_run = header.siman_run
    saved_copy = header.copy_to_cluster_flag
    if run_mode == 0:
        header.siman_run = 0
        header.copy_to_cluster_flag = False
        # An older structure entry can put the job outside the benchmark folder.
        up = 'up2' if force or (existing is not None and not _in_bench(existing)) else 'up3'
        copy_to_server = False
    else:
        header.siman_run = 1
        header.copy_to_cluster_flag = True
        up = 'up2'
        copy_to_server = True
    try:
        add(
            it, purpose, version, it_folder=FOLDER,
            input_st=input_st, cluster=cluster, up=up,
            run=run_mode, copy_to_server=copy_to_server,
            update_set=update_set, check_job=0 if run_mode == 0 or overwrite else 1,
        )
    finally:
        header.siman_run = saved_run
        header.copy_to_cluster_flag = saved_copy


def remote_project(cluster_name):
    home = header.CLUSTERS[cluster_name]['homepath'].rstrip('/')
    project = header.PATH2PROJECT.strip('/')
    return home + '/' + project


def use_remote_project(text, cluster_name):
    remote = remote_project(cluster_name)
    return text.replace('//siman_develop', remote).replace('/home/a.burov//siman_develop', remote)


def write_all_versions(run_file, versions):
    """Keep one batch header and run every POSCAR version before the footer."""
    path = Path(run_file)
    text = path.read_text()
    marker = '#Basic run:\n'
    footer_marker = '#Footer section:'
    if marker not in text or footer_marker not in text:
        return
    head, rest = text.split(marker, 1)
    footer_at = rest.rfind(footer_marker)
    body = rest[:footer_at]
    footer = rest[footer_at + len(footer_marker):]
    if marker in body:
        body = body.split(marker, 1)[0]
    body = re.sub(r'(?m)^cp \d+\.(INCAR|KPOINTS) \1\n', '', body)
    match = re.search(r'cp (\d+)\.POSCAR POSCAR', body)
    if not match:
        return
    source = match.group(1)
    folder = path.parent
    blocks = []
    for version in versions:
        block = body.replace(f'cp {source}.POSCAR POSCAR', f'cp {version}.POSCAR POSCAR')
        block = block.replace(f'.{source}.log', f'.{version}.log')
        block = re.sub(rf'(?m)^mv (\S+) {source}\.', rf'mv \1 {version}.', block)
        extra = ''
        for filename in ('INCAR', 'KPOINTS'):
            if (folder / f'{version}.{filename}').is_file():
                extra += f'cp {version}.{filename} {filename}\n'
        block = block.replace(
            f'cp {version}.POSCAR POSCAR\n',
            f'cp {version}.POSCAR POSCAR\n{extra}',
        )
        blocks.append(marker + block.rstrip() + '\n')
    path.write_text(head + '\n'.join(blocks) + '\n' + footer_marker + footer)


def _submit_combined(it, purpose, cluster):
    """Upload one already combined batch file and submit it once."""
    name = calc_name(it, purpose)
    local = os.path.join(FOLDER, name)
    saved_copy = header.copy_to_cluster_flag
    header.copy_to_cluster_flag = True
    try:
        choose_cluster(cluster, None, None)
        files = [
            os.path.join(local, filename)
            for filename in os.listdir(local)
            if os.path.isfile(os.path.join(local, filename))
        ]
        remote_dir = header.project_path_cluster.rstrip('/') + '/' + local
        push_to_server(files, remote_dir, header.cluster_address)
        batch = remote_dir + '/' + name + '.run'
        run_on_server('chmod +x ' + batch, header.cluster_address)
        prepare_run()
        make_run(local, header.schedule_system, batch)
        complete_run()
        run_on_server('./run', header.CLUSTER_ADDRESS)
    finally:
        header.copy_to_cluster_flag = saved_copy


def _snapshot_inputs(it, purpose, version):
    """One folder is shared, so keep each version's INCAR and KPOINTS."""
    folder = Path(FOLDER) / calc_name(it, purpose)
    for filename in ('INCAR', 'KPOINTS'):
        source = folder / filename
        if source.is_file():
            (folder / f'{version}.{filename}').write_text(source.read_text())


def launch(run_mode):
    """run_mode 0 writes the jobs. run_mode 2 submits them."""
    load_project()
    built = structures()
    fermi = solvent_fermi()
    for it, purpose, versions, cluster in CALCS:
        if purpose == 'cp' and fermi is None:
            print(f'{calc_name(it, purpose)} waits for the Fermi level of Li_slab.sol')
            continue
        update = {'TARGETMU': fermi} if purpose == 'cp' else None
        # Several versions share one batch file. Write them all, then submit once.
        # run=2 overwrites finished and running jobs instead of leaving them.
        write_mode = 0 if len(versions) > 1 else run_mode
        for version in versions:
            _add_one(
                it, purpose, version, cluster,
                built[(it, purpose)][version], write_mode, update_set=update,
                force=len(versions) > 1, overwrite=run_mode == 2,
            )
            if len(versions) > 1:
                _snapshot_inputs(it, purpose, version)
        run_file = Path(FOLDER) / calc_name(it, purpose) / (calc_name(it, purpose) + '.run')
        if run_file.is_file():
            run_file.write_text(use_remote_project(run_file.read_text(), cluster))
            if len(versions) > 1:
                write_all_versions(run_file, versions)
        if len(versions) > 1 and run_mode == 2:
            _submit_combined(it, purpose, cluster)
    _drop_placeholder()
    write_database()
    _drop_placeholder()


def _drop_placeholder():
    """Calculation() also stores the latest object under (None, '0', 1)."""
    placeholder = (None, '0', 1)
    if dict.__contains__(header.db, placeholder):
        dict.__delitem__(header.db, placeholder)
    import dbm
    import shelve
    database = dbm.open(header.PATH2SHELVE_DB, 'c')
    with shelve.Shelf(database, protocol=3) as stored:
        key = str(placeholder)
        if key in stored:
            del stored[key]


def fermi_from_outcar():
    """Fermi level from the finished solvent OUTCAR when the database has none."""
    for folder in (RESULT_DIR, FOLDER):
        path = Path(folder) / 'Li_slab.sol' / '1.OUTCAR'
        if not path.is_file():
            continue
        value = None
        for line in path.read_text(errors='replace').splitlines():
            if 'E-fermi' in line:
                value = float(line.split()[2])
            elif line.strip().startswith('Fermi energy'):
                value = float(line.split()[-1])
        if value is not None:
            return value
    return None


def relaunch(purposes):
    """Rewrite and submit selected jobs, including ones that already finished."""
    from siman.calc_manage import smart_structure_read

    load_project()
    fermi = solvent_fermi()
    if fermi is None:
        fermi = fermi_from_outcar()
    sources = {
        'phonon': 'siman_bench/Li2.phonon/1.POSCAR',
        'sol': 'siman_bench/Li_slab.sol/1.POSCAR',
        'cp': 'siman_bench/Li_slab.sol/1.POSCAR',
    }
    for it, purpose, versions, cluster in CALCS:
        if purpose not in purposes:
            continue
        if purpose == 'cp' and fermi is None:
            print(f'{calc_name(it, purpose)} waits for the Fermi level of Li_slab.sol')
            continue
        structure = smart_structure_read(sources[purpose])
        update = {'TARGETMU': fermi} if purpose == 'cp' else None
        for version in versions:
            _add_one(
                it, purpose, version, cluster, structure, 2,
                update_set=update, overwrite=True,
            )
        run_file = Path(FOLDER) / calc_name(it, purpose) / (calc_name(it, purpose) + '.run')
        if run_file.is_file():
            run_file.write_text(use_remote_project(run_file.read_text(), cluster))
    _drop_placeholder()
    write_database()
    _drop_placeholder()


def read_results():
    load_project()
    for it, purpose, versions, cluster in CALCS:
        if _db_calc(it, purpose, versions[0]) is None:
            print(f'{calc_name(it, purpose)} has not been created')
            continue
        try:
            res(it, purpose, versions, up='up2', check_job=1, show='', cluster="razor128_res")
        except Exception as error:
            print(f'{calc_name(it, purpose)} was not read: {error}')
        for version in versions:
            cl = _db_calc(it, purpose, version)
            if cl is None:
                continue
            energy = getattr(cl, 'e0', None)
            fermi = getattr(cl, 'efermi', None)
            print(f'{calc_name(it, purpose)}.{version}  state={cl.state}  E={energy}  Fermi={fermi}')


def _has_output(path):
    try:
        names = os.listdir(path)
    except FileNotFoundError:
        return False
    markers = ('OUTCAR', 'OSZICAR', 'CONTCAR', 'DOSCAR', 'vasprun', 'PHI', 'XDATCAR', 'AECCAR')
    return any(any(marker in name for marker in markers) for name in names)


def copy_results():
    load_project()
    os.makedirs(RESULT_DIR, exist_ok=True)
    copied = set()
    for it, purpose, versions, _cluster in CALCS:
        name = calc_name(it, purpose)
        cl = _db_calc(it, purpose, versions[0])
        if cl is None or not getattr(cl, 'dir', None):
            print(f'{name} has not been created')
            continue
        source = cl.dir.rstrip('/')
        if source in copied:
            continue
        copied.add(source)
        if not _has_output(source):
            print(f'{name} has no output in {source}')
            continue
        destination = os.path.join(RESULT_DIR, name)
        shutil.copytree(source, destination, dirs_exist_ok=True)
        print(f'{source} -> {destination}')


if __name__ == '__main__':
    if dry_run:
        launch(0)
    if run:
        launch(2)
    if read:
        read_results()
    if copy:
        copy_results()
