"""VASP file checks. These read OUTCAR, vasprun.xml, and the other saved VASP files.

A finished calculation has a closed OUTCAR, an OSZICAR ionic or MD summary whose
final E0 matches that OUTCAR, and the extra file that belongs to the job:
DOSCAR points for a DOS, a line-mode vasprun for bands, an XDATCAR frame per MD
step, or a NGX NGY NGZ grid for Bader. An unfinished or missing calculation is
a failure. A method whose own extra file was not written, such as LOCPOT when
LVTOT is off, is a skip.
"""

import math
import re
import traceback
from pathlib import Path


HERE = Path(__file__).resolve().parent


# Li.eos is scaled from -15% to +15% in ten steps. linspace does not include 0%.
EOS_VERSIONS = range(1, 11)


def planned():
    """One result per calculator, each EOS and polaron image, and each VASP method."""
    count = len(CALCULATORS) + len(tuple(EOS_VERSIONS)) + 2 + len(VASP_METHODS) + len(LOCAL_VASP)
    folder = HERE / 'Li.eos'
    if all((folder / f'{version}.OUTCAR').is_file() for version in EOS_VERSIONS):
        count += 1
    return count


# Analysis methods. Each one is recorded even when its VASP files are not here yet.
# label, job folder, required files, what the report line checks
VASP_METHODS = [
    ('read_results()', 'Li.sp', ['1.OUTCAR'], 'parse the static total energy from OUTCAR'),
    ('read_oszicar()', 'Li.sp', ['1.OSZICAR'], 'parse the electronic and ionic steps'),
    ('dos()', 'Li.dos', ['1.OUTCAR', '1.DOSCAR'], 'plot the total DOS from DOSCAR'),
    ('get_band_info()', 'Li.band', ['1.OUTCAR'], 'read the band edges from OUTCAR'),
    ('band()', 'Li.band', ['1.vasprun.xml', 'KPOINTS'], 'plot the line-mode band structure'),
    ('plot_energy_force()', 'Li2.relax', ['1.OUTCAR'], 'plot energy against the maximum force'),
    ('plot_energy_step()', 'Li2.relax', ['1.OUTCAR'], 'plot the energy change at each ionic step'),
    ('plot_energy_conv()', 'Li.sp', ['1.OUTCAR'], 'plot the electronic convergence'),
    ('plot_energy_total_step()', 'Li2.aimd', ['1.OUTCAR'], 'plot the total MD energy'),
    ('plot_md_energies()', 'Li2.aimd', ['1.OUTCAR'], 'estimate the MD energy plateau'),
    ('e0_fu()', 'Li.sp', ['1.OUTCAR'], 'energy per formula unit'),
    ('e0_at()', 'Li.sp', ['1.OUTCAR'], 'energy per atom'),
    ('show_force()', 'Li2.relax', ['1.OUTCAR'], 'print the maximum force series'),
    ('plot_locpot()', 'Li_slab.sol', ['1.OUTCAR', '1.LOCPOT'], 'plot the planar electrostatic potential'),
    ('locpot_z_profile()', 'Li_slab.sol', ['1.OUTCAR', '1.LOCPOT'], 'average LOCPOT along c'),
    ('bader()', 'Li2.bader', ['1.OUTCAR', '1.CHGCAR'], 'run the local Bader analysis on CHGCAR'),
    ('get_occ_mat()', 'MnO.occ', ['1.OUTCAR'], 'return one onsite occupation matrix'),
    ('msd()', 'Li2.aimd', ['1.OUTCAR', '1.XDATCAR'], 'mean-square displacement of Li'),
    ('plot_cp_convergence()', 'Li_slab.cp', ['1.OUTCAR', '1.OSZICAR'], 'plot electron count, mu, and GCE'),
]

# Writers and parsers that use the Li POSCAR and do not need a finished VASP run.
LOCAL_VASP = [
    ('make_incar()', 'write an INCAR for the Li cell'),
    ('make_kpoints_file()', 'write a KPOINTS file for the Li cell'),
    ('update_incar()', 'change ENCUT in that INCAR'),
    ('add_potcar()', 'copy the Li POTCAR into the job folder'),
    ('calculate_nbands()', 'count valence electrons from the Li POTCAR'),
    ('check_kpoints()', 'choose a k-point mesh from KSPACING'),
    ('set_occ_mat()', 'store a spin-polarized d occupation matrix'),
    ('get_occ_mat()', 'read that occupation matrix back'),
    ('write_occmatrix()', 'write OCCMATRIX for one d shell'),
    ('vasp_dipole_center()', 'dipole center before a charge density exists'),
    ('determine_filenames()', 'name CHGCAR and AECCAR next to OUTCAR'),
    ('visualize_grid() vesta', 'convert a charge grid to a cube and open it in VESTA'),
    ('visualize_grid() jmol', 'convert a charge grid to a cube and open it in Jmol'),
]

# job folder, files that must exist before parsing, what the report line checks
CALCULATORS = [
    ('Li.sp', ['1.OUTCAR', '1.CONTCAR', '1.OSZICAR'], 'static energy, final geometry, and OSZICAR'),
    ('Li.dos', ['1.OUTCAR', '1.DOSCAR'], 'DOS total energy and the DOSCAR file'),
    ('Li.band', ['1.OUTCAR', '1.vasprun.xml', 'KPOINTS'], 'band-structure energy, vasprun, and the k-path'),
    ('Li2.relax', ['1.OUTCAR', '1.CONTCAR'], 'relaxed energy and the final CONTCAR'),
    ('Li2.aimd', ['1.OUTCAR', '1.XDATCAR'], 'AIMD energy series and the trajectory'),
    ('Li2.bader', ['1.OUTCAR', '1.CHGCAR', '1.AECCAR0', '1.AECCAR2'], 'charge density files used by Bader'),
    ('Li2.phonon', ['1.OUTCAR'], 'finite-difference phonon energy'),
    ('MnO.occ', ['1.OUTCAR'], 'DFT+U onsite occupation matrices'),
    ('MnO.ncl', ['1.OUTCAR'], 'noncollinear moments along x, y, and z'),
    ('Li_slab.sol', ['1.OUTCAR', '1.OSZICAR'], 'solvent energy, Fermi level, and solvation correction'),
    ('Li_slab.cp', ['1.OUTCAR', '1.OSZICAR'], 'constant-potential mu and grand-canonical energy'),
]


def _nonempty(path):
    return path.is_file() and path.stat().st_size > 0


def _refusal(path):
    """The sentence inside the last VASP error box."""
    lines = path.read_text(errors='replace').splitlines()
    start = None
    for index, line in enumerate(lines):
        if 'VASP internal routines' in line or 'Gamma-point-only' in line or 'I REFUSE' in line:
            start = index
    if start is None:
        return 'OUTCAR has no timing summary'
    prose = []
    for line in lines[start:]:
        if set(line.strip()) <= set('-') and len(line.strip()) > 20 and prose:
            break
        text = ' '.join(line.strip().strip('|').split())
        if not text or '###' in text:
            continue
        if text[:1] in 'ER' and text.count(' ') > 8:
            continue
        prose.append(text)
    return ' '.join(prose)[:500]


def finished_outcar(folder, version):
    """Return the completed OUTCAR, or a reason the calculation is not usable."""
    path = folder / f'{version}.OUTCAR'
    if _nonempty(path) and 'General timing' in path.read_text(errors='replace'):
        return path, None
    unfinished = folder / f'{version}.OUTCAR_unfinished'
    broken = path if _nonempty(path) else unfinished if unfinished.is_file() else None
    if broken is not None:
        return None, _refusal(broken)
    if not folder.is_dir():
        return None, 'folder is missing'
    return None, 'no finished OUTCAR'


def _last_ecorr(text):
    value = None
    for line in text.splitlines():
        if 'Ecorr_sol' in line:
            value = float(line.split()[-1])
    return value


_IONIC = re.compile(r'^\s*(\d+)\s+F=\s+\S+\s+E0=\s+([-+.\dEe]+)')
_MD = re.compile(
    r'^\s*(\d+)\s+T=\s+([-+.\dEe]+)\s+E=\s+\S+\s+F=\s+\S+\s+E0=\s+([-+.\dEe]+)'
)


def _is_int_token(token):
    try:
        int(token)
    except ValueError:
        return False
    return True


def oszicar_records(path):
    """Ionic or MD summaries from a finished OSZICAR.

    A static or relaxed run ends each ionic step with ``N F= ... E0=``.
    An AIMD run ends each step with ``N T= ... E0=``. The step numbers
    start at 1 and have no gaps.
    """
    text = path.read_text(errors='replace')
    if 'DAV:' not in text and 'RMM:' not in text:
        raise AssertionError('OSZICAR has no electronic steps')
    md = []
    ionic = []
    for line in text.splitlines():
        matched = _MD.match(line)
        if matched:
            md.append((int(matched.group(1)), float(matched.group(2)), float(matched.group(3))))
            continue
        matched = _IONIC.match(line)
        if matched:
            ionic.append((int(matched.group(1)), float(matched.group(2))))
    steps = md or ionic
    if not steps:
        raise AssertionError('OSZICAR has no ionic summary')
    numbers = [step[0] for step in steps]
    if numbers != list(range(1, len(numbers) + 1)):
        raise AssertionError(f'OSZICAR steps are {numbers}')
    return ('md' if md else 'ionic'), steps


def outcar_e0(path):
    """Last energy(sigma->0) written in a finished OUTCAR."""
    value = None
    for line in path.read_text(errors='replace').splitlines():
        if 'energy(sigma->0)' in line:
            value = float(line.split()[-1])
    if value is None:
        raise AssertionError('OUTCAR has no energy(sigma->0)')
    return value


def contcar_natom(path):
    """Atom count from a VASP 4 or VASP 5 CONTCAR."""
    lines = [line.strip() for line in path.read_text(errors='replace').splitlines() if line.strip()]
    if len(lines) < 7:
        raise AssertionError('CONTCAR is shorter than a structure header')
    species = lines[5].split()
    counts = species if all(_is_int_token(token) for token in species) else lines[6].split()
    if not counts or not all(_is_int_token(token) for token in counts):
        raise AssertionError('CONTCAR atom counts are not integers')
    total = sum(int(token) for token in counts)
    if total <= 0:
        raise AssertionError('CONTCAR has no atoms')
    return total


def charge_grid(path):
    """NGX NGY NGZ line of a CHGCAR or AECCAR."""
    for line in path.read_text(errors='replace').splitlines():
        parts = line.split()
        if len(parts) == 3 and all(part.isdigit() for part in parts):
            grid = tuple(int(part) for part in parts)
            if all(size > 0 for size in grid):
                return grid
    raise AssertionError(f'{path.name} has no NGX NGY NGZ line')


def doscar_summary(path):
    """NEDOS and Fermi energy from the total-DOS header, then that many points."""
    lines = path.read_text(errors='replace').splitlines()
    header = None
    for index, line in enumerate(lines):
        parts = line.split()
        if len(parts) != 5:
            continue
        try:
            nedos = int(float(parts[2]))
            efermi = float(parts[3])
        except ValueError:
            continue
        if nedos > 10:
            header = (index, nedos, efermi)
            break
    if header is None:
        raise AssertionError('DOSCAR has no NEDOS header')
    start, nedos, efermi = header
    points = [line for line in lines[start + 1:start + 1 + nedos] if line.strip()]
    if len(points) != nedos:
        raise AssertionError(f'DOSCAR has {len(points)} points, header says {nedos}')
    return nedos, efermi


def xdatcar_frames(path):
    """Number of sequential Direct configurations in an XDATCAR."""
    numbers = [int(value) for value in re.findall(r'Direct configuration=\s*(\d+)', path.read_text(errors='replace'))]
    if not numbers or numbers != list(range(1, len(numbers) + 1)):
        raise AssertionError(f'XDATCAR configurations are {numbers}')
    return len(numbers)


def vasprun_text(path):
    """A finished vasprun.xml starts as XML and closes the modeling element."""
    text = path.read_text(errors='replace')
    if not text.lstrip().startswith('<?xml'):
        raise AssertionError('vasprun.xml does not start with an XML declaration')
    if '<calculation>' not in text or '</modeling>' not in text:
        raise AssertionError('vasprun.xml has no closed calculation')
    return text


def _potim(folder):
    """MD time step in fs. A bare OUTCAR parse does not fill the SIMAN set."""
    for name in ('INCAR', '1.INCAR'):
        path = folder / name
        if not path.is_file():
            continue
        for line in path.read_text(errors='replace').splitlines():
            if line.strip().upper().startswith('POTIM'):
                return float(line.split('=', 1)[1].split()[0])
    return None


def parse_outcar(folder, version):
    from siman.calculators.vasp import CalculationVasp
    out = folder / f'{version}.OUTCAR'
    cl = CalculationVasp(iid=(folder.name, 'bench', version), output=str(out))
    cl.dir = str(folder) + '/'
    cl.read_results()
    osz = folder / f'{version}.OSZICAR'
    if osz.is_file():
        cl.read_oszicar(path=str(osz), warn=False)
    return cl


def _force_values(cl):
    values = []
    for item in getattr(cl, 'maxforce_list', []) or []:
        if isinstance(item, (list, tuple)):
            values.append(float(item[-1]))
        else:
            values.append(float(item))
    return values


def _remember(problems, action):
    try:
        action()
    except (AssertionError, ValueError, IndexError, OSError) as error:
        problems.append(str(error))


def _problems(name, cl, outcar):
    """Physical checks that a finished OUTCAR still has to satisfy."""
    text = outcar.read_text(errors='replace')
    folder = outcar.parent
    version = outcar.name.split('.', 1)[0]
    problems = []
    energy = getattr(cl, 'e0', None)
    try:
        if energy is None or not math.isfinite(float(energy)):
            problems.append('no energy')
        elif abs(float(energy) - outcar_e0(outcar)) > 1e-4:
            problems.append('parsed energy does not match energy(sigma->0)')
    except (TypeError, ValueError, AssertionError) as error:
        problems.append(str(error))
    if 'aborting loop because EDIFF is reached' not in text:
        problems.append('the electronic loop did not reach EDIFF')
    if 'General timing' not in text:
        problems.append('OUTCAR has no timing summary')

    contcar = folder / f'{version}.CONTCAR'
    if _nonempty(contcar):
        def geometry():
            count = contcar_natom(contcar)
            atoms = getattr(cl.end, 'natom', None)
            if atoms is not None and count != atoms:
                raise AssertionError(f'CONTCAR has {count} atoms, OUTCAR has {atoms}')
        _remember(problems, geometry)

    oszicar = folder / f'{version}.OSZICAR'
    if _nonempty(oszicar):
        def ionic():
            kind, steps = oszicar_records(oszicar)
            e0 = steps[-1][-1]
            if abs(e0 - float(cl.e0)) > 1e-4:
                raise AssertionError(f'OSZICAR E0 {e0} does not match OUTCAR {cl.e0}')
            if name == 'Li2.aimd':
                if kind != 'md' or len(steps) != 10:
                    raise AssertionError(f'AIMD OSZICAR has {len(steps)} {kind} steps')
            elif name == 'Li2.relax':
                if kind != 'ionic' or len(steps) < 2:
                    raise AssertionError(f'relaxation OSZICAR has {len(steps)} {kind} steps')
            elif kind != 'ionic' or len(steps) != 1:
                raise AssertionError(f'static OSZICAR has {len(steps)} {kind} steps')
        _remember(problems, ionic)

    if name == 'Li2.relax':
        if 'reached required accuracy' not in text:
            problems.append('the ionic relaxation did not reach the required accuracy')
        forces = _force_values(cl)
        if len(forces) < 2 or not forces[-1] < forces[0] or forces[-1] > 50:
            problems.append(f'forces did not fall below 50 meV/A: {forces}')
    if name == 'Li2.aimd':
        steps = len(getattr(cl, 'list_e_sigma0', []) or [])
        if steps != 10:
            problems.append(f'{steps} MD steps, expected 10')
        trajectory = folder / f'{version}.XDATCAR'
        if _nonempty(trajectory):
            def frames():
                count = xdatcar_frames(trajectory)
                if count != 10:
                    raise AssertionError(f'XDATCAR has {count} frames, expected 10')
            _remember(problems, frames)
    if name == 'Li.dos':
        doscar = folder / f'{version}.DOSCAR'
        if _nonempty(doscar):
            def dos_points():
                nedos, efermi = doscar_summary(doscar)
                if nedos != 2001:
                    raise AssertionError(f'DOSCAR NEDOS is {nedos}, INCAR asked for 2001')
                parsed = getattr(cl, 'efermi', None)
                if parsed is not None and abs(efermi - float(parsed)) > 1e-3:
                    raise AssertionError(f'DOSCAR Fermi {efermi} does not match OUTCAR {parsed}')
            _remember(problems, dos_points)
    if name == 'Li.band':
        def bands():
            xml = vasprun_text(folder / f'{version}.vasprun.xml')
            if '<eigenvalues>' not in xml:
                raise AssertionError('vasprun.xml has no eigenvalues')
            kpoints = (folder / 'KPOINTS').read_text(errors='replace')
            if 'Line-mode' not in kpoints:
                raise AssertionError('KPOINTS is not line-mode')
        _remember(problems, bands)
    if name == 'Li2.bader':
        def grids():
            for filename in (f'{version}.CHGCAR', f'{version}.AECCAR0', f'{version}.AECCAR2'):
                charge_grid(folder / filename)
        _remember(problems, grids)
    if name == 'MnO.occ':
        matrices = getattr(cl, 'occ_matrices', None) or {}
        elements = list(cl.end.get_elements()) if getattr(cl, 'end', None) is not None else []
        manganese = elements.count('Mn')
        if not matrices:
            problems.append('OUTCAR has no occupation matrices')
        elif manganese and len(matrices) != manganese:
            problems.append(f'{len(matrices)} occupation matrices for {manganese} Mn')
        else:
            for atom, matrix in matrices.items():
                if len(matrix) != 10:
                    problems.append(f'atom {atom} occupation matrix has {len(matrix)} rows, expected 10')
                    break
    if name == 'MnO.ncl':
        moments = getattr(cl, 'magmom_dir', {}) or {}
        atoms = getattr(cl.end, 'natom', None)
        missing = [
            axis for axis in ('x', 'y', 'z')
            if moments.get(axis) is None or atoms is None or len(moments.get(axis)) != atoms
        ]
        z = moments.get('z')
        if missing:
            problems.append('magmom ' + ', '.join(missing) + ' does not have one value per atom')
        elif max(abs(float(value)) for value in z) < 0.5:
            problems.append('noncollinear moments were not read')
    if name == 'Li_slab.sol':
        ecorr = _last_ecorr(text)
        if ecorr is None or ecorr >= -0.1:
            problems.append(f'solvation correction was not applied (Ecorr_sol={ecorr})')
        if getattr(cl, 'efermi', None) is None:
            problems.append('no Fermi level')
        if re.search(r'LVTOT\s*=\s*T\b', text) and not _nonempty(folder / f'{version}.LOCPOT'):
            problems.append('LVTOT is set but LOCPOT was not saved')
    return problems


def _parsed(report, name, cl, outcar):
    detail = (
        f'E={getattr(cl, "e0", None)} '
        f'Fermi={getattr(cl, "efermi", None)} '
        f'atoms={getattr(cl.end, "natom", None)}'
    )
    if name == 'MnO.ncl':
        z = (getattr(cl, 'magmom_dir', {}) or {}).get('z')
        largest = max(abs(float(value)) for value in z) if z is not None and len(z) else 0
        detail += f' max|m_z|={largest:.3f}'
    if name == 'MnO.occ':
        detail += f' occ_matrices={len(getattr(cl, "occ_matrices", {}) or {})}'
    if name == 'Li2.aimd':
        detail += f' md_steps={len(getattr(cl, "list_e_sigma0", []) or [])}'
    if name == 'Li_slab.sol':
        detail += f' Ecorr_sol={_last_ecorr(outcar.read_text(errors="replace"))}'
    if name == 'Li_slab.cp':
        detail += f' mu={getattr(cl, "mu_e", None)} GCE={getattr(cl, "gce", None)}'
    note = dict((item[0], item[2]) for item in CALCULATORS).get(name, '')
    if note:
        detail = f'{note}. {detail}'
    problems = _problems(name, cl, outcar)
    if problems:
        report.add('FAIL', name, detail + '. ' + '; '.join(problems))
    else:
        report.add('PASS', name, detail)


def _one(report, folder, version, label):
    out, problem = finished_outcar(folder, version)
    if problem:
        report.add('FAIL', label, problem)
        return None
    try:
        cl = parse_outcar(folder, version)
    except Exception:
        report.add('FAIL', label, traceback.format_exc().strip().splitlines()[-1])
        return None
    return cl


def run(report):
    import matplotlib
    matplotlib.use('Agg')
    report.section('VASP input')
    for label, note in LOCAL_VASP:
        try:
            detail = call_local(label)
        except (Exception, SystemExit) as error:
            report.add('FAIL', label, f'{note}. {type(error).__name__}: {error}')
            continue
        report.add('PASS', label, f'{note}. {detail}')

    report.section('VASP')
    for name, required, note in CALCULATORS:
        folder = HERE / name
        out, problem = finished_outcar(folder, 1)
        if problem:
            report.add('FAIL', name, f'{note}. {problem}')
            continue
        missing = [item for item in required if not _nonempty(folder / item)]
        if missing:
            report.add('FAIL', name, f'{note}. missing or empty ' + ', '.join(missing))
            continue
        cl = _one(report, folder, 1, name)
        if cl is not None:
            _parsed(report, name, cl, out)

    report.section('VASP Li.eos')
    folder = HERE / 'Li.eos'
    found = []
    for version in EOS_VERSIONS:
        cl = _one(report, folder, version, f'Li.eos.{version}')
        if cl is None:
            continue
        out, _problem = finished_outcar(folder, version)
        energy = getattr(cl, 'e0', None)
        volume = float(cl.end.vol)
        found.append((energy, volume))
        problems = _problems('Li.eos', cl, out)
        if problems:
            report.add('FAIL', f'Li.eos.{version}', '; '.join(problems))
        else:
            report.add(
                'PASS', f'Li.eos.{version}',
                f'static energy of one scaled volume. E={energy} V={volume:.3f}')
    if len(found) == len(tuple(EOS_VERSIONS)) and all(energy is not None for energy, _volume in found):
        energies = [energy for energy, _volume in found]
        volumes = [volume for _energy, volume in found]
        problems = []
        if any(volumes[i] >= volumes[i + 1] for i in range(len(volumes) - 1)):
            problems.append('the volumes are not in increasing order')
        if energies.index(min(energies)) in (0, len(energies) - 1):
            problems.append('the lowest energy is at an endpoint, so the volume window missed the minimum')
        static = parse_outcar(HERE / 'Li.sp', 1)
        static_volume = float(static.end.vol)
        if not volumes[0] < static_volume < volumes[-1]:
            problems.append(
                f'static volume {static_volume:.3f} is outside {volumes[0]:.3f}..{volumes[-1]:.3f}')
        detail = (
            'ten volumes from -15% to +15%, covering the static cell, '
            'with the minimum inside the window. '
            + ' '.join(f'{value:.4f}' for value in energies)
        )
        if problems:
            report.add('FAIL', 'Li.eos ten energies', detail + '. ' + '; '.join(problems))
        else:
            report.add('PASS', 'Li.eos ten energies', detail)

    report.section('VASP MnO.polaron')
    folder = HERE / 'MnO.polaron'
    images = {}
    for version in (1, 2):
        cl = _one(report, folder, version, f'MnO.polaron.{version}')
        if cl is None:
            continue
        out, _problem = finished_outcar(folder, version)
        images[version] = cl
        problems = _problems('MnO.polaron', cl, out)
        if version == 2 and 1 in images:
            if images[1].end.natom != cl.end.natom:
                problems.append('the two images have different atom counts')
            if abs(float(images[1].e0) - float(cl.e0)) < 1e-6:
                problems.append('the two distorted images have the same energy')
        if problems:
            report.add('FAIL', f'MnO.polaron.{version}', '; '.join(problems))
        else:
            report.add(
                'PASS', f'MnO.polaron.{version}',
                f'energy of polaron image {version}. E={cl.e0} atoms={cl.end.natom}')

    report.section('VASP methods')
    for label, job, required, note in VASP_METHODS:
        folder = HERE / job
        missing = [item for item in required if not _nonempty(folder / item)]
        if missing and label in ('plot_locpot()', 'locpot_z_profile()'):
            report.add('SKIP', label, f'{note}. LVTOT is off, so this job has no LOCPOT')
            continue
        if missing and label == 'plot_cp_convergence()':
            report.add('FAIL', label, f'{note}. the constant-potential job was not run')
            continue
        if missing:
            report.add('SKIP', label, f'{note}. missing ' + ', '.join(missing))
            continue
        try:
            detail = call_vasp_method(label, folder)
        except (Exception, SystemExit):
            report.add('FAIL', label, f'{note}. ' + traceback.format_exc().strip().splitlines()[-1])
            continue
        if isinstance(detail, tuple) and detail[0] == 'SKIP':
            report.add('SKIP', label, f'{note}. {detail[1]}')
            continue
        report.add('PASS', label, f'{note}. {detail}')


def project_setting(name):
    """Read one string setting from project_conf.py, or from the example."""
    import ast
    conf = HERE.parent / 'project_conf.py'
    if not conf.is_file():
        conf = HERE.parent / 'project_conf.example.py'
    if not conf.is_file():
        return None
    for line in conf.read_text(encoding='utf-8').splitlines():
        code = line.split('#', 1)[0].strip()
        if '=' not in code:
            continue
        key, _, raw = code.partition('=')
        if key.strip() != name:
            continue
        try:
            value = ast.literal_eval(raw.strip())
        except (SyntaxError, ValueError):
            return None
        return value if isinstance(value, str) and value else None
    return None


def call_vasp_method(label, folder):
    """Run one VASP analysis method. Matplotlib must not open a window."""
    import matplotlib
    matplotlib.use('Agg')
    if label == 'read_oszicar()':
        from siman.calculators.vasp import CalculationVasp
        cl = CalculationVasp(iid=(folder.name, 'bench', 1), output=str(folder / '1.OUTCAR'))
        cl.dir = str(folder) + '/'
        cl.read_oszicar(path=str(folder / '1.OSZICAR'), warn=False)
        kind, steps = oszicar_records(folder / '1.OSZICAR')
        if kind != 'ionic' or len(steps) != 1:
            raise AssertionError(f'static OSZICAR has {len(steps)} {kind} steps')
        if abs(steps[-1][1] - outcar_e0(folder / '1.OUTCAR')) > 1e-4:
            raise AssertionError('OSZICAR E0 does not match OUTCAR')
        electronic = len(getattr(cl, 'oszicar_electronic', []) or [])
        if electronic < 1:
            raise AssertionError('no electronic steps were parsed')
        return f'electronic={electronic} E0={steps[-1][1]}'
    cl = parse_outcar(folder, 1)
    if label == 'read_results()':
        if abs(float(cl.e0) - outcar_e0(folder / '1.OUTCAR')) > 1e-4:
            raise AssertionError('parsed energy does not match energy(sigma->0)')
        return f'E={cl.e0}'
    if label == 'dos()':
        from siman.dos_functions import plot_dos
        nedos, efermi = doscar_summary(folder / '1.DOSCAR')
        if abs(efermi - float(cl.efermi)) > 1e-3:
            raise AssertionError(f'DOSCAR Fermi {efermi} does not match OUTCAR {cl.efermi}')
        plot_dos(cl, dostype='total', show=0, savefile=True, path=str(folder))
        if not list(folder.glob('*.pdf')):
            raise AssertionError('DOS plot was not written')
        return f'NEDOS={nedos} Fermi={efermi:.4f}'
    if label == 'get_band_info()':
        cl.get_band_info(method='siman')
        xml = vasprun_text(folder / '1.vasprun.xml') if (folder / '1.vasprun.xml').is_file() else ''
        if xml and '<eigenvalues>' not in xml:
            raise AssertionError('vasprun.xml has no eigenvalues')
        return 'eigenvalues were read'
    if label == 'band()':
        xml = vasprun_text(folder / '1.vasprun.xml')
        if '<eigenvalues>' not in xml:
            raise AssertionError('vasprun.xml has no eigenvalues')
        kpoints = (folder / 'KPOINTS').read_text(errors='replace')
        if 'Line-mode' not in kpoints:
            raise AssertionError('KPOINTS is not line-mode')
        cl.path['kpoints'] = str(folder / 'KPOINTS')
        try:
            cl.band()
        except AttributeError as error:
            if 'savefig' not in str(error):
                raise
            raise AssertionError(
                'vasprun.xml has a line-mode band structure; '
                'band() saves the axes instead of the figure'
            ) from error
        return 'plotted'
    if label == 'plot_energy_force()':
        cl.plot_energy_force()
        return 'plotted'
    if label == 'plot_energy_step()':
        cl.plot_energy_step()
        return 'plotted'
    if label == 'plot_energy_conv()':
        cl.plot_energy_conv()
        return 'plotted'
    if label == 'plot_energy_total_step()':
        cl.plot_energy_total_step()
        return 'plotted'
    if label == 'plot_md_energies()':
        # read_results() does not copy POTIM into the set, which this plot uses as the time step.
        cl.set.params['POTIM'] = _potim(folder)
        frames = xdatcar_frames(folder / '1.XDATCAR')
        result = cl.plot_md_energies(show=False, smooth_window=1, tail_fraction=1.0)
        if not isinstance(result, dict) or len(result.get('time', [])) != frames:
            raise AssertionError('MD plot does not cover every XDATCAR frame')
        return f'frames={frames} slope={float(result["slope"]):.3e}'
    if label == 'e0_fu()':
        per_formula = float(cl.e0_fu(n_fu=1))
        if abs(per_formula - float(cl.e0)) > 1e-6:
            raise AssertionError(f'{per_formula} != {cl.e0}')
        return f'E_fu={per_formula}'
    if label == 'e0_at()':
        per_atom = float(cl.e0_at)
        expected = float(cl.e0) / cl.end.natom
        if abs(per_atom - expected) > 1e-6:
            raise AssertionError(f'{per_atom} != {expected}')
        return f'E_at={per_atom}'
    if label == 'show_force()':
        cl.show_force()
        forces = _force_values(cl)
        if len(forces) < 2 or not forces[-1] < forces[0] or forces[-1] > 50:
            raise AssertionError(f'forces did not fall below 50 meV/A: {forces}')
        return f'steps={len(forces)} last={forces[-1]:.0f} meV/A'
    if label == 'plot_locpot()':
        cl.plot_locpot(filename='benchmark_locpot')
        return 'plotted'
    if label == 'locpot_z_profile()':
        from siman.chg.chg_func import locpot_z_profile
        coord, values = locpot_z_profile(
            cl, filename=str(folder / '1.LOCPOT'), plot=False, return_numpy=True)
        return f'points={len(coord)} last={float(values[-1]):.4f}'
    if label == 'bader()':
        import shutil
        from siman import header
        binary = getattr(header, 'PATH2BADER', None)
        if not binary or not Path(str(binary)).is_file():
            binary = project_setting('PATH2BADER')
        if not binary or not Path(str(binary)).is_file():
            binary = shutil.which('bader')
        if not binary:
            return ('SKIP', 'the bader executable is not installed')
        header.PATH2BADER = binary
        cl.determine_filenames()
        cl.bader()
        acf = Path(cl.dir) / 'ACF.dat'
        if not acf.is_file() or acf.stat().st_size == 0:
            raise AssertionError('ACF.dat was not written')
        charges = []
        for line in acf.read_text(errors='replace').splitlines():
            parts = line.split()
            if parts and parts[0].isdigit() and len(parts) >= 5:
                charges.append(float(parts[4]))
        if len(charges) != cl.end.natom:
            raise AssertionError(f'{len(charges)} Bader charges for {cl.end.natom} atoms')
        return f'charges={len(charges)} sum={sum(charges):.3f}'
    if label == 'get_occ_mat()':
        atoms = list(getattr(cl, 'occ_matrices', {}) or {})
        if not atoms:
            raise AssertionError('OUTCAR has no occupation matrices')
        matrix = cl.get_occ_mat(atoms[0])
        rows = list(matrix)
        flat = [float(value) for row in rows for value in row]
        if len(rows) != 10 or not flat or not all(math.isfinite(value) for value in flat):
            raise AssertionError(f'occupation matrix has {len(rows)} rows')
        return f'atom={atoms[0]} rows={len(rows)}'
    if label == 'msd()':
        try:
            import mlyzed  # noqa: F401
        except ImportError:
            return ('SKIP', 'mlyzed is not installed')
        frames = xdatcar_frames(folder / '1.XDATCAR')
        cl.path['xdatcar'] = str(folder / '1.XDATCAR')
        cl.msd(el='Li')
        return f'frames={frames}'
    if label == 'plot_cp_convergence()':
        return str(cl.plot_cp_convergence())
    raise RuntimeError(f'unknown VASP method {label}')


def _li_potcar():
    """POTCAR from the archived job. A fresh download keeps only OUTCAR and CONTCAR."""
    for path in (
        HERE / 'Li.sp' / 'POTCAR',
        HERE.parent / 'siman_bench' / 'Li.sp' / 'POTCAR',
    ):
        if path.is_file():
            return path
    return None


def call_local(label):
    """Exercise VASP writers on the primitive Li cell. Nothing is submitted."""
    import os
    import tempfile

    import numpy as np

    from siman import header
    from siman.calc_manage import add_des, smart_structure_read
    from siman.calculators.vasp import CalculationVasp
    from siman.set_functions import InputSet

    potcar = _li_potcar()
    poscar = HERE.parent / 'poscar' / 'Li_mp-135.POSCAR-1'
    if not poscar.is_file():
        raise FileNotFoundError(poscar)
    structure = smart_structure_read(str(poscar))
    structure.poscar_atom_order = list(range(structure.natom))
    folder = tempfile.mkdtemp(prefix='siman_vasp_input_')
    inset = InputSet()
    inset.ise = 'bench'
    inset.vasp_params = {
        'NSW': 0, 'IBRION': -1, 'ENCUT': 400, 'ISMEAR': 1, 'SIGMA': 0.1,
        'PREC': 'Normal', 'EDIFF': 1e-5, 'KSPACING': 0.5,
    }
    if potcar is not None:
        inset.path_to_potcar = str(potcar)
    else:
        inset.potdir = {3: 'Li'}
        library = project_setting('PATH2POTENTIALS')
        if library:
            inset.path2pot = library
    cl = CalculationVasp(inset, iid=('LiBench', 'unit', 1))
    cl.init = structure
    cl.end = structure.copy()
    cl.des = 'benchmark'
    cl.dir = folder + os.sep
    os.makedirs(cl.dir, exist_ok=True)
    add_des(header.struct_des, 'LiBench', 'siman_bench', des='benchmark', override=True)

    if label == 'make_incar()':
        written = cl.make_incar()
        text = open(written[0], encoding='utf-8').read()
        if 'ENCUT' not in text:
            raise AssertionError(text)
        return written[0]
    if label == 'make_kpoints_file()':
        written = cl.make_kpoints_file()
        return str(written)
    if label == 'update_incar()':
        cl.make_incar()
        here = os.getcwd()
        try:
            os.chdir(cl.dir)
            cl.update_incar('ENCUT', 450, write=True, run=True)
        finally:
            os.chdir(here)
        text = open(os.path.join(cl.dir, 'INCAR'), encoding='utf-8').read()
        if '450' not in text:
            raise AssertionError(text)
        return 'ENCUT=450'
    if label == 'add_potcar()':
        path = cl.add_potcar()
        if not path or not os.path.isfile(path):
            raise AssertionError(str(path))
        return path
    if label == 'calculate_nbands()':
        source = str(potcar) if potcar is not None else cl.add_potcar()
        cl.calculate_nbands(cl.set, path_to_potcar=source)
        electrons = sum(cl.init.zval)
        if electrons <= 0:
            raise AssertionError(cl.init.zval)
        return f'valence electrons={electrons}'
    if label == 'check_kpoints()':
        mesh = cl.check_kpoints([4, 4, 4])
        if list(mesh) != [4, 4, 4]:
            raise AssertionError(mesh)
        return str(mesh)
    matrix = np.eye(5).tolist() + np.eye(5).tolist()
    cl.occ_matrices = {}
    if label == 'set_occ_mat()':
        stored = cl.set_occ_mat(0, matrix)
        if len(stored.occ_matrices[0]) != 10:
            raise AssertionError(len(stored.occ_matrices[0]))
        return '10 rows'
    if label == 'get_occ_mat()':
        stored = cl.set_occ_mat(0, matrix)
        read = stored.get_occ_mat(0)
        if len(read) != 10:
            raise AssertionError(len(read))
        return '10 rows'
    if label == 'write_occmatrix()':
        cl.occ_matrices = {0: matrix}
        from siman.inout import write_occmatrix
        path = write_occmatrix(cl.occ_matrices, folder)
        text = open(path, encoding='utf-8').read()
        if 'spin 1' not in text or 'spin 2' not in text:
            raise AssertionError(text)
        return path
    if label == 'vasp_dipole_center()':
        center = cl.vasp_dipole_center()
        return f'center={center}'
    if label == 'determine_filenames()':
        cl.path['output'] = cl.dir + '1.OUTCAR'
        cl.determine_filenames()
        if not str(cl.path.get('chgcar', '')).endswith('1.CHGCAR'):
            raise AssertionError(cl.path.get('chgcar'))
        return cl.path['chgcar']
    if label in ('visualize_grid() vesta', 'visualize_grid() jmol'):
        from viewer_paths import apply_viewer_paths, close_new_viewers, viewer_pids
        paths = apply_viewer_paths()
        program = 'vesta' if label.endswith('vesta') else 'jmol'
        key = 'PATH2VESTA' if program == 'vesta' else 'PATH2JMOL'
        grid_name, grid_dir = _grid_for_view(cl.dir)
        cl.dir = grid_dir
        cl.path['output'] = os.path.join(grid_dir, '1.OUTCAR')
        before = viewer_pids()
        cube = cl.visualize_grid(grid_name, program=program)
        closed = close_new_viewers(before)
        if not os.path.isfile(cube):
            raise AssertionError(cube)
        header_line = open(cube, encoding='utf-8').readline()
        if not header_line.startswith('JVXL'):
            raise AssertionError('cube header is not readable by Jmol: ' + header_line.strip())
        if header.PATH2VESTA != paths['PATH2VESTA'] or header.PATH2JMOL != paths['PATH2JMOL']:
            raise AssertionError('viewer path was replaced')
        if not closed:
            raise AssertionError(f'{program} did not start')
        return f'{key}={paths[key]} closed'
    raise RuntimeError(f'unknown local VASP method {label}')


def _grid_for_view(fallback_dir):
    """Use a finished VASP grid when one is present, otherwise a two-point LOCPOT."""
    import os
    for folder, name in (
        (HERE / 'Li_slab.sol', 'LOCPOT'),
        (HERE / 'Li_slab.sol', 'PHI'),
        (HERE / 'Li2.bader', 'CHGCAR'),
        (HERE / 'Li.sp', 'CHGCAR'),
    ):
        if (folder / f'1.{name}').is_file() or (folder / f'1.{name}.gz').is_file():
            return name, str(folder) + '/'
    path = os.path.join(fallback_dir, '1.LOCPOT')
    with open(path, 'w', encoding='utf-8') as locpot:
        locpot.write(
            'tiny\n1.0\n'
            '  2.0 0.0 0.0\n  0.0 2.0 0.0\n  0.0 0.0 2.0\n'
            'Li\n1\nDirect\n  0.0 0.0 0.0\n\n'
            '2 2 2\n 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8\n'
        )
    return 'LOCPOT', fallback_dir
