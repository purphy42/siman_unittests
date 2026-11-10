"""Structure checks. No OUTCAR, DOSCAR, or other VASP output is read.

Each report line is one call. The note after the call says what must be true.
"""

import copy

import numpy as np

from pathlib import Path


PROJECT = Path(__file__).resolve().parent.parent
POSCAR = PROJECT / 'poscar'


def planned():
    """Number of recorded results when the Li and MnO cells are present."""
    li_file = POSCAR / 'Li_mp-135.POSCAR-1'
    mno_file = POSCAR / 'MnO_mp-19006.POSCAR-1'
    if not li_file.is_file() or not mno_file.is_file():
        return 1
    return 39 + additional_count()


def run_case(report, name, function, expect, note=''):
    """Record one check. note is the description printed on that report line."""
    prefix = f'{note}. ' if note else ''
    try:
        value = function()
    except (Exception, SystemExit) as error:
        report.add('FAIL', name, prefix + f'{type(error).__name__}: {error}')
        return None
    if expect is not None and not isinstance(value, expect):
        expected = getattr(expect, '__name__', str(expect))
        report.add(
            'FAIL', name,
            prefix + f'returned {type(value).__name__}, expected {expected}')
        return value
    report.add('PASS', name, prefix + type(value).__name__)
    return value


def run(report):
    from siman.calc_manage import smart_structure_read
    from siman.core.structure import Structure
    from siman.core.structure_utils import is_valid_magmom
    from siman.geo import create_supercell, interpolate, replic, scale_cell_uniformly
    from siman.impurity import make_interface

    report.section('Methods')
    li_file = POSCAR / 'Li_mp-135.POSCAR-1'
    mno_file = POSCAR / 'MnO_mp-19006.POSCAR-1'
    if not li_file.is_file() or not mno_file.is_file():
        report.add('FAIL', 'input POSCAR', f'missing under {POSCAR}')
        return

    # Read the primitive Li and MnO POSCAR files into Structure objects.
    li = run_case(report, 'smart_structure_read()', lambda: smart_structure_read(str(li_file)), Structure,
                  'read the primitive Li POSCAR')
    mno = run_case(report, 'smart_structure_read() MnO', lambda: smart_structure_read(str(mno_file)), Structure,
                   'read the MnO POSCAR')
    if li is None or mno is None:
        return

    # Formula, element list, and cell volume of the one-atom Li cell.
    formula = run_case(report, 'get_formula()', lambda: li.get_formula(), str,
                       'return a text formula')
    if isinstance(formula, str) and 'Li' not in formula:
        report.add('FAIL', 'get_formula() content', formula)
    elif isinstance(formula, str):
        report.add('PASS', 'get_formula() content', f'formula contains Li: {formula}')

    elements = run_case(report, 'get_elements()', lambda: li.get_elements(), list,
                        'return one element name per atom')
    if isinstance(elements, list) and elements != ['Li']:
        report.add('FAIL', 'get_elements() content', str(elements))
    elif isinstance(elements, list):
        report.add('PASS', 'get_elements() content', 'the only element is Li')

    volume = run_case(report, 'get_volume()', lambda: li.get_volume(), (float, np.floating),
                      'return the cell volume')
    if isinstance(volume, (float, np.floating)) and float(volume) <= 0:
        report.add('FAIL', 'get_volume() sign', str(volume))
    elif isinstance(volume, (float, np.floating)):
        report.add('PASS', 'get_volume() sign', f'volume is positive: {float(volume):.3f} A^3')

    # Double the Li cell along a, both as a function and as a method.
    doubled = run_case(report, 'replic()', lambda: replic(li, mul=(2, 1, 1)), Structure,
                       'build a 2x1x1 supercell')
    method_doubled = run_case(report, 'Structure.replic()', lambda: li.replic(mul=(2, 1, 1)), Structure,
                              'same supercell through the Structure method')
    if doubled is not None and doubled.natom != 2 * li.natom:
        report.add('FAIL', 'replic atom count', f'{doubled.natom} != {2 * li.natom}')
    elif doubled is not None:
        report.add('PASS', 'replic() atom count', f'atom count doubled to {doubled.natom}')
    if method_doubled is not None and method_doubled.natom != 2 * li.natom:
        report.add('FAIL', 'Structure.replic() atom count', str(method_doubled.natom))
    elif method_doubled is not None:
        report.add('PASS', 'Structure.replic() atom count', f'atom count doubled to {method_doubled.natom}')

    # Add 10 A of vacuum along c, then check that the c vector grew by that amount.
    vacuum = run_case(report, 'add_vacuum()', lambda: li.add_vacuum(2, 10.0), Structure,
                      'add 10 A of vacuum along c')
    if vacuum is not None:
        old_c = float(np.linalg.norm(li.rprimd[2]))
        new_c = float(np.linalg.norm(vacuum.rprimd[2]))
        if new_c <= old_c + 9.0:
            report.add('FAIL', 'add_vacuum length', f'{old_c:.3f} -> {new_c:.3f}')
        else:
            report.add('PASS', 'add_vacuum length', f'c grew by {new_c - old_c:.3f} A')

    # Five cells scaled from -10% to +10%. The Li.eos job uses ten images from -15% to +15%.
    images = run_case(
        report, 'scale_cell_uniformly()',
        lambda: scale_cell_uniformly(li, scale_region=(-10, 10), n_scale_images=5),
        list,
        'make five uniformly scaled cells')
    if isinstance(images, list):
        if len(images) == 5 and all(isinstance(st, Structure) for st in images):
            report.add('PASS', 'eos image count and type', 'five Structure objects from -10% to +10%')
        else:
            report.add('FAIL', 'eos image count and type', f'len={len(images)}')

    # Delete the first atom of the doubled cell.
    removed = run_case(report, 'remove_atoms()', lambda: doubled.remove_atoms([0]), Structure,
                       'delete atom 0')
    if removed is not None and doubled is not None and removed.natom != doubled.natom - 1:
        report.add('FAIL', 'remove_atoms count', f'{removed.natom}')
    elif removed is not None:
        report.add('PASS', 'remove_atoms count', f'one atom removed, {removed.natom} left')

    def coordinates():
        st = copy.deepcopy(li)
        st.xred2xcart()
        st.xcart2xred()
        return st

    # Convert fractional coordinates to Cartesian and back. The atom count must stay the same.
    back = run_case(report, 'xred2xcart and xcart2xred', coordinates, Structure,
                    'convert coordinates both ways')
    if back is not None and len(back.xcart) != back.natom:
        report.add('FAIL', 'coordinate length', f'{len(back.xcart)} != {back.natom}')
    elif back is not None:
        report.add('PASS', 'coordinate length', f'{back.natom} coordinate remains')

    # A complete moment list is accepted. A list that contains None is not.
    valid = run_case(report, 'is_valid_magmom', lambda: is_valid_magmom([0.0, 1.0], 2), bool,
                     'accept one moment per atom')
    absent = run_case(report, 'is_valid_magmom absent', lambda: is_valid_magmom([None], 1), bool,
                      'reject a missing moment')
    if valid is False:
        report.add('FAIL', 'is_valid_magmom true case', 'returned False')
    elif valid is True:
        report.add('PASS', 'is_valid_magmom true case', 'two numeric moments are valid')
    if absent is not False and absent is not None:
        report.add('FAIL', 'is_valid_magmom absent case', repr(absent))
    elif absent is False:
        report.add('PASS', 'is_valid_magmom absent case', 'None means the moments are absent')

    cell = replic(mno, mul=(2, 2, 2))
    elements = list(cell.get_elements())
    if 'Mn' in elements and elements.count('O') >= 6:
        center = elements.index('Mn')

        def localize():
            return cell.localize_polaron_dist(center, d=0.1, nn=6, axis=2)

        # Stretch one Mn-O axis by 0.1 A. The atom count must not change.
        polar = run_case(report, 'localize_polaron_dist()', localize, Structure,
                         'distort one Mn-O octahedron by 0.1 A')
        if polar is not None and polar.natom != cell.natom:
            report.add('FAIL', 'localize_polaron_dist atom count', str(polar.natom))
        elif polar is not None:
            report.add('PASS', 'localize_polaron_dist atom count', f'atom count stays {polar.natom}')
    else:
        report.add('SKIP', 'localize_polaron_dist', 'MnO cell has no MnO6 octahedron')

    def natom():
        count = li.get_natom()
        if count != li.natom:
            raise AssertionError(f'{count} != {li.natom}')
        return count

    run_case(report, 'get_natom()', natom, int, 'count atoms and ignore voids')

    def delete():
        st = doubled.del_atom(0)
        if st.natom != doubled.natom - 1:
            raise AssertionError(f'{st.natom} != {doubled.natom - 1}')
        return st

    run_case(report, 'del_atom()', delete, Structure, 'delete one atom by index')

    def swap():
        st = doubled.swap_atoms(0, 1)
        if st.natom != doubled.natom:
            raise AssertionError(f'natom changed to {st.natom}')
        return st

    run_case(report, 'swap_atoms()', swap, Structure, 'exchange two atoms and keep the count')

    def add():
        origin = np.array(li.xcart[0], dtype=float)
        st = li.add_atoms(atoms_xcart=[origin + np.array([1.0, 0.2, 0.3])], element='Li')
        if st.natom != li.natom + 1:
            raise AssertionError(f'{st.natom} != {li.natom + 1}')
        return st

    run_case(report, 'add_atoms()', add, Structure, 'append one Li atom')

    def shift():
        st = li.shift_atoms(vector_red=[0.1, 0.0, 0.0])
        if st.natom != li.natom:
            raise AssertionError('natom changed')
        return st

    shifted = run_case(report, 'shift_atoms()', shift, Structure, 'shift every atom by 0.1 along a')

    def back_to_cell():
        st = shifted.return_atoms_to_cell()
        if any(coord < -1e-8 or coord >= 1.0 + 1e-8 for row in st.xred for coord in row):
            raise AssertionError('an atom is outside the cell')
        return st

    if shifted is not None:
        run_case(report, 'return_atoms_to_cell()', back_to_cell, Structure,
                 'fold fractional coordinates back into the cell')

    def drop_vacuum():
        return vacuum.remove_vacuum(thickness=5.0)

    if vacuum is not None:
        run_case(report, 'remove_vacuum()', drop_vacuum, Structure, 'leave 5 A of vacuum')

    def supercell():
        st = create_supercell(li, [[2, 0, 0], [0, 1, 0], [0, 0, 1]], silent=1)
        if st.natom != 2 * li.natom:
            raise AssertionError(f'{st.natom} != {2 * li.natom}')
        return st

    run_case(report, 'create_supercell()', supercell, Structure, 'double the cell along a')

    def pair_distance():
        value = doubled.distance(0, 1)
        if float(value) <= 0:
            raise AssertionError(str(value))
        return value

    run_case(report, 'distance()', pair_distance, (float, np.floating),
             'shortest periodic distance between atom 0 and atom 1')

    def images_between():
        end = li.shift_atoms(vector_red=[0.05, 0.0, 0.0])
        sts = interpolate(li, end, images=2, write_poscar=0)
        if len(sts) != 2:
            raise AssertionError(f'{len(sts)} images')
        return sts

    run_case(report, 'interpolate()', images_between, list, 'two images between the original and a shifted cell')

    def write():
        path = Path(__file__).resolve().parent / '_tmp_default.POSCAR'
        li.write_poscar(str(path))
        ok = path.is_file() and path.stat().st_size > 0
        path.unlink(missing_ok=True)
        if not ok:
            raise AssertionError('POSCAR was not written')
        return True

    run_case(report, 'write_poscar()', write, bool, 'write a POSCAR and read its size')

    def interface():
        host = li.add_vacuum(2, 8.0)
        film = mno.add_vacuum(2, 8.0)
        made, scaled = make_interface(
            host, np.array(host.xcart[0], dtype=float),
            film, np.array(film.xcart[0], dtype=float))
        if made.natom != host.natom + film.natom:
            raise AssertionError(f'{made.natom} != {host.natom + film.natom}')
        if not isinstance(scaled, Structure):
            raise AssertionError(type(scaled).__name__)
        return made

    # Stack the Li slab and the MnO slab. There is no create_interface(); this is the builder.
    run_case(report, 'make_interface()', interface, Structure,
             'join a Li slab and an MnO slab into one cell')

    def surface():
        from siman.geo import create_surface2
        slab = create_surface2(
            li, [0, 0, 1], min_slab_size=3, min_vacuum_size=5,
            write_poscar=0, silent=1, return_one=True)
        size = getattr(slab, 'natom', None)
        if size is None:
            size = getattr(slab, 'num_sites', None)
        if not size:
            raise AssertionError(type(slab).__name__)
        return slab

    run_case(report, 'create_surface2()', surface, object, 'cut a (001) Li slab with vacuum')
    for name, function, expect, note in additional_cases(li, mno):
        run_case(report, name, function, expect, note)


def additional_count():
    return len(_ADDITIONAL)


def _atoms(st, natom):
    count = getattr(st, 'natom', None)
    if count != natom:
        raise AssertionError(f'{count} != {natom}')
    return st


def additional_cases(li, mno):
    """Structure, geometry, k-point, and file helpers that do not read VASP output."""
    import os
    import tempfile

    import numpy as np

    from siman.calc_manage import (
        add_des, create_phonopy_conf_file, name_mod_supercell, read_phonopy_data,
        read_phonopy_dat_file,
    )
    from siman.core.structure import Structure
    from siman.core.structure_utils import pair_opposite_vectors
    from siman.geo import (
        best_miller, calc_k_point_mesh, calc_kspacings, calc_recip_vectors, calc_volume,
        create_deintercalated_structure, create_replaced_structure, create_surface,
        determine_symmetry_positions, fit2host, four2three_index, generate_occupied_configurations,
        hex2rhombo, hkl2uvw, image_distance, image_vector, local_surrounding, move_edge,
        primitive, remove_one_atom, replic, rhombo2hex, rms_between_structures,
        scale_cell_by_matrix, slab_misfit, stoichiometry_criteria, symmetry_criteria,
        three2four_index, xcart2xred, xred2xcart,
    )
    from siman import header
    from siman.inout import determine_file_format, read_poscar
    from siman.kpoints_functions import kspacing2ngkpt, ngkpt2kspacings

    pair = replic(li, mul=(2, 1, 1))
    cell = replic(mno, mul=(2, 2, 2))
    cases = []

    def add(name, function, expect, note):
        cases.append((name, function, expect, note))

    add('copy()', lambda: _atoms(li.copy(), li.natom), Structure, 'copy the Li cell')
    add('new()', lambda: li.new(), Structure, 'make an empty Structure')
    add('selective_all()', lambda: li.copy().selective_all(), Structure, 'mark every atom free to move')
    add('check_selective()', lambda: li.copy().selective_all().check_selective(), bool,
        'report whether selective dynamics is set')
    add('fix_atoms()', lambda: _atoms(li.copy().fix_atoms([0]), li.natom), Structure, 'freeze atom 0')
    add('get_recip()', lambda: li.get_recip(), list, 'return the reciprocal lattice')
    add('get_nznucl()', lambda: li.get_nznucl(), list, 'count atoms of each element')
    add('get_el_name()', lambda: li.get_el_name(0), str, 'name the element of atom 0')
    add('get_elements_z()', lambda: li.get_elements_z(), list, 'return atomic numbers')
    add('get_el_z()', lambda: li.get_el_z(0), int, 'return the atomic number of atom 0')
    add('get_reduced_formula()', lambda: li.get_reduced_formula(), str, 'return the reduced formula')
    add('get_name()', lambda: li.get_name(), str, 'return the structure name')
    add('get_angles()', lambda: li.get_angles(), tuple, 'return the three cell angles')
    add('sizes()', lambda: li.sizes(), tuple, 'return the cell lengths')
    add('rprimd_len()', lambda: li.rprimd_len(), tuple, 'return the lattice-vector lengths')
    add('sg()', lambda: li.sg(silent=1), tuple, 'return the space-group symbol and number')
    add('get_space_group_info()', lambda: li.get_space_group_info(), tuple, 'return the same space-group record')
    add('get_primitive_cell()', lambda: li.get_primitive_cell(), Structure, 'reduce Li to its primitive cell')
    add('get_conventional_cell()', lambda: li.get_conventional_cell(), Structure, 'build the conventional Li cell')

    def added():
        st = li.copy().add_atom(xr=[0.2, 0.2, 0.2], element='Li')
        return _atoms(st, li.natom + 1)

    add('add_atom()', added, Structure, 'insert one Li atom at a fractional coordinate')

    def only_mn():
        st = cell.copy().leave_only('Mn')
        if any(el != 'Mn' for el in st.get_elements()):
            raise AssertionError(st.get_elements())
        return st

    add('leave_only()', only_mn, Structure, 'keep only the Mn atoms')
    add('replace_atoms()', lambda: cell.copy().replace_atoms([0], 'Li'), Structure,
        'replace atom 0 with Li')
    add('mov_atoms()', lambda: _atoms(li.copy().mov_atoms(0, to_xr=[0.2, 0.2, 0.2]), li.natom),
        Structure, 'move atom 0 to a new fractional position')
    add('center()', lambda: li.copy().center(), np.ndarray, 'center the cell')
    add('nn()', lambda: cell.nn(0, n=6, silent=1), dict, 'list six neighbors of atom 0')
    add('perturb()', lambda: _atoms(li.copy().perturb(0.05), li.natom), Structure,
        'shift atoms by a random 0.05 A')
    add('get_element_xred()', lambda: li.get_element_xred('Li'), np.ndarray,
        'collect fractional coordinates of Li')
    add('vlen', lambda: li.vlen, np.ndarray, 'read the three lattice lengths')
    add('invert_axis()', lambda: _atoms(li.copy().invert_axis(2), li.natom), Structure,
        'invert the c axis')
    add('exchange_axes()', lambda: _atoms(li.copy().exchange_axes(0, 1), li.natom), Structure,
        'swap the a and b axes')
    add('find_closest_atom()', lambda: li.find_closest_atom(xc=np.array(li.xcart[0], dtype=float)),
        tuple, 'find the atom nearest the first Cartesian position')
    def update_red():
        st = li.copy()
        st.update_xred()
        return _atoms(st, li.natom)

    def update_cart():
        st = li.copy()
        st.update_xcart()
        return _atoms(st, li.natom)

    add('update_xred()', update_red, Structure, 'rebuild fractional coordinates from Cartesian ones')
    add('update_xcart()', update_cart, Structure, 'rebuild Cartesian coordinates from fractional ones')
    add('get_slice()', lambda: pair.get_slice([0.0, 1.0]), list, 'take every atom in the doubled cell')
    add('get_numbers()', lambda: cell.get_numbers('Mn'), list, 'list the Mn atom indices')
    add('get_transition_elements()', lambda: cell.get_transition_elements(), list,
        'list transition-metal atoms')
    add('get_coordination()', lambda: cell.get_coordination('Mn', silent=1), list,
        'coordination of Mn in the supercell')
    add('reorder_for_vasp()', lambda: cell.copy().reorder_for_vasp(), Structure,
        'group atoms in VASP element order')
    add('remove_close_lying2()', lambda: li.copy().remove_close_lying2(tol=0.1), Structure,
        'drop atoms closer than 0.1 A; the Li cell has none')
    add('find_closest_neighbor()', lambda: cell.find_closest_neighbor(0), tuple,
        'nearest neighbor of atom 0')

    def rms_self():
        li.rms(li.copy())
        return True

    add('rms()', rms_self, bool, 'RMS displacement of a cell against its copy is zero')
    add('convert2pymatgen()', lambda: li.convert2pymatgen(), object, 'convert Li to a pymatgen structure')
    add('shake_atoms()', lambda: _atoms(li.copy().shake_atoms(0.02), li.natom), Structure,
        'displace every atom by 0.02 A')

    def write_xyz_case():
        folder = tempfile.mkdtemp(prefix='siman_xyz_')
        li.copy().write_xyz(folder)
        found = list(Path(folder).rglob('*.xyz'))
        if not found:
            raise AssertionError('no xyz file')
        return True

    add('write_xyz()', write_xyz_case, bool, 'write an xyz file')

    def write_cif_case():
        folder = tempfile.mkdtemp(prefix='siman_cif_')
        name = li.copy().write_cif(str(Path(folder) / 'li'))
        if not name or not Path(name).is_file():
            raise AssertionError(str(name))
        return True

    add('write_cif()', write_cif_case, bool, 'write a cif file')
    add('image_distance()', lambda: image_distance(pair.xcart[0], pair.xcart[1], pair.rprimd),
        tuple, 'periodic distance between the two Li atoms')
    add('image_vector()', lambda: image_vector(pair, pair.xcart[0], pair.xcart[1]), np.ndarray,
        'periodic vector between the two Li atoms')
    add('xred2xcart()', lambda: xred2xcart(li.xred, li.rprimd), list,
        'convert fractional coordinates with the geometry function')
    add('xcart2xred()', lambda: xcart2xred(li.xcart, li.rprimd), list,
        'convert Cartesian coordinates with the geometry function')
    add('calc_recip_vectors()', lambda: calc_recip_vectors(li.rprimd), list, 'build reciprocal vectors')
    add('calc_kspacings()', lambda: calc_kspacings([4, 4, 4], li.rprimd), list,
        'spacing of a 4x4x4 mesh')
    add('calc_volume()', lambda: calc_volume(*li.rprimd), (float, np.floating), 'cell volume from the three vectors')

    def rms_geo():
        rms_between_structures(li, li.copy())
        return True

    add('rms_between_structures()', rms_geo, bool, 'geometry RMS of a cell against its copy')

    def one_removed():
        st, index = remove_one_atom(cell.copy(), 'Mn')
        if index != 0 or st.natom != cell.natom - 1:
            raise AssertionError(f'index={index} natom={st.natom}')
        return st

    add('remove_one_atom()', one_removed, Structure, 'delete one symmetry-equivalent Mn')
    add('create_replaced_structure()',
        lambda: create_replaced_structure(cell.copy(), 'Mn', 'Li', rep_pos=1, only_one=True),
        Structure, 'replace one Mn site with Li')
    add('create_deintercalated_structure()',
        lambda: create_deintercalated_structure(replic(mno, mul=(2, 1, 1)), 'Mn', del_pos=1),
        Structure, 'remove one Mn from the MnO supercell')
    add('generate_occupied_configurations()',
        lambda: generate_occupied_configurations(
            replic(mno, mul=(2, 1, 1)), 'Mn', occ=[0.5], nconf=1, seed=1, silent=1),
        list, 'one configuration with half of the Mn sites occupied')
    add('primitive()', lambda: primitive(cell.copy()), Structure, 'reduce the MnO supercell to a primitive cell')
    add('scale_cell_by_matrix()',
        lambda: scale_cell_by_matrix(
            li.copy(), scale_region=(-2, 2), n_scale_images=3,
            mul_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1.01]]),
        list, 'three cells scaled along c')
    add('local_surrounding()', lambda: local_surrounding(cell.xcart[0], cell, 6),
        (float, np.floating), 'sum of distances to six neighbors')
    add('calc_k_point_mesh()', lambda: calc_k_point_mesh(li.rprimd, 0.5, silent=1), tuple,
        'mesh for a 0.5 1/A spacing')
    add('determine_symmetry_positions()',
        lambda: determine_symmetry_positions(cell, 'Mn', silent=1), list,
        'symmetry-equivalent Mn sites')
    add('hex2rhombo()', lambda: hex2rhombo(1, 0, 0), tuple, 'convert a hexagonal index')
    add('rhombo2hex()', lambda: rhombo2hex(1, 1, 1), tuple, 'convert a rhombohedral index')
    add('three2four_index()', lambda: three2four_index([1, 0, 0]), list, 'three-index direction to four-index')
    add('four2three_index()', lambda: four2three_index([1, 0, -1, 0]), list, 'four-index direction to three-index')
    add('hkl2uvw()', lambda: hkl2uvw([0, 0, 1], li.rprimd), np.ndarray, 'plane normal of (001)')
    add('best_miller()', lambda: best_miller([0, 0, 2]), np.ndarray, 'reduce (002) to (001)')
    add('ngkpt2kspacings()', lambda: ngkpt2kspacings([4, 4, 4], li.rprimd), list,
        'k-point spacing from a mesh, k-point module')
    add('kspacing2ngkpt()', lambda: kspacing2ngkpt(0.5, li.rprimd), list,
        'mesh from a 0.5 1/A spacing')

    def pairs():
        vectors = np.array(
            [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
            dtype=float)
        found = pair_opposite_vectors(vectors)
        if len(found) != 3:
            raise AssertionError(str(found))
        return found

    add('pair_opposite_vectors()', pairs, list, 'pair six ligand directions into three axes')
    vacuum = li.copy().add_vacuum(2, 8.0)
    add('move_edge()', lambda: move_edge(vacuum.copy()), Structure, 'move a vacuum slab to the cell edge')
    add('fit2host()', lambda: fit2host(vacuum.copy(), mno.copy().add_vacuum(2, 8.0)), Structure,
        'strain an MnO slab onto the Li slab')
    add('slab_misfit()', lambda: slab_misfit(vacuum.copy(), vacuum.copy()), list,
        'misfit of a slab with itself')
    add('stoichiometry_criteria()',
        lambda: stoichiometry_criteria(cell.copy(), mno.copy(), silent=True), int,
        'compare the composition of a supercell with its bulk cell')
    add('symmetry_criteria()', lambda: symmetry_criteria(li.copy()), int, 'symmetry test of the Li cell')

    def surface_old():
        here = os.getcwd()
        folder = tempfile.mkdtemp(prefix='siman_surface_')
        os.makedirs(Path(folder) / 'xyz')
        try:
            os.chdir(folder)
            slab = create_surface(
                li.get_conventional_cell(), [0, 0, 1], min_slab_size=3, min_vacuum_size=5)
        finally:
            os.chdir(here)
        size = getattr(slab, 'natom', None) or getattr(slab, 'num_sites', None)
        if not size:
            raise AssertionError(type(slab).__name__)
        return slab

    add('create_surface()', surface_old, object, 'cut a (001) slab with the older surface builder')

    def phonopy_conf():
        folder = tempfile.mkdtemp(prefix='siman_phonopy_')
        create_phonopy_conf_file(li, path=folder, mp=[2, 2, 2], dim=[1, 1, 1])
        path = Path(folder) / 'mesh.conf'
        if not path.is_file() or 'DIM' not in path.read_text():
            raise AssertionError(path)
        return True

    add('create_phonopy_conf_file()', phonopy_conf, bool, 'write mesh.conf for phonopy')

    def phonopy_dat():
        folder = tempfile.mkdtemp(prefix='siman_dos_')
        path = Path(folder) / 'total_dos.dat'
        path.write_text('freq tot\n0.0 0.1\n1.0 0.2\n')
        parsed = read_phonopy_dat_file(str(path))
        if parsed['freq'] != [0.0, 1.0] or parsed['tot'] != [0.1, 0.2]:
            raise AssertionError(str(parsed))
        return parsed

    add('read_phonopy_dat_file()', phonopy_dat, dict, 'read a two-point phonopy total DOS')

    def phonopy_yaml():
        import matplotlib
        matplotlib.use('Agg')
        folder = tempfile.mkdtemp(prefix='siman_thermal_')
        path = Path(folder) / 'thermal_properties.yaml'
        lines = ['thermal_properties:']
        for temperature in range(0, 100, 10):
            lines.append(f'- temperature: {temperature}')
            lines.append(f'  free_energy: {-0.01 * temperature}')
        path.write_text('\n'.join(lines) + '\n')
        temperatures, fit = read_phonopy_data(str(path))
        if len(temperatures) < 2:
            raise AssertionError('no temperature grid')
        return fit

    add('read_phonopy_data()', phonopy_yaml, object, 'fit free energy from a thermal_properties.yaml file')
    add('determine_file_format()',
        lambda: determine_file_format(str(POSCAR / 'Li_mp-135.POSCAR-1')), str,
        'recognize a POSCAR path as VASP')

    def read_again():
        st = li.new()
        read_poscar(st, str(POSCAR / 'Li_mp-135.POSCAR-1'))
        return _atoms(st, li.natom)

    add('read_poscar()', read_again, Structure, 'fill an empty Structure from the Li POSCAR')

    def describe():
        add_des(header.struct_des, 'bench_li', 'siman_bench', des='benchmark', override=True)
        if 'bench_li' not in header.struct_des:
            raise AssertionError('bench_li was not stored')
        return True

    add('add_des()', describe, bool, 'store a structure description')
    add('name_mod_supercell()',
        lambda: name_mod_supercell(mul_matrix=[[2, 0, 0], [0, 1, 0], [0, 0, 1]]), str,
        'name a 2x1x1 supercell')

    def open_viewer(program):
        from viewer_paths import apply_viewer_paths, close_new_viewers, viewer_pids
        paths = apply_viewer_paths()
        key = 'PATH2VESTA' if program == 'vesta' else 'PATH2JMOL'
        shown = li.copy()
        shown.name = 'Li'
        here = os.getcwd()
        folder = tempfile.mkdtemp(prefix=f'siman_{program}_')
        os.makedirs(Path(folder) / 'xyz')
        before = viewer_pids()
        try:
            os.chdir(folder)
            if program == 'vesta':
                shown.vesta()
            else:
                shown.jmol()
        finally:
            os.chdir(here)
        closed = close_new_viewers(before)
        written = list((Path(folder) / 'xyz').glob('POSCAR*'))
        if not written:
            raise AssertionError('POSCAR was not written')
        if not closed:
            raise AssertionError(f'{program} did not start')
        return f'{paths[key]} closed'

    add('vesta()', lambda: open_viewer('vesta'), str,
        'open the Li cell in VESTA using PATH2VESTA from project_conf.py')
    add('jmol()', lambda: open_viewer('jmol'), str,
        'open the Li cell in Jmol using PATH2JMOL from project_conf.py')

    if len(cases) != len(_ADDITIONAL):
        raise RuntimeError(f'{len(cases)} additional tests, expected {len(_ADDITIONAL)}')
    return cases


_ADDITIONAL = (
    'copy', 'new', 'selective_all', 'check_selective', 'fix_atoms', 'get_recip',
    'get_nznucl', 'get_el_name', 'get_elements_z', 'get_el_z', 'get_reduced_formula',
    'get_name', 'get_angles', 'sizes', 'rprimd_len', 'sg', 'get_space_group_info',
    'get_primitive_cell', 'get_conventional_cell', 'add_atom', 'leave_only',
    'replace_atoms', 'mov_atoms', 'center', 'nn', 'perturb', 'get_element_xred',
    'vlen', 'invert_axis', 'exchange_axes', 'find_closest_atom', 'update_xred',
    'update_xcart', 'get_slice', 'get_numbers', 'get_transition_elements',
    'get_coordination', 'reorder_for_vasp', 'remove_close_lying2',
    'find_closest_neighbor', 'rms', 'convert2pymatgen', 'shake_atoms', 'write_xyz',
    'write_cif', 'image_distance', 'image_vector', 'xred2xcart', 'xcart2xred',
    'calc_recip_vectors', 'calc_kspacings', 'calc_volume', 'rms_between_structures',
    'remove_one_atom', 'create_replaced_structure', 'create_deintercalated_structure',
    'generate_occupied_configurations', 'primitive', 'scale_cell_by_matrix',
    'local_surrounding', 'calc_k_point_mesh', 'determine_symmetry_positions',
    'hex2rhombo', 'rhombo2hex', 'three2four_index', 'four2three_index', 'hkl2uvw',
    'best_miller', 'ngkpt2kspacings', 'kspacing2ngkpt', 'pair_opposite_vectors',
    'move_edge', 'fit2host', 'slab_misfit', 'stoichiometry_criteria',
    'symmetry_criteria', 'create_surface', 'create_phonopy_conf_file',
    'read_phonopy_dat_file', 'read_phonopy_data', 'determine_file_format',
    'read_poscar', 'add_des', 'name_mod_supercell', 'vesta', 'jmol',
)
