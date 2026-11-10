# SIMAN unittests

These tests check a SIMAN installation against a small set of finished VASP calculations: bulk Li (static, DOS, bands, and a 10-point equation of state), a two-atom Li cell (relaxation, a 10-step AIMD run, Bader charges, and phonons), MnO (polaron, DFT+U occupations, and noncollinear moments), and a Li slab with solvent and constant-potential settings.

The saved outputs live in `siman_unittests/`. The tests read those files. They do not submit calculations unless you turn that on in `siman_unittests.py`.

## Run

Install SIMAN, then from this directory:

```bash
python3 siman_unittests/benchmark.py --path=/path/to/siman --level=default
python3 siman_unittests/benchmark.py --path=/path/to/siman --level=methods
```

`--level=default` checks structure methods on `poscar/Li_mp-135.POSCAR-1` and `poscar/MnO_mp-19006.POSCAR-1`. `--level=methods` and `--level=all` also parse the saved VASP outputs and call the calculator methods. The report is written to `siman_unittests/report.log`.

## Configure

Copy `project_conf.example.py` to `project_conf.py` and set the paths on this computer:

- `PATH2POTENTIALS` — VASP POTCAR library. POTCAR files are not in this repository. The input-writer tests build a Li POTCAR from this directory when `siman_unittests/Li.sp/POTCAR` is absent.
- `PATH2BADER` — Henkelman `bader` executable used by the Bader test.
- `PATH2JMOL` and `PATH2VESTA` — viewer programs. Those two tests open a window and close it.
- `CLUSTERS` — only needed to write or submit jobs. `siman_unittests.py` uses the cluster names `razor128` and `razor128_cp`.

`project_conf.py` is gitignored. Leave API keys and home-directory paths there.

## Replace the saved calculations

In `siman_unittests.py`, set one flag and run `python3 siman_unittests.py`:

- `dry_run` writes the inputs under `siman_bench/` and does not submit.
- `run` submits every benchmark job again.
- `read` downloads and parses the outputs.
- `copy` copies finished folders into `siman_unittests/`.

The constant-potential job is submitted only after the solvent job has a Fermi level.
