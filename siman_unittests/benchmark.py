"""
Run the SIMAN benchmarks.

    python3 benchmark.py --path=/path/to/siman --level=all

--level=default runs structure methods only.
--level=methods adds the VASP calculator methods, including DOS and band plots.
--level=all runs every calculator. VASP is the only calculator in this suite.
The report is written to report.log in the project directory.
"""

import argparse
import sys
from datetime import datetime
from pathlib import Path


HERE = Path(__file__).resolve().parent

# default: benchmark_default.py, structure edits only.
# methods: that file plus benchmark_methods_vasp.py (DOS, bands, OUTCAR, OSZICAR).
# all: every calculator. VASP is the only one registered here.
LEVELS = {
    'default': ('default',),
    'methods': ('default', 'vasp'),
    'all': ('default', 'vasp'),
}


def siman_on_path(path):
    """Put the requested SIMAN tree ahead of every other copy."""
    path = Path(path).expanduser().resolve()
    if (path / '__init__.py').is_file() and path.name == 'siman':
        root = path.parent
    elif (path / 'siman' / '__init__.py').is_file():
        root = path
    else:
        raise SystemExit(
            f'{path} is not a siman package and does not contain one')
    sys.path.insert(0, str(root))
    for name in list(sys.modules):
        if name == 'siman' or name.startswith('siman.'):
            del sys.modules[name]
    return root


class Report:
    def __init__(self, total):
        self.total = total
        self.done = 0
        self.header = []
        self.benchmarks = []
        self.current = None

    def section(self, title):
        self.current = {
            'index': len(self.benchmarks) + 1,
            'name': title,
            'tests': [],
        }
        self.benchmarks.append(self.current)

    def add(self, status, name, detail=''):
        if self.current is None:
            self.section('Tests')
        number = len(self.current['tests']) + 1
        record = {
            'number': number,
            'name': name,
            'status': status,
            'detail': detail,
        }
        self.current['tests'].append(record)
        self.done += 1
        self._progress(record)

    def _progress(self, record):
        total = max(self.total, self.done)
        width = 24
        filled = int(width * self.done / total) if total else width
        bar = '#' * filled + '-' * (width - filled)
        bench = f'{self.current["index"]}. {self.current["name"]}'
        mark = '!!! FAIL !!!' if record['status'] == 'FAIL' else record['status']
        print(
            f'[{self.done:3d}/{total}] [{bar}] {bench} | '
            f'{record["number"]}. {record["name"]} | {mark}',
            flush=True,
        )

    def text(self):
        lines = list(self.header)
        for bench in self.benchmarks:
            lines.append('')
            lines.append(f'Benchmark {bench["index"]}. {bench["name"]}')
            lines.append('-' * len(lines[-1]))
            for test in bench['tests']:
                detail = f'  {test["detail"]}' if test['detail'] else ''
                if test['status'] == 'FAIL':
                    lines.append(
                        f'!!! {test["number"]}. {test["name"]}  FAIL{detail} !!!'
                    )
                else:
                    lines.append(
                        f'{test["number"]}. {test["name"]}  {test["status"]}{detail}'
                    )

        lines.append('')
        lines.append('Compiled')
        lines.append('--------')
        total_tests = sum(len(bench['tests']) for bench in self.benchmarks)
        successful = self._count('PASS')
        failed = self._count('FAIL')
        skipped = self._count('SKIP')
        lines.append(f'Total tests: {total_tests}')
        lines.append(f'Successful: {successful}')
        lines.append(f'Failed: {failed}')
        lines.append(f'Skipped: {skipped}')
        for bench in self.benchmarks:
            tests = bench['tests']
            passed = sum(1 for test in tests if test['status'] == 'PASS')
            broken = sum(1 for test in tests if test['status'] == 'FAIL')
            left = sum(1 for test in tests if test['status'] == 'SKIP')
            grouped = {}
            for test in tests:
                method = test['name']
                call = method.find('()')
                if call != -1:
                    method = method[:call + 2]
                grouped.setdefault(method, []).append(test['status'])
            successful_methods = [
                method for method, statuses in grouped.items() if 'FAIL' not in statuses and 'PASS' in statuses
            ]
            failed_methods = [
                method for method, statuses in grouped.items() if 'FAIL' in statuses
            ]
            lines.append('')
            lines.append(f'Benchmark {bench["index"]}. {bench["name"]}')
            lines.append(f'  tests: {len(tests)}')
            lines.append(f'  successful: {passed}')
            lines.append(f'  failed: {broken}')
            lines.append(f'  skipped: {left}')
            lines.append('  methods: ' + ', '.join(grouped))
            lines.append('  successful methods: ' + (', '.join(successful_methods) or 'none'))
            lines.append('  failed methods: ' + (', '.join(failed_methods) or 'none'))

        lines.append('')
        lines.append('Failed')
        lines.append('------')
        failures = [
            (bench, test)
            for bench in self.benchmarks
            for test in bench['tests']
            if test['status'] == 'FAIL'
        ]
        if not failures:
            lines.append('No failed tests.')
        for bench, test in failures:
            lines.append(
                f'!!! Benchmark {bench["index"]}. {bench["name"]} | '
                f'{test["number"]}. {test["name"]} !!!'
            )
            lines.append(f'    reason: {test["detail"] or "failed"}')
        return '\n'.join(lines) + '\n'

    def _count(self, status):
        return sum(
            1
            for bench in self.benchmarks
            for test in bench['tests']
            if test['status'] == status
        )


def main():
    parser = argparse.ArgumentParser(description='Run the SIMAN benchmarks')
    parser.add_argument('--path', required=True, help='path to the siman package or the directory that contains it')
    parser.add_argument('--level', required=True, choices=sorted(LEVELS), help='all, methods, or default')
    args = parser.parse_args()
    selected = LEVELS[args.level]
    root = siman_on_path(args.path)
    import siman
    total = 0
    if 'default' in selected:
        from benchmark_default import planned as planned_default
        total += planned_default()
    if 'vasp' in selected:
        from benchmark_methods_vasp import planned as planned_vasp
        total += planned_vasp()
    report = Report(total)
    report.header.extend([
        f'SIMAN file: {siman.__file__}',
        f'SIMAN path: {root}',
        f'Time: {datetime.now().isoformat(timespec="seconds")}',
        f'Level: {args.level}',
    ])
    if 'default' in selected:
        from benchmark_default import run as run_default
        run_default(report)
    if 'vasp' in selected:
        from benchmark_methods_vasp import run as run_vasp
        run_vasp(report)
    text = report.text()
    destination = HERE.parent / 'report.log'
    destination.write_text(text)
    print()
    print(text, end='')
    print(f'Report written to {destination}')
    return 1 if report._count('FAIL') else 0


if __name__ == '__main__':
    sys.exit(main())
