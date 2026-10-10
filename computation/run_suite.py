"""Run checks sequentially, optionally regenerating all n=5 evidence.

The regeneration mode runs in a temporary directory and compares exact
certificate data. It never overwrites the supplied certificates. RSS values
in the receipt use the operating system's child-process high-water mark.
"""
import argparse
import hashlib
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LIMIT = 1_500_000_000


def run(script, args=(), base=ROOT):
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', PYTHONHASHSEED='0')
    command = [sys.executable, '-u', str(base/'computation'/script), *args]
    label = script.replace('.py', '') + ('-' + '-'.join(args) if args else '')
    log = ROOT/'results'/f'{label}.log'
    started = time.monotonic()
    print('Running', script, *args, flush=True)
    with log.open('w') as output:
        result = subprocess.run(command, cwd=base, env=env,
                                stdout=output, stderr=subprocess.STDOUT)
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if platform.system() != 'Darwin':
        rss *= 1024
    receipt = dict(script=script, arguments=list(args), returncode=result.returncode,
                   seconds=round(time.monotonic()-started, 3),
                   cumulative_child_peak_rss_bytes=rss)
    print(receipt, flush=True)
    if result.returncode:
        print(log.read_text(), file=sys.stderr)
        raise RuntimeError(f'{script} failed; see results/{log.name}')
    if rss >= LIMIT:
        raise RuntimeError('Child peak RSS exceeds the 1.5 GB budget.')
    return receipt


def without_timing(value):
    if isinstance(value, list):
        return [without_timing(v) for v in value]
    if isinstance(value, dict):
        return {k: without_timing(v) for k, v in value.items() if k != 'seconds'}
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--regenerate', action='store_true')
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('Run Python without -O: the checkers use assertions.')
    (ROOT/'results').mkdir(exist_ok=True)
    receipts = []
    if args.regenerate:
        with tempfile.TemporaryDirectory(prefix='negami-reproduction-') as temp:
            base = Path(temp)
            shutil.copytree(ROOT/'computation', base/'computation',
                            ignore=shutil.ignore_patterns('__pycache__'))
            (base/'artifacts').mkdir()
            (base/'results').mkdir()
            for script, extra in [
                ('compute.py', ('--mode', 'main')),
                ('exact_check.py', ()),
                ('diagonal_symbolic.py', ()),
                ('lift_certificate.py', ()),
                ('sample_torus.py', ()),
                ('verify_certificates.py', ()),
                ('verify_n5.py', ())
            ]:
                receipts.append(run(script, extra, base))
            compared = []
            for source in sorted((ROOT/'artifacts').glob('*.json')):
                if source.name == 'jordan-laurent.json':
                    continue
                generated = base/'artifacts'/source.name
                assert without_timing(json.loads(source.read_text())) == \
                       without_timing(json.loads(generated.read_text())), source.name
                compared.append(source.name)
            for source in (base/'results').glob('*'):
                shutil.copy2(source, ROOT/'results'/source.name)
            for name in ('rational-checks.json',):
                shutil.copy2(base/'artifacts'/name, ROOT/'results'/name)
            print('Regenerated certificates match:', ', '.join(compared), flush=True)
    else:
        receipts.append(run('verify_n5.py'))
    receipts.append(run('verify_jordan.py'))
    report = dict(python=platform.python_version(), platform=platform.system(),
                  mode='regenerate' if args.regenerate else 'verify',
                  single_computation_process=True, rss_limit_bytes=LIMIT,
                  commands=receipts,
                  certificate_sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                                      for f in sorted((ROOT/'artifacts').glob('*.json'))})
    filename = 'regeneration-run.json' if args.regenerate else 'verification-run.json'
    (ROOT/'results'/filename).write_text(json.dumps(report, indent=2)+'\n')
    print('ALL CHECKS PASSED', flush=True)


if __name__ == '__main__':
    main()
