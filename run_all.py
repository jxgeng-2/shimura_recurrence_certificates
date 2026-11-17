"""One entrypoint: all certificates, with only the small under-200 sanity check."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from build_report import build, REPORTS, validate_reports, require

ROOT=Path(__file__).resolve().parent
STAGES=[('general_identities','certificates/general_identities.py'),
    ('families','certificates/families.py'),('notebooks','certificates/notebooks.py'),
    ('coefficients','certificates/coefficients.py'),('appendix_c','certificates/appendix_c.py'),
    ('nelson_inputs','certificates/nelson_inputs.py'),('numerical','certificates/numerical.py'),
    ('cm_claims','certificates/cm_claims.py'),('cm_sanity','certificates/cm_sanity.sage'),
    ('integration','certificates/integration.py')]

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def source_files():
    return [p for p in sorted(ROOT.rglob('*')) if p.is_file()
        and 'results' not in p.relative_to(ROOT).parts and '__pycache__' not in p.parts
        and not ({'.git', '.venv', 'venv', '.ipynb_checkpoints'} & set(p.relative_to(ROOT).parts))
        and not p.name.endswith(('.pyc','.sage.py')) and p.name not in ('REPORT.html',)]

def receipt(mode,stages):
    report_paths=set(REPORTS.values())|{'results/summary.json','results/recovered_parameters.json',
        'results/shimura_parameter_recovery.ipynb','results/shimura_schwarzian_14_15_21.ipynb',
        'results/coefficients_first20.md','REPORT.html'}
    result=dict(status='PASS_REQUESTED_SCOPE',mode=mode,created_utc=datetime.now(timezone.utc).isoformat(),
        python=sys.version,platform=platform.platform(),full_under200_run=False,stages=stages,
        source_sha256={str(p.relative_to(ROOT)).replace(os.sep,'/'):digest(p) for p in source_files()},
        output_sha256={name:digest(ROOT/name) for name in sorted(report_paths)})
    (ROOT/'results'/'run_record.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

def check_recorded():
    record=json.loads((ROOT/'results'/'run_record.json').read_text(encoding='utf-8'))
    require(record['full_under200_run'] is False,'unexpected full computation flag')
    for group in ('source_sha256','output_sha256'):
        for path,sha in record[group].items():
            require((ROOT/path).is_file() and digest(ROOT/path)==sha,'Recorded file changed: '+path)
    validate_reports()
    print('PASS: saved source/output hashes and all recorded certificate results match. No computations rerun.')

def main():
    if not __debug__ or sys.flags.optimize:
        raise RuntimeError('Run without -O or PYTHONOPTIMIZE: exact certificates use assertions.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-recorded',action='store_true',help='Verify saved result/file integrity only.')
    parser.add_argument('--collect-recorded',action='store_true',help='Assemble previously executed per-certificate reports; does not claim a new run.')
    args=parser.parse_args()
    if args.check_recorded:
        check_recorded();return
    (ROOT/'results').mkdir(exist_ok=True)
    if args.collect_recorded:
        build();receipt('assembled from separately executed certificate reports',[])
        print('PASS: previously executed reports assembled. Full CM computation not run.');return
    import sage.all
    initial={str(p):digest(p) for p in source_files()}
    stages=[]
    for label,path in STAGES:
        print('Running '+label+' ...',flush=True)
        started=time.monotonic()
        if path.endswith('.sage'):
            code="from pathlib import Path; from sage.repl.preparse import preparse_file; p=Path("+repr(str(ROOT/path))+"); ns={'__file__':str(p),'__name__':'__main__'}; exec('from sage.all import *',ns); exec(compile(preparse_file(p.read_text()),str(p),'exec'),ns)"
            command=[sys.executable,'-c',code]
        else:
            command=[sys.executable,str(ROOT/path)]
        with (ROOT/'results'/(label+'.log')).open('w',encoding='utf-8') as log:
            child=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        row=dict(certificate=label,returncode=child.returncode,seconds=round(time.monotonic()-started,3))
        stages.append(row)
        if child.returncode:
            (ROOT/'results'/'run_record.json').write_text(json.dumps(dict(status='INCOMPLETE',stages=stages),indent=2))
            raise RuntimeError(label+' failed; inspect results/'+label+'.log')
        print('PASS '+label,flush=True)
    # Families writes a generated report in data/; exclude that one output.
    for path,sha in initial.items():
        if Path(path)!=ROOT/'data'/'families_certificate.json':
            require(digest(Path(path))==sha,'Input changed during computation: '+path)
    build();receipt('fresh sequential execution of all listed certificates',stages)
    print('PASS: all certificates complete. CM table: sanity only. Open REPORT.html.')

if __name__=='__main__': main()
