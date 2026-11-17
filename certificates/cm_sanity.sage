"""Only the explicitly requested small CM-table sanity check; never run_all."""
from pathlib import Path
import os, json, time
root = Path(__file__).resolve().parents[1]
os.chdir(root / 'twelve_under200_reproduce')
started = time.time()
load('sanity_check.sage')
report = sanity_check()
assert report['passed'] and report['full_computation_run'] is False
report['seconds'] = round(time.time() - started, 3)
report['scope'] = 'Exact supplied CM input checks and 20-term setups only; no numerical CM image table reproduced.'
def clean(value):
    if isinstance(value, dict): return {str(k): clean(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [clean(v) for v in value]
    if isinstance(value, (str,bool,int,float)) or value is None: return value
    if value in ZZ: return int(value)
    return str(value)
(root / 'results' / 'cm_sanity.json').write_text(json.dumps(clean(report), indent=2) + '\n')
print('PASS: CM sanity only; full under-200 run was not invoked.')
