"""Execute the two supplied exact Sage notebooks from clean kernels."""
import argparse
import json
from pathlib import Path
import time
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]

def run_certificate():
    import nbformat
    from nbclient import NotebookClient
    from jupyter_client import KernelManager
    from jupyter_client.kernelspec import KernelSpecManager
    # Bind the notebook to the same Sage Python used to launch this script,
    # rather than an unrelated global Jupyter kernel registration.
    temp = tempfile.TemporaryDirectory(prefix='shimura-kernel-')
    kernel_dir = Path(temp.name) / 'certificate-sage'
    kernel_dir.mkdir()
    (kernel_dir / 'kernel.json').write_text(json.dumps(dict(
        argv=[sys.executable, '-m', 'sage.repl.ipython_kernel', '-f', '{connection_file}'],
        display_name='Certificate Sage', language='sage')))
    records = []
    recovered = {}
    for source in sorted((ROOT / 'inputs').glob('*.ipynb')):
        started = time.time()
        notebook = nbformat.read(source, as_version=4)
        for cell in notebook.cells:
            if cell.cell_type == 'code':
                cell.outputs = []
                cell.execution_count = None
        if 'parameter_recovery' in source.name:
            export_code = '''import json as _json
_certificate_data = {label: dict(
    H=[str(t) for t in case['H'].list()],
    R=str(results[label]['R']).replace('^','**'),
    L=str(results[label]['solution']['L']),
    u0=str(case['u0']),u1=str(results[label]['solution']['u1']),
    c=[str(t) for t in case['c']]) for label,case in cases.items()}
print('CERTIFICATE_DATA='+_json.dumps(_certificate_data))'''
        else:
            export_code = '''import json as _json
_certificate_data = {case['label']: dict(
    H=[str(t) for t in Pu(case['H']).list()],
    R=str(case['R']).replace('^','**'),L=str(case['L']),
    u0=str(case['u0']),u1=str(case['u1']),
    c=['1','1',str(case['c2']),str(case['c3'])]) for case in cases}
print('CERTIFICATE_DATA='+_json.dumps(_certificate_data))'''
        notebook.cells.append(nbformat.v4.new_markdown_cell(
            '## Package integration record\nThis cell exports the values just derived above for exact comparison with the other certificates.'))
        notebook.cells.append(nbformat.v4.new_code_cell(export_code))
        if notebook.nbformat_minor < 5:
            for cell in notebook.cells:
                cell.pop('id', None)
        output = ROOT / 'results' / source.name
        def progress(cell_index, **kwargs):
            print(f'{source.stem}: cell {cell_index + 1}/{len(notebook.cells)}', flush=True)
        manager = KernelManager(kernel_name='certificate-sage',
            kernel_spec_manager=KernelSpecManager(kernel_dirs=[temp.name]))
        client = NotebookClient(notebook, km=manager, timeout=900, startup_timeout=600,
            resources={'metadata': {'path': str(ROOT)}}, allow_errors=False,
            on_cell_executed=progress)
        try:
            client.execute()
        finally:
            nbformat.write(notebook, output)
        errors = [o for c in notebook.cells if c.cell_type == 'code'
                  for o in c.outputs if o.output_type == 'error']
        assert not errors, errors
        assert all(c.execution_count is not None for c in notebook.cells if c.cell_type == 'code')
        nbformat.validate(notebook)
        for output_cell in notebook.cells[-1].outputs:
            if output_cell.output_type == 'stream':
                for line in output_cell.text.splitlines():
                    if line.startswith('CERTIFICATE_DATA='):
                        recovered.update(json.loads(line.split('=',1)[1]))
        records.append(dict(notebook=source.name, passed=True,
            code_cells=sum(c.cell_type == 'code' for c in notebook.cells),
            seconds=round(time.time() - started, 3),
            evidence='Fresh Sage kernel; every code cell including answer comparisons executed.'))
    assert len(records) == 2
    assert len(recovered) == 12
    (ROOT / 'results' / 'recovered_parameters.json').write_text(json.dumps(recovered, indent=2) + '\n')
    report = dict(passed=True, kind='exact algebra conditional on stated initial coefficients and geometry',
                  notebooks=records)
    (ROOT / 'results' / 'notebooks.json').write_text(json.dumps(report, indent=2) + '\n')
    temp.cleanup()
    return report

if __name__ == '__main__':
    print(json.dumps(run_certificate(), indent=2))
