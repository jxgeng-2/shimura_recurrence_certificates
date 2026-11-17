"""Assemble checked certificate reports into a readable package report."""
from decimal import Decimal
import hashlib
import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPORTS = {
    'General identities': 'results/general_identities.json',
    'Table 3 families': 'data/families_certificate.json',
    'Section 4 notebooks': 'results/notebooks.json',
    'All twelve expansions': 'results/coefficients.json',
    'Appendix C': 'results/appendix_c.json',
    'Nelson exact inputs': 'results/nelson_inputs.json',
    'Nelson numerical agreement': 'results/numerical.json',
    'Section 5 selected polynomials': 'results/cm_claims.json',
    'CM source sanity': 'results/cm_sanity.json',
    'Cross-certificate consistency': 'results/integration.json',
}

def require(condition, message):
    if not condition:
        raise ValueError(message)

def validate_reports():
    reports={}
    for label,path in REPORTS.items():
        value=json.loads((ROOT/path).read_text(encoding='utf-8'))
        require(value.get('passed') is True or str(value.get('status','')).startswith('PASS'),
                label+': no passing report')
        reports[label]=value
    num=reports['Nelson numerical agreement']
    require(len(num['forms'])==12 and num['interval_certified'] is False, 'numerical scope')
    require(num['lattice_cutoffs']==[144,196], 'two expected numerical cutoffs')
    require(len(reports['All twelve expansions']['forms'])==12, 'coefficient coverage')
    require(reports['Appendix C']['coordinate_row_count']==24, 'Appendix C coverage')
    require(reports['CM source sanity']['full_computation_run'] is False, 'full CM run is outside requested scope')
    require(len(reports['CM source sanity']['forms'])==12, 'CM sanity coverage')
    require(len(reports['Cross-certificate consistency']['forms'])==12, 'integration coverage')
    for dep in reports['Cross-certificate consistency']['dependencies']:
        require(hashlib.sha256((ROOT/dep['path']).read_bytes()).hexdigest()==dep['sha256'],
                'Integration dependency changed: '+dep['path'])
    for name in ('shimura_parameter_recovery.ipynb','shimura_schwarzian_14_15_21.ipynb'):
        nb=json.loads((ROOT/'results'/name).read_text(encoding='utf-8'))
        cells=[c for c in nb['cells'] if c['cell_type']=='code']
        require(all(c['execution_count'] is not None for c in cells),'unexecuted notebook cell')
        require(not any(o['output_type']=='error' for c in cells for o in c['outputs']), 'notebook error output')
    for form in num['forms']:
        require([p['cutoff'] for p in form['passes']]==[144,196], 'missing numerical pass')
        for p in form['passes']:
            require(all(Decimal(c['observed_error'])<Decimal('1e-40') for c in p['coefficients']), 'coefficient tolerance')
            require(Decimal(p['normalized_alpha_observed_error'])<Decimal('1e-40'),'alpha tolerance')
            require(Decimal(p['theta_observed_error'])<Decimal('1e-40'),'theta tolerance')
    return reports

def build():
    reports=validate_reports()
    num=reports['Nelson numerical agreement']
    rows=[]
    for f in num['forms']:
        low,high=f['passes']
        rows.append(dict(form=f['form'],
            coefficient_max_cutoff144=str(max(Decimal(c['observed_error']) for c in low['coefficients'])),
            coefficient_max_cutoff196=str(max(Decimal(c['observed_error']) for c in high['coefficients'])),
            theta_max_cutoff196=high['theta_observed_error'],
            normalized_alpha_max_cutoff196=high['normalized_alpha_observed_error']))
    summary=dict(status='PASS_REQUESTED_SCOPE',full_under200_run=False,
        exact_checks='All exact and conditional algebraic checks passed.',
        numerical_checks='Observed finite-sum agreement, not a rigorous enclosure of analytic values.',
        forms=12,coefficient_comparisons=43,cm_coordinate_rows=24,cm_input_divisors=5616,
        max_coefficient_discrepancy_cutoff196=str(max(Decimal(r['coefficient_max_cutoff196']) for r in rows)),
        max_theta_discrepancy_cutoff196=str(max(Decimal(r['theta_max_cutoff196']) for r in rows)),
        max_normalized_alpha_discrepancy_cutoff196=str(max(Decimal(r['normalized_alpha_max_cutoff196']) for r in rows)),
        reports=REPORTS,numerical_by_form=rows)
    # Derive displayed powers of ten from this run, so a future rerun cannot
    # inherit a stronger threshold merely because an older run achieved it.
    summary['observed_uniform_thresholds']={}
    for label,key in [('coefficients','max_coefficient_discrepancy_cutoff196'),
                      ('theta','max_theta_discrepancy_cutoff196'),
                      ('normalized_alpha','max_normalized_alpha_discrepancy_cutoff196')]:
        value=Decimal(summary[key])
        require(value>0,'Cannot infer an agreement threshold from a rounded zero.')
        threshold=Decimal(1).scaleb(value.adjusted()+1)
        require(value<threshold,'Displayed threshold must exceed observed discrepancy.')
        summary['observed_uniform_thresholds'][label]='1e'+str(value.adjusted()+1)
    (ROOT/'results'/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    esc=lambda x:html.escape(str(x))
    sections=['<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
        '<title>Shimura recurrence certificate</title><style>body{font:16px/1.6 system-ui,sans-serif;color:#172c36;max-width:1100px;margin:50px auto;padding:0 24px;background:#f7f9fa}h1{font-size:36px;line-height:1.15}h2{margin-top:42px}a{color:#185b86}table{border-collapse:collapse;width:100%;background:white;margin:20px 0}th,td{text-align:left;padding:10px 14px;border-bottom:1px solid #d8e1e5;overflow-wrap:anywhere}th{background:#e8eff2}code{font-size:.88em;overflow-wrap:anywhere}details{background:white;border:1px solid #d8e1e5;margin:12px 0;padding:12px 18px}summary{font-weight:650;cursor:pointer}.tag{color:#246746;font-weight:650}.note{border-left:4px solid #247c86;padding:10px 18px;background:#edf4f4}.numeric{font-variant-numeric:tabular-nums}footer{margin-top:48px;font-size:14px;color:#506674}@media print{body{background:white;max-width:none;margin:0}details{break-inside:avoid}}</style>',
        '<body><p class="tag">COMPLETED · ALL REQUESTED CHECKS PASS</p><h1>Code certificate for<br>Shimura-curve recurrence formulas</h1>',
        '<p>Twelve eigenforms · nine curves · Tables 1–7 · Sections 2–5 and Appendices A–C.</p>',
        '<p><a href="README.md">Run instructions</a> · <a href="docs/coverage.md">Claim-by-claim coverage</a> · <a href="inputs/draft.pdf">Reference draft</a> · <a href="results/summary.json">Machine-readable summary</a></p>',
        '<div class="note">The exact checks certify algebraic identities and consequences of the stated geometric and coefficient inputs. Nelson comparisons establish observed numerical agreement. The twelve-under-200 source was checked only with its small sanity routine; the full image computation was not run.</div>',
        '<h2>Numerical agreement achieved</h2><p>90 decimal digits; weight target 80 digits; lattice cutoffs 144 and 196. Across the 43 assumed coefficients, the final cutoff gives:</p><table><tr><th>Quantity</th><th>Largest observed discrepancy</th><th>Uniform observed threshold</th></tr>']
    for label,key,threshold_key in [('cₙ','max_coefficient_discrepancy_cutoff196','coefficients'),('Θ','max_theta_discrepancy_cutoff196','theta'),('Θ/(πΩ²)','max_normalized_alpha_discrepancy_cutoff196','normalized_alpha')]:
        exponent=summary['observed_uniform_thresholds'][threshold_key].split('e')[1]
        sections.append(f'<tr><td>{label}</td><td>{Decimal(summary[key]):.8E}</td><td>&lt; 10<sup>{exponent}</sup></td></tr>')
    sections.extend(['</table><p>These thresholds describe the finite computations; they are not interval-certified bounds for the infinite Nelson sum. See <a href="docs/numerical_scope.md">numerical scope</a>.</p>',
        '<table><tr><th>Form</th><th>max |Δc|, cutoff 144</th><th>max |Δc|, cutoff 196</th><th>|ΔΘ/(πΩ²)|, cutoff 196</th></tr>'])
    for r in rows:
        sections.append('<tr><td>'+esc(r['form'])+'</td>'+''.join(f'<td>{Decimal(r[k]):.3E}</td>' for k in ('coefficient_max_cutoff144','coefficient_max_cutoff196','normalized_alpha_max_cutoff196'))+'</tr>')
    sections.append('</table><h2>Executed certificates</h2><table><tr><th>Certificate</th><th>Evidence</th></tr>')
    for label,path in REPORTS.items():
        sections.append(f'<tr><td>{esc(label)}</td><td><a href="{path}">Passing report</a></td></tr>')
    sections.extend(['</table><p><a href="results/shimura_parameter_recovery.ipynb">Executed parameter-recovery notebook</a> · <a href="results/shimura_schwarzian_14_15_21.ipynb">Executed Schwarzian notebook</a></p>',
        '<h2>First 20 coefficients of every form</h2><p>For w=(z−τ)/(z−conjugate(τ)), the appendix normalization is f(z)=(1−w)² Σ cₙ(Θw)ⁿ/n!. Each list gives c₀,…,c₁₉. The API also returns ordinary q coefficients and coefficients after multiplying by (1−w)².</p>'])
    for form,exp in reports['All twelve expansions']['expansions'].items():
        sections.append(f'<details><summary>{esc(form)} — c₀ through c₁₉</summary><p>τ = <code>{esc(exp["tau"])}</code><br>Θ/(πΩ²) = <code>{esc(exp["theta_over_pi_Omega_D_squared"])}</code></p><table><tr><th>n</th><th>cₙ</th></tr>')
        sections.extend(f'<tr><td>{n}</td><td><code>{esc(c)}</code></td></tr>' for n,c in enumerate(exp['c']))
        sections.append('</table></details>')
    sections.append('<footer>Original external inputs remain unchanged. This package retains the computation sources, exact inputs, fresh results and references needed for the supplied draft. No full CM-image reproduction was performed.</footer></body></html>')
    (ROOT/'REPORT.html').write_text('\n'.join(sections),encoding='utf-8')
    return summary

if __name__=='__main__':
    build()
    print('PASS: all recorded reports validated; REPORT.html and results/summary.json generated.')
