"""Section5: Phi : Jac(X_N) -> E for all twelve forms.
For Recurrence formulas for Taylor expansions at CM points on genus-one and genus-two Shimura curves.

From this folder: sage "Compute Phi and CM Images.sage" --sanity       (small checks only)
                 sage "Compute Phi and CM Images.sage" --full         (compute all |D|<200 images and write the HTML table)
                 sage "Compute Phi and CM Images.sage" --table        (rebuild HTML from existing results)
Or load('Compute Phi and CM Images.sage'); run_all(labels=('14a1',)). Loading never starts a batch.
Numerical recognition is not an interval proof of analytic CM images.
"""
from sage.all import *
from pathlib import Path
from sage.repl.load import load as _phi_load
import json,sys

REPRO_ROOT=Path(globals().get('PHI_ROOT',Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd())).resolve()
if not (REPRO_ROOT/'Compute CM Expansion Series.sage').exists() and (REPRO_ROOT/'Section_5 Compute Phi and CM Images/Compute CM Expansion Series.sage').exists(): REPRO_ROOT=REPRO_ROOT/'Section_5 Compute Phi and CM Images'
assert (REPRO_ROOT/'Compute CM Expansion Series.sage').exists(), 'Start in the Section_5 Compute Phi and CM Images folder or set PHI_ROOT before loading.'
_sides=json.loads((REPRO_ROOT/'Fundamental Domain Side Data.json').read_text())
def side_data(level): return dict(_sides[str(level)])
for _module in ('Integrate Divisors and Recognize Images.sage','Compute CM Expansion Series.sage','Match Periods and Normalizations.sage','CM Quaternion and Embedding Data.sage','Compute and Validate CM Images.sage'):
    _phi_load(str(REPRO_ROOT/_module),globals())

CM_LABELS=('14a1','15a1','21a1','26a1','26b1','33a1','34a1','38a1','38b1','46a1','58a1','58b1')
CM_GENUS_TWO={'26a1','26b1','38a1','38b1','58a1','58b1'}
CM_BOUND=200
CM_PAPER='Recurrence formulas for Taylor expansions at CM points on genus-one and genus-two Shimura curves'
CM_TABLE5=[('14a1', -163, 'w2', ['144420233907863660158612249', '150338180418782421668964075']), ('15a1', -127, 'w3', ['563606776197313', '-590593001260921', '305940630600433', '-84764507545775', '11149788901250', '5166015625']), ('58a1', -131, 'w29', ['-88687', '211985', '-185426', '139246', '-156927', '83521'])]
CM_PASS_SPECS={
    '14a1':(1000,1200,1000),'15a1':(1000,1200,1000),
    '21a1':(2400,2600,1000),'26a1':(1000,1200,1000),
    '26b1':(3000,3200,1600),'33a1':(4000,4400,1000),
    '34a1':(2400,2800,1000),'38a1':(1600,1800,1000),
    '38b1':(3000,3200,1000),'46a1':(4000,4400,1000),
    '58a1':(4000,4400,1000),'58b1':(4000,4400,1000)}
CM_NATIVE_TARGETS={
    '14a1':('14a2',2),'15a1':('15a1',1),'21a1':('21a2',2),
    '26a1':('26a1',1),'26b1':('26b2',7),'33a1':('33a1',1),
    '34a1':('34a3',3),'38a1':('38a1',1),'38b1':('38b2',5),
    '46a1':('46a2',2),'58a1':('58a1',1),'58b1':('58b2',5)}



def cm_output_dir(output_dir=None):
    return (REPRO_ROOT/'results' if output_dir is None else Path(output_dir)).resolve()

def cm_result_path(label,output_dir=None):
    if label not in CM_LABELS: raise ValueError('Unknown target: '+str(label))
    return cm_output_dir(output_dir)/('x%s_cm_images_%s_result.sobj'%(label[:-2],label))

def load_result(label,output_dir=None):
    return load(str(cm_result_path(label,output_dir)))

def run_all(labels=CM_LABELS,pass_specs=None,output_dir=None):
    """Full calculation, or resume its generated checkpoints; never called by sanity_check."""
    labels=(labels,) if isinstance(labels,str) else tuple(labels)
    if not labels or len(labels)!=len(set(labels)) or not set(labels)<=set(CM_LABELS):
        raise ValueError('Choose distinct labels from CM_LABELS.')
    specs=dict(CM_PASS_SPECS);specs.update(pass_specs or {})
    destination=cm_output_dir(output_dir);destination.mkdir(parents=True,exist_ok=True)
    results={};coverage=expected_coverage()
    for label in labels:
        if label in CM_GENUS_TWO:
            results[label]=run_genus_two(labels=(label,),out_dir=str(destination),pass_specs={label:specs[label]})[label]
        else:
            N=int(label[:-2])
            results[label]=run_six_curves(levels=(N,),output_dir=str(destination),pass_specs={N:specs[label]})[N]
        require_complete_result(label,results[label])
        validate_result(label,results[label],coverage[str(int(label[:-2]))])
        _check_truncations(results[label],destination)
    return results

def cm_divisor_point(label,D,operator,orientation_index=1,class_index=1,output_dir=None):
    """Point over Q(x,y), coordinate fields, inclusion into H_D, and branch."""
    result=load_result(label,output_dir)
    key=(ZZ(D),ZZ(orientation_index),operator)
    item=next(r for r in result['numeric_rows'] if
              (r['D'],r['orientation_index'],r['operator'])==key and r['class_index']==class_index)
    row=next(r for r in result['rows'] if (r['D'],r['orientation_index'],r['operator'])==key)
    if row['status']=='identity':return dict(status='identity',point=EllipticCurve(label)(0),field=QQ)
    if row['status']!='recognized':raise ValueError('Unresolved entry: '+str(row['status']))
    exact=next(r for r in result['exact_rows'] if (r['D'],r['orientation_index'],r['operator'])==key)
    match=next(m for m in exact['numeric_matches'] if m['orbit_index']==item['orbit_index'])
    point=exact['points'][int(match['point_index'])-1]
    return dict(point,status='recognized',ambient_class_field=exact['ambient_class_field'],
                complex_embedding=exact['complex_embedding'],numerical_match=match)


def compact_point(entry):
    """Optionally simplify Q(x,y)'s presentation, preserving its inclusion."""
    if entry['status']=='identity':return dict(entry)
    P=entry['point'];L=entry['point_field']['field']
    if L is QQ:
        K=QQ;Q=P;inclusion=entry['point_field']['inclusion']
    else:
        K,fromK,toK=L.optimized_representation()
        Q=EllipticCurve(K,[QQ(a) for a in P.curve().ainvs()])([toK(P[0]),toK(P[1])])
        inclusion=entry['point_field']['inclusion']*fromK
    PH=entry['ambient_point']
    assert inclusion(Q[0])==PH[0] and inclusion(Q[1])==PH[1]
    return dict(point=Q,field=K,inclusion=inclusion,ambient_point=PH,
                complex_embedding=entry['complex_embedding']*inclusion)


def sanity_check():
    """20 terms for 12 forms, exact CM inputs and the 3 printed polynomial claims."""
    from sage.repl.preparse import preparse_file
    for source in REPRO_ROOT.glob('*.sage'):
        compile(preparse_file(source.read_text()),str(source),'exec')
    inputs=check_cm_inputs()
    checks=[]
    for N in (14,15,21,33,34,46):
        data=six_load_inputs(N);setup=data['builder'](data['side_path'],20,128)
        assert len(setup['b'])==20 and setup['b'][0]==1 and setup['tau0'].imag()>0
        record=data['records'][0];CC=ComplexField(128);RF=RealField(128)
        M=data['matrix'](record['mu_coords'],RF);tau=_cm_fixed_point_upper_half_plane(M,CC)
        assert tau.imag()>0 and abs(_six_mob(M,tau)-tau)<RF(2)^(-60)
        checks.append(str(N)+'a1')
    six33_exact_al_audit(six_load_inputs(33))
    E=EllipticCurve('33a1');T=E(0,0)
    assert six33_translate_point(E,{},source_is_identity=True)==T
    assert six33_translate_point(E,{'x':T[0],'y':T[1]},target_is_identity=True).is_zero()
    EC=E.change_ring(QQbar);P=EC.lift_x(QQbar(1))
    assert six33_translate_point(EC,{'x':P[0],'y':P[1]})==P+EC(T)
    for label in ('26a1','26b1','38a1','38b1','58a1','58b1'):
        genus_two_sanity_check(label,Nper=20,prec_bits=128);checks.append(label)
    assert set(checks)==set(CM_LABELS)
    R=PolynomialRing(QQ,'x')
    # Table5 claims; coefficients are constant-first. These checks do not
    # recompute the analytic divisor images whose recognition is claimed.
    printed=CM_TABLE5
    for label,D,operator,coefficients in printed:
        f=R(coefficients)
        assert f.degree()==QuadraticField(D).class_number() and f.is_irreducible()
        assert gcd([ZZ(c) for c in f.list()])==1 and f.leading_coefficient()>0
        assert operator in cm_inputs(int(label[:-2]))['CM_AL_COORDS']
    print('PASS: 20-term setups for all 12 forms; exact CM inputs for 5616 divisors; 3 Table 5 polynomial sanity checks.')
    print('Full under-200 images, period matching and recognition were NOT run.')
    return dict(forms=checks,CM_inputs=inputs,printed_polynomials=3,full_under200_run=False)

def write_cm_table(results=None,output_dir=None,path=None):
    """Validate all saved images and write the standalone, five-column HTML table."""
    from html import escape
    from sage.env import SAGE_VERSION
    import hashlib
    destination=cm_output_dir(output_dir)
    coverage=expected_coverage(); R=PolynomialRing(QQ,'x'); rows=[]; summaries=[]
    for label in CM_LABELS:
        result=load_result(label,output_dir) if results is None else results[label]
        require_complete_result(label,result)
        expected=coverage[str(int(label[:-2]))]
        summaries.append(validate_result(label,result,expected))
        _check_truncations(result,destination)
        b2,b4,b6,b8=EllipticCurve(label).b_invariants()
        torsion=R([b6,2*b4,b2,4])
        for item in expected:
            D,h=item['D'],item['h']
            for w in sorted(cm_inputs(int(label[:-2]))['CM_AL_COORDS'],key=lambda w:int(w[1:])):
                copies=[r for r in result['rows'] if int(r['D'])==D and r['operator']==w]
                assert len(copies)==item['orientations']
                signatures={(r['status'],tuple(r['polynomial'].list()) if r['status']=='recognized' else ()) for r in copies}
                assert len(signatures)==1,(label,D,w,'Repeated images differ')
                status,coefficients=signatures.pop(); annotation=''; terms=[]; plain=[]
                if status=='identity':
                    value='<i>P</i> = <span class="origin">𝒪</span>'; text='P = O (point at infinity)'
                else:
                    f=R(coefficients)
                    assert f.degree()>=1 and f.leading_coefficient()>0 and f.is_irreducible()
                    assert all(c in ZZ for c in coefficients) and gcd(coefficients)==1
                    if f.degree()==1:
                        q=-f[0]/f[1]; a,b=abs(q.numerator()),q.denominator()
                        value='<i>x</i> = '+('&minus;' if q<0 else '')
                        value+=(str(a) if b==1 else '<span class="fraction"><span class="numerator">%s</span><span class="denominator">%s</span></span>'%(a,b))
                        text='x = '+str(q)
                    else:
                        for n in reversed(range(len(coefficients))):
                            c=coefficients[n]
                            if not c: continue
                            sign=('-' if c<0 else '') if not terms else (' - ' if c<0 else ' + ')
                            number='' if n and abs(c)==1 else str(abs(c))
                            terms.append('<span class="term">'+sign.replace('-','&minus;')+number+('<i>x</i>' if n else '')+('<sup>%s</sup>'%n if n>1 else '')+'</span>')
                            plain.append(sign+number+('*' if n and number else '')+('x' if n else '')+('^%s'%n if n>1 else ''))
                        value='<wbr>'.join(terms); text=''.join(plain)
                    if f.divides(torsion):
                        annotation='2-torsion'; value+=' <span class="annotation">(2-torsion)</span>'; text+=' (2-torsion)'
                rows.append(dict(E=label,D=D,h=h,w=int(w[1:]),status=status,
                    coefficients=[str(c) for c in coefficients],images=h*len(copies),
                    value=value,text=text,annotation=annotation))
    assert len(rows)==720 and sum(r['images'] for r in rows)==5616
    for label,D,w,coefficients in CM_TABLE5:
        row=next(r for r in rows if (r['E'],r['D'],r['w'])==(label,D,int(w[1:])))
        assert row['status']=='recognized' and row['coefficients']==coefficients,(label,'Table 5')
    def curve(label): return label[:-1]+'<sub>'+label[-1]+'</sub>'
    parts=[]
    for label in CM_LABELS:
        parts.append('<tbody id="curve-%s">'%label)
        selected=[r for r in rows if r['E']==label]
        for i,r in enumerate(selected):
            classes=('new-discriminant' if i%3==0 else '')+(' new-curve' if i==0 else '')
            parts.append('<tr class="%s" data-curve="%s" data-D="%s" data-w="%s" data-images="%s">'%(classes,label,r['D'],r['w'],r['images']))
            if i%3==0:
                parts.append('<td rowspan="3" class="curve">%s</td><td rowspan="3" class="discriminant">&minus;%s</td><td rowspan="3" class="class-number">%s</td>'%(curve(label),abs(r['D']),r['h']))
            long=' long-expression' if any(len(c.lstrip('-'))>60 for c in r['coefficients']) else ''
            parts.append('<td class="operator"><i>w</i><sub>%s</sub></td><td class="expression%s" aria-label="%s">%s</td></tr>'%(r['w'],long,escape(r['text'],quote=True),r['value']))
        parts.append('</tbody>')
    navigation=' · '.join('<a href="#curve-%s">%s</a>'%(s,curve(s)) for s in CM_LABELS)
    digest=hashlib.sha256(json.dumps([(r['E'],r['D'],r['h'],r['w'],r['status'],r['coefficients']) for r in rows],separators=(',',':')).encode()).hexdigest()
    document='''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Section 5 — Table of all CM images</title><style>
:root { color-scheme: light; }
  * { box-sizing: border-box; }
  body { margin: 0; padding: 36px 30px 48px; color: #191919; background: #fff;
         font-family: Georgia, 'Times New Roman', serif; }
  main { max-width: 1500px; margin: 0 auto; }
  .section { margin: 0 0 10px; font: 12px/1.4 Arial, sans-serif; letter-spacing: .13em; text-transform: uppercase; color: #555; }
  h1 { font-size: 30px; line-height: 1.2; font-weight: normal; margin: 0 0 12px; }
  .description { max-width: 1180px; font-size: 15px; line-height: 1.65; margin: 0 0 6px; }
  .navigation { margin: 18px 0 24px; font-size: 17px; line-height: 1.8; }
  .navigation a { color: #243f63; text-decoration-thickness: 1px; text-underline-offset: 3px; }
  .navigation > span { margin: 0 7px; color: #999; }
  .math, table, .curve-label { font-family: 'Cambria Math', 'Times New Roman', serif; }
  .table-wrap { overflow-x: auto; }
  table { border-collapse: collapse; width: 100%; min-width: 900px; table-layout: fixed; border: 1px solid #777; }
  col.e { width: 88px; } col.d { width: 78px; } col.h { width: 44px; } col.w { width: 66px; }
  th { font-size: 20px; font-weight: normal; padding: 9px 12px; border: 1px solid #777;
       border-bottom: 3px double #444; background: #fff; }
  th:last-child { text-align: left; }
  td { border: 1px solid #b0b0b0; padding: 8px 12px; vertical-align: middle; }
  td.curve, td.discriminant, td.class-number, td.operator { text-align: center; font-size: 18px; white-space: nowrap; }
  .expression { font-size: 18px; line-height: 1.6; overflow-wrap: anywhere; }
  .long-expression { font-size: 16px; }
  .term { display: inline-block; max-width: 100%; overflow-wrap: anywhere; }
  sub { font-size: .72em; vertical-align: -.2em; }
  sup { font-size: .72em; vertical-align: .45em; line-height: 0; }
  .fraction { display: inline-flex; flex-direction: column; text-align: center; vertical-align: middle;
              line-height: 1.1; margin: 1px 2px; font-size: .94em; }
  .numerator { border-bottom: 1px solid currentColor; padding: 0 3px 2px; }
  .denominator { padding: 2px 3px 0; }
  .annotation { font-size: .92em; white-space: nowrap; }
  .origin { font-family: 'Cambria Math', 'Times New Roman', serif; }
  .new-discriminant td { border-top-color: #777; }
  .new-curve td { border-top: 3px double #333; }
  tbody { scroll-margin-top: 18px; }
  footer { margin-top: 24px; max-width: 1200px; font-size: 13px; line-height: 1.65; color: #444; }
  footer p { margin: 6px 0; }
  .hash { overflow-wrap: anywhere; font: 11px/1.6 monospace; color: #666; }
  @media (max-width: 700px) { body { padding: 22px 14px; } h1 { font-size: 25px; } }
  @page { size: A4 landscape; margin: 12mm; }
  @media print {
    body { padding: 0; color: #000; }
    h1 { font-size: 20pt; }
    .section { font-size: 8pt; }
    .description { font-size: 10pt; }
    .navigation { display: none; }
    .table-wrap { overflow: visible; }
    table { min-width: 0; }
    thead { display: table-header-group; }
    th { font-size: 13pt; padding: 5px 8px; }
    td { padding: 4px 8px; }
    td.curve, td.discriminant, td.class-number, td.operator { font-size: 11pt; }
    .expression { font-size: 11pt; line-height: 1.5; }
    .long-expression { font-size: 10pt; }
    tr { break-inside: avoid; }
    footer { font-size: 9pt; }
    a { color: #000; text-decoration: none; }
  }
</style></head><body><main>
<p class="section">Section 5 · Numerical CM images</p><h1>Images of CM divisors</h1>
<p class="description">All admissible negative fundamental discriminants with <span class="math">|<i>D</i>| &lt; 200</span>.
The 720 rows below represent all 5,616 computed divisor images on twelve elliptic curves; repeated entries are listed once.</p>
<p class="description math"><i>P</i><sub><i>D</i>,<i>w</i></sub> = Φ([<i>w</i>τ<sub><i>D</i></sub>] − [τ<sub><i>D</i></sub>]).
Here <i>h</i> is the class number, <i>w</i><sub><i>m</i></sub> is the Atkin–Lehner involution, and 𝒪 is the point at infinity.</p>
<nav class="navigation" aria-label="Jump to curve">__NAVIGATION__</nav>
<div class="table-wrap"><table aria-label="All computed CM divisor images">
<colgroup><col class="e"><col class="d"><col class="h"><col class="w"><col></colgroup>
<thead><tr><th scope="col"><i>E</i></th><th scope="col"><i>D</i></th><th scope="col"><i>h</i></th><th scope="col"><i>w</i></th>
<th scope="col"><i>x</i>(<i>P</i><sub><i>D</i>,<i>w</i></sub>) or its minimal polynomial</th></tr></thead>
__ROWS__
</table></div><footer>
<p>Curve names are Cremona labels. Rational values are displayed exactly; higher-degree entries are primitive integral multiples of minimal polynomials, with positive leading coefficient. The “2-torsion” annotation is checked using the curve’s 2-division polynomial.</p>
<p>Companion to <cite>__PAPER__</cite>, by Jingxuan Geng. Analytic CM images are numerically recognized; the saved results pass exact curve and field checks and comparisons at two truncations. Minimum observed coordinate agreement: __BITS__ bits, using |z_num−z_exact| / max(1, |z_num|, |z_exact|). This is a numerical agreement measure.</p>
<p>Generated by <code>Compute Phi and CM Images.sage --table</code> using SageMath __VERSION__. Use <code>--full</code> to compute the results and table from scratch.</p>
<p class="hash">Ordered table data SHA-256: __DIGEST__</p>
</footer></main></body></html>
'''
    for key,value in dict(NAVIGATION=navigation,ROWS='\n'.join(parts),PAPER=escape(CM_PAPER),
        BITS=str(min(s['min_xy_match_bits'] for s in summaries)),VERSION=escape(SAGE_VERSION),DIGEST=digest).items():
        document=document.replace('__'+key+'__',value)
    path=REPRO_ROOT/'Section_5 Table of All CM Images.html' if path is None else Path(path)
    path.write_text(document,encoding='utf-8')
    print('PASS: 720 table rows, 5616 images, all 3 current Table 5 examples; wrote',path)
    return path


if __name__=='__main__':
    if '--sanity' in sys.argv: sanity_check()
    elif '--full' in sys.argv: write_cm_table(run_all())
    elif '--table' in sys.argv: write_cm_table()
