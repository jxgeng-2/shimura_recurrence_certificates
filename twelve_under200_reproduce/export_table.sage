"""Full-range table export; no historical results or caches are required."""
from sage.all import *
from pathlib import Path
import csv,json,hashlib
if "REPRO_ROOT" not in globals():
    load("reproduce.sage")

def expected_coverage():
    coverage={}
    for N in sorted({int(s[:-2]) for s in CM_LABELS}):
        primes=ZZ(N).prime_divisors(); records=[]
        for d in range(3,CM_BOUND):
            D=ZZ(-d)
            if fundamental_discriminant(D)!=D:continue
            chi=[ZZ(kronecker(D,p)) for p in primes]
            if 1 in chi:continue
            h=ZZ(D).class_number(); orientations=prod(1-c for c in chi)
            records.append(dict(D=int(D),h=int(h),orientations=int(orientations),
                                cm_points=int(h*orientations)))
        coverage[str(N)]=records
    return coverage


def _twelve_key(r):
    return (int(r['D']),int(r['orientation_index']),r['operator'])


def validate_result(label,result,expected):
    """Check full coverage, stable recognition, exact fields, both coordinates and AL identities."""
    N=int(label[:-2]); E=EllipticCurve(label); primes=ZZ(N).prime_divisors()
    assert E.optimal_curve()==E
    operators=['w'+str(d) for d in (*primes,N)]
    rows={_twelve_key(r):r for r in result['rows']}
    exact={_twelve_key(r):r for r in result['exact_rows']}
    assert len(rows)==len(result['rows'])
    wanted={(r['D'],ori,cl,w) for r in expected
            for ori in range(1,r['orientations']+1)
            for cl in range(1,r['h']+1) for w in operators}
    actual={(int(r['D']),int(r['orientation_index']),int(r['class_index']),r['operator'])
            for r in result['numeric_rows']}
    assert actual==wanted,(label,'incomplete CM coverage',len(actual),len(wanted))
    assert len(actual)==len(result['numeric_rows'])
    class_numbers={r['D']:r['h'] for r in expected}
    assert all(int(r['h'])==class_numbers[int(r['D'])] for r in result['numeric_rows'])
    assert set(rows)=={(D,o,w) for D,o,c,w in wanted}
    assert all(r['E']==label and r['status'] in ('recognized','identity') for r in rows.values())
    stable_rows=0
    for r in rows.values():
        if r['status']!='recognized':continue
        rec=r['recognition']; windows=rec['algdep_stability_windows']
        assert len(windows)>=2 and rec['algdep_holdout_bits']>=30
        assert all(len(w)>=3 and all(int(w[i+1])-int(w[i])==1 for i in range(len(w)-1)) for w in windows)
        centers=[QQ(w[0]+w[-1])/2 for w in windows]
        assert max(centers)-min(centers)>=20
        assert rec['class_field_verified'] and rec['polynomial']==r['polynomial']
        stable_rows+=1
    matches=[]; points={}; finite=origin=0
    for item in result['numeric_rows']:
        key=_twelve_key(item); row=rows[key]
        pkey=(int(item['D']),int(item['orbit_index']),item['operator'])
        if row['status']=='identity':
            assert item['is_infinity']
            origin+=1; points[pkey]=None; continue
        assert not item['is_infinity']
        entry=exact[key]; H=entry['ambient_class_field']
        assert H.absolute_degree()==2*item['h']
        match=next(m for m in entry['numeric_matches'] if m['orbit_index']==item['orbit_index'])
        assert int(match['xy_match_bits'])>=max(50,int(row['recognition']['algdep_use_bits'])+30)
        point=entry['points'][int(match['point_index'])-1]
        P=point['point']; PH=point['ambient_point']; L=point['point_field']['field']
        assert list(P.curve().ainvs())==[L(a) for a in E.ainvs()]
        assert P.curve().is_on_curve(P[0],P[1]) and PH.curve().is_on_curve(PH[0],PH[1])
        inclusion=point['point_field']['inclusion']
        assert inclusion(P[0])==PH[0] and inclusion(P[1])==PH[1]
        assert PH[0].minpoly().monic()==row['polynomial'].change_ring(QQ).monic()
        emb=entry['complex_embedding']; coord_bits=[]
        for j,name in enumerate(('x','y')):
            v=emb(PH[j]); a=item[name]
            error=abs(v-a)/max(1,abs(v),abs(a))
            coord_bits.append(emb.codomain().precision() if error==0 else int(floor(-error.log2())))
        assert min(coord_bits)>=int(match['xy_match_bits'])
        matches.append(min(coord_bits)); finite+=1; points[pkey]=PH
    checked_triplets=0
    if label in CM_GENUS_TWO:
        kernel={'26a1':'w13','26b1':'w2','38a1':'w19','38b1':'w2',
                '58a1':'w58','58b1':'w2'}[label]
        other=[w for w in operators if w!=kernel]
        for D,idx in {(p[0],p[1]) for p in points}:
            assert points[D,idx,kernel] is None
            a=points[D,idx,other[0]];b=points[D,idx,other[1]]
            assert (a is None and b is None) or (a is not None and b is not None and a==b)
            checked_triplets+=1
    else:
        constant={14:'w7',15:'w5',21:'w3',33:'w11',34:'w2',46:'w23'}[N]
        constant_x={14:None,15:3,21:2,33:0,34:-2,46:None}[N]
        other=next('w'+str(p) for p in primes if 'w'+str(p)!=constant)
        for D,idx in {(p[0],p[1]) for p in points}:
            T=points[D,idx,constant]; P=points[D,idx,other]; Q=points[D,idx,'w'+str(N)]
            if constant_x is None:
                assert T is None
                assert (P is None and Q is None) or (P is not None and Q is not None and P==Q)
            else:
                assert T is not None and T[0]==constant_x and (2*T).is_zero()
                parent=T.curve(); PP=parent(0) if P is None else P; QQp=parent(0) if Q is None else Q
                assert QQp==PP+T
            checked_triplets+=1
    preconditioned=[r for r in result['numeric_rows'] if r.get('kernel_preconditioning',{}).get('selected') not in (None,'identity')]
    for item in preconditioned:
        info=item['kernel_preconditioning']
        assert info['original_tau']==item['tauD']
        assert min(info['original_xy_agreement_bits'])>=60
        for j,name in enumerate(('x','y')):
            a,b=info['original'][name],item[name]
            error=abs(a-b)/max(1,abs(a),abs(b))
            actual_bits=a.parent().precision() if error==0 else max(0,int(floor(-error.log2())))
            assert actual_bits>=info['original_xy_agreement_bits'][j]
    return dict(E=label,N=N,discriminants=len(expected),divisors=finite+origin,
                finite=finite,identity=origin,unresolved=0,min_xy_match_bits=min(matches),
                stable_recognition_rows=stable_rows,exact_AL_triplets=checked_triplets,
                kernel_preconditioned_evaluations=len(preconditioned),
                exact_curves_and_descents_checked=True,complete_class_orientation_counts=True)


def _twelve_write_csv(path,rows,columns):
    with open(path,'w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
        writer.writerows({k:r.get(k,'') for k in columns} for r in rows)


def _twelve_tex_polynomial(coeffs,max_chars=110):
    terms=[]
    for power in range(len(coeffs)-1,-1,-1):
        c=coeffs[power]
        if not c:continue
        sign=('-' if c<0 else '') if not terms else ('- ' if c<0 else '+ ')
        term=str(abs(c)) if power==0 or abs(c)!=1 else ''
        if power:term+='x' if power==1 else 'x^{%s}'%power
        terms.append(sign+term)
    lines=[]; line=''
    for term in terms:
        if line and len(line)+len(term)+1>max_chars:lines.append(line);line=term
        else:line+=(' ' if line else '')+term
    if line:lines.append(line)
    return lines[0] if len(lines)==1 else r'\begin{aligned}&'+r'\\&'.join(lines)+r'\end{aligned}'



EXPECTED_TABLE_SHA256='78e0bdc8e256bcb18e805156cddb89ab07f1207565b84bc90c5a312ac7b38bea'


def table_fingerprint(rows):
    canonical=[{k:r[k] for k in ('E','D','h','operator','status','polynomial_coefficients')} for r in rows]
    canonical.sort(key=lambda r:(r['E'],abs(r['D']),int(r['operator'][1:])))
    return hashlib.sha256(json.dumps(canonical,sort_keys=True,separators=(',',':'),default=int).encode()).hexdigest()


def _check_truncations(result,output_dir):
    meta=result['metadata']
    pairs={}
    for kind in ('low','high'):
        path=Path(meta['numeric_'+kind+'_cache'])
        payload=load(str(path if path.is_absolute() else output_dir/path))
        pairs[kind]={(int(r['D']),int(r['orbit_index']),r['operator']):r for r in payload['rows']}
        assert len(pairs[kind])==len(payload['rows'])
    expected={(int(r['D']),int(r['orbit_index']),r['operator']) for r in result['numeric_rows']}
    assert set(pairs['low'])==set(pairs['high'])==expected
    for r in result['numeric_rows']:
        key=(int(r['D']),int(r['orbit_index']),r['operator'])
        a,b=pairs['low'][key],pairs['high'][key]
        assert r['x']==b['x'] and r['y']==b['y']
        if r['is_infinity']:
            assert min(a['identity_bits'],b['identity_bits'])>=70
            continue
        measured=[]
        for name in ('x','y'):
            x,y=a[name],b[name]
            error=abs(x-y)/max(1,abs(x),abs(y))
            measured.append(x.parent().precision() if error==0 else max(0,int(floor(-error.log2()))))
        assert min(measured)>=r['agreement_bits']
        assert r['usable_bits']<=min(measured)-15


def _write_tex(curves,output_dir,prefix):
    lines=[r'\documentclass[10pt]{article}',
        r'\usepackage[a4paper,landscape,margin=15mm]{geometry}',
        r'\usepackage{amsmath,amssymb,longtable,array}',
        r'\setlength{\tabcolsep}{4pt}\renewcommand{\arraystretch}{1.16}',
        r'\begin{document}\begin{center}\Large CM divisor images on optimal elliptic curves\end{center}',
        r'All admissible negative fundamental discriminants $|D|<200$, all ideal classes and local orientations, and all three nonidentity Atkin--Lehner operators. ',
        r'Exact points and field inclusions are checked; the analytic identification is numerical, not interval-certified.']
    for i,c in enumerate(curves):
        if i:lines.append(r'\clearpage')
        label=c['E'];display=label[:-2]+label[-2].upper()+'_{'+label[-1]+'}'
        lines.append(r'\section*{$X(%s,1)\longrightarrow %s$}'%(c['N'],display))
        lines.append(r'Base discriminant $%s$; %s and %s terms at %s bits. '%(c['base_discriminant'],*c['truncations'],c['prec_bits']))
        lines.append(r'Native quotient %s; degree-%s isogeny to this optimal model. '%(c['native_quotient'],c['isogeny_degree']))
        lines.append(r'\begin{longtable}{|r|r|r|c|>{\raggedright\arraybackslash}p{205mm}|}\hline $E$&$D$&$h$&$\tau$&Minimal polynomial of $x(P(D,\tau))$\\\hline\endhead')
        previous=None
        for j,r in enumerate(c['rows']):
            new=r['D']!=previous
            cols=[('$'+str(v)+'$') if v!='' else '' for v in (display if j==0 else '',r['D'] if new else '',r['h'] if new else '')]
            if r['status']=='identity':value=r'\text{point at infinity}'
            elif 'rational_numerator' in r:value='x='+latex(QQ(r['rational_numerator'])/r['rational_denominator'])
            else:value=_twelve_tex_polynomial(r['polynomial_coefficients'])
            if r['annotation']:value+=r'\quad(\text{'+r['annotation']+'})'
            lines.append('&'.join(cols)+r'&$\tau_{%s}$&$%s$\\\hline'%(r['tau'],value))
            previous=r['D']
        lines.append(r'\end{longtable}')
    lines.append(r'\end{document}')
    (output_dir/(prefix+'.tex')).write_text('\n'.join(lines)+'\n',encoding='utf-8')


def export_table(labels=CM_LABELS,output_dir=None,validate=True,check_reference=True):
    """Write one CSV/JSON/LaTeX table from newly computed result files.

    A complete twelve-target export is compared to the canonical fingerprint
    of the previously computed720-row table. No old table/results are loaded.
    Exact Sage points and all field embeddings remain in the result SOBJ files.
    """
    labels=(labels,) if isinstance(labels,str) else tuple(labels)
    if not labels or len(labels)!=len(set(labels)) or not set(labels)<=set(CM_LABELS):
        raise ValueError('Choose one or more distinct target labels from CM_LABELS.')
    destination=cm_output_dir(output_dir)
    coverage=expected_coverage();curves=[];compact=[];individuals=[];checks=[]
    for label in labels:
        N=int(label[:-2]);r=load_result(label,destination);meta=r['metadata']
        if validate:
            check=validate_result(label,r,coverage[str(N)])
            _check_truncations(r,destination);checks.append(check)
        rowmap={_twelve_key(a):a for a in r['rows']}
        exactmap={_twelve_key(a):a for a in r['exact_rows']};groups={}
        for row in r['rows']:
            if row['status'] not in ('recognized','identity'):
                raise ValueError('Unresolved entry: '+str((label,row['D'],row['operator'],row['status'])))
            groups.setdefault((int(row['D']),row['operator']),[]).append(row)
        crows=[]
        for (D,w),group in sorted(groups.items(),key=lambda kv:(abs(kv[0][0]),int(kv[0][1][1:]))):
            assert len({(a['status'],str(a.get('polynomial'))) for a in group})==1
            source=group[0];f=source.get('polynomial');status=source['status']
            row=dict(E=label,N=N,D=D,h=int(source['h']),operator=w,tau=int(w[1:]),status=status,
                orientation_count=len(group),polynomial='' if f is None else str(f),
                polynomial_coefficients=[] if f is None else [int(c) for c in f.list()],
                annotation=source.get('annotation',''),
                point_field_degrees=sorted({int(d) for a in group for d in a.get('point_field_degrees',[])}),
                min_xy_match_bits=min(int(a['xy_match_bits']) for a in group) if status=='recognized' else None)
            if f is not None and f.degree()==1:
                q=-QQ(f[0])/QQ(f[1]);row.update(rational_numerator=int(q.numerator()),rational_denominator=int(q.denominator()))
            crows.append(row);compact.append(row)
        cindividuals=[]
        for item in r['numeric_rows']:
            key=_twelve_key(item);row=rowmap[key]
            out=dict(E=label,D=int(item['D']),h=int(item['h']),operator=item['operator'],
                orientation_index=int(item['orientation_index']),class_index=int(item['class_index']),
                orbit_index=int(item['orbit_index']),status=row['status'],
                x_minpoly='' if row['polynomial'] is None else str(row['polynomial']),
                independent_agreement_bits=int(item.get('agreement_bits',0)))
            if row['status']=='recognized':
                entry=exactmap[key]
                match=next(m for m in entry['numeric_matches'] if m['orbit_index']==item['orbit_index'])
                point=entry['points'][int(match['point_index'])-1]
                out.update(x_field_degree=int(point['x_field']['degree']),y_field_degree=int(point['y_field']['degree']),
                    point_field_degree=int(point['point_field']['degree']),
                    point_field_polynomial=str(point['point_field']['defining_polynomial']),
                    exact_point=str(point['point']),xy_match_bits=int(match['xy_match_bits']))
            else:out.update(point_field_degree=1,exact_point='(0 : 1 : 0)')
            cindividuals.append(out);individuals.append(out)
        origin=sum(a['status']=='identity' for a in cindividuals)
        native,degree=CM_NATIVE_TARGETS[label]
        curves.append(dict(N=N,E=label,source_genus=2 if label in CM_GENUS_TWO else 1,
            ainvs=[int(a) for a in EllipticCurve(label).ainvs()],base_discriminant=int(meta['base_discriminant']),
            truncations=[int(meta['Nper_lo']),int(meta['Nper_hi'])],prec_bits=int(meta['prec_bits']),
            evaluation_method=meta.get('evaluation_method','endpoint_reduction'),
            native_quotient=native,isogeny_degree=degree,
            discriminants=[a['D'] for a in coverage[str(N)]],rows=crows,
            individual_divisors=len(cindividuals),finite_divisors=len(cindividuals)-origin,
            identity_divisors=origin,orientation_rows=len(r['rows']),compact_rows=len(crows),
            unresolved=0,min_xy_match_bits=min(a['min_xy_match_bits'] for a in crows if a['min_xy_match_bits'] is not None)))
    complete=set(labels)==set(CM_LABELS) and len(labels)==len(CM_LABELS)
    summary=dict(curves=len(curves),compact_rows=len(compact),individual_divisors=len(individuals),
        finite_divisors=sum(c['finite_divisors'] for c in curves),identity_divisors=sum(c['identity_divisors'] for c in curves),
        min_xy_match_bits=min(c['min_xy_match_bits'] for c in curves))
    digest=table_fingerprint(compact)
    if complete:
        assert (summary['individual_divisors'],summary['compact_rows'])==(5616,720)
        if check_reference:
            assert digest==EXPECTED_TABLE_SHA256,'The regenerated table differs from the reference; inspect before accepting.'
    payload=dict(scope='all admissible negative fundamental discriminants with |D|<200; no class-number cutoff',
        target_convention='Sage X_0(N)-optimal targets',complete=complete,curves=curves,summary=summary,
        validated=bool(validate),reference_table_sha256=digest,
        reference_table_matches=bool(complete and digest==EXPECTED_TABLE_SHA256),
        proof_boundary='Exact curve equations and field inclusions; analytic identification numerical, not interval-certified.')
    prefix='cm_divisors_twelve_under200' if complete else 'cm_divisors_subset'
    destination.mkdir(parents=True,exist_ok=True)
    (destination/(prefix+'_data.json')).write_text(json.dumps(payload,indent=2,default=int)+'\n',encoding='utf-8')
    _twelve_write_csv(destination/(prefix+'.csv'),compact,
        ['E','D','h','operator','status','polynomial','annotation','point_field_degrees','orientation_count','min_xy_match_bits'])
    _twelve_write_csv(destination/(prefix+'_all_points.csv'),individuals,
        ['E','D','h','operator','orientation_index','class_index','orbit_index','status','x_minpoly',
         'x_field_degree','y_field_degree','point_field_degree','point_field_polynomial','exact_point','xy_match_bits','independent_agreement_bits'])
    _write_tex(curves,destination,prefix)
    if validate:
        (destination/(prefix+'_validation.json')).write_text(json.dumps(dict(summary=summary,curves=checks,
            reference_table_matches=payload['reference_table_matches']),indent=2,default=int)+'\n',encoding='utf-8')
    print(summary,flush=True)
    return payload
