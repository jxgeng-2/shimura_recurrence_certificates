"""Fresh, resumable computation for all fundamental |D|<200 on six factors.

Run from this directory. Inputs are only the adjacent source/geometry files,
../cm_data exports, and ../common helpers. No numerical or recognized results
are supplied or required. Generated files go exclusively into out_dir.

The high-order guarded series is generated once; its prefix supplies the lower
independent truncation. Recognition still requires separated algdep windows,
30 held-out bits, exact H_D containment, both-coordinate matching and descent.
"""
from sage.all import *
from sage.repl.load import load as _g2_source_load
from pathlib import Path
import hashlib,json,os,time

G2_DIR=Path.cwd().resolve()
for _g2_module in ('../common/heegner_image_from_series.sage',
                  '../common/cm_image_table_package.sage','periods.sage',
                  'x26_cm11_setup.sage','x38_cm11_setup.sage','x58_cm11_setup.sage'):
    load(_g2_module)

G2_PRESETS={'26a1':(1000,1200,1000),'26b1':(3000,3200,1600),
           '38a1':(1600,1800,1000),'38b1':(3000,3200,1000),
           '58a1':(4000,4400,1000),'58b1':(4000,4400,1000)}

def _g2_save(value,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=str(path)+'.pending.sobj';save(value,temporary);os.replace(temporary,path)

def _g2_progress(out,label,stage,**values):
    path=out/(label+'_progress.json');temporary=path.with_name(path.name+'.pending')
    record=dict(E=label,stage=stage,updated_unix=time.time(),**values)
    temporary.write_text(json.dumps(record,indent=2,default=str),encoding='utf-8');os.replace(temporary,path)

def _g2_data(N,check_coverage=True):
    ns=dict(globals());_g2_source_load(str(G2_DIR.parent/'cm_data'/('x%s_cm_data.sage'%N)),ns)
    data={k:ns[k] for k in ('CM_DB','CM_RAMIFIED_PRIMES','CM_ABS_D_MIN','CM_ABS_D_MAX',
        'CM_ALGEBRA_I2','CM_ALGEBRA_J2','CM_ORDER_BASIS','CM_BASE_COORDS','CM_AL_COORDS')}
    data['points']=[dict(r) for r in ns['cm_points_all_orientations']]
    assert data['CM_DB']==N and data['CM_ABS_D_MIN']==3 and data['CM_ABS_D_MAX']==199
    if check_coverage:
        primes=ZZ(N).prime_divisors()
        expected=[-d for d in range(3,200) if fundamental_discriminant(-d)==-d and all(kronecker(-d,p)!=1 for p in primes)]
        assert set(expected)=={r['D'] for r in data['points']}
        for D in expected:
            rows=[r for r in data['points'] if r['D']==D];h=QuadraticField(D,'s').class_number()
            assert all(r['h']==h for r in rows)
            assert len(rows)==h*prod(1-kronecker(D,p) for p in primes)
            assert len({r['orbit_index'] for r in rows})==len(rows)
            orientations=prod(1-kronecker(D,p) for p in primes)
            assert {r['orientation_index'] for r in rows}==set(range(1,int(orientations)+1))
            for orientation in range(1,int(orientations)+1):
                subset=[r for r in rows if r['orientation_index']==orientation]
                assert len(subset)==h and {r['class_index'] for r in subset}==set(range(1,int(h)+1))
    return data

def _g2_fingerprint(label,spec):
    N=int(label[:-2]);digest=hashlib.sha256(repr((label,tuple(map(int,spec)))).encode())
    paths=[G2_DIR/'run_genus_two.sage',G2_DIR/'periods.sage',G2_DIR/('x%s_cm11_setup.sage'%N),
           G2_DIR/('x%s_cm11_side_data.json'%N),G2_DIR.parent/'cm_data'/('x%s_cm_data.sage'%N),
           G2_DIR.parent/'common'/'heegner_image_from_series.sage',G2_DIR.parent/'common'/'cm_image_table_package.sage']
    for path in paths:digest.update(path.read_bytes())
    return digest.hexdigest()

def _g2_agreement(a,b):
    err=abs(a-b)/max(1,abs(a),abs(b))
    return min(a.parent().precision(),b.parent().precision()) if err==0 else max(0,int(floor(-err.log2())))

def _g2_setup_pair(label,lo,hi,bits,out,stamp):
    N=int(label[:-2]);component=label[-2]
    path=out/('%s_setup_%s_%s.sobj'%(label,hi,bits))
    high=load(str(path)) if path.exists() else None
    if high is None or high.get('source_fingerprint')!=stamp:
        _g2_progress(out,label,'generating_series',Nper=hi,prec_bits=bits)
        high=globals()['build_x%s_cm11_setup'%N](label[:-1] if N==26 else component,hi,bits)
        high['source_fingerprint']=stamp;_g2_save(high,path)
    low=dict(high,b=list(high['b'][:lo]),Nper=lo,derived_from_Nper=hi)
    _g2_save(low,out/('%s_setup_%s_%s.sobj'%(label,lo,bits)))
    assert len(low['b'])==lo and high['b'][:lo]==low['b']
    return low,high

def _g2_geometry_and_kernel(data):
    N=int(data['CM_DB']);B=QuaternionAlgebra(QQ,data['CM_ALGEBRA_I2'],data['CM_ALGEBRA_J2'])
    assert B.discriminant()==N
    inverse=matrix(QQ,data['CM_ORDER_BASIS']).inverse()
    order_coords=lambda q:vector(QQ,q)*inverse
    AL={k:B(c) for k,c in data['CM_AL_COORDS'].items()}
    for name,w in AL.items():
        assert w.reduced_norm()>0
        for row in data['CM_ORDER_BASIS']:
            assert all(c in ZZ for c in order_coords(w*B(row)*(~w)))
        square=w*w/w.reduced_norm()
        assert square.reduced_norm()==1 and all(c in ZZ for c in order_coords(square))
    for record in data['points']:
        q=B(record['mu_coords'])
        assert q.reduced_trace()**2-4*q.reduced_norm()==record['D']
    v=AL['w2'];commutators={}
    for name,w in AL.items():
        g=w*v*(~(v*w));coordinates=order_coords(g)
        assert g.reduced_norm()==1 and all(c in ZZ for c in coordinates)
        commutators[name]=dict(quaternion_coords=list(map(str,g)),order_coords=list(map(str,coordinates)),norm=1)
    return dict(kernel='w2',normalizes_order=True,commutators=commutators)

def _g2_al_eigenvalues(label,setup,AL,reducer):
    """Fresh checks of f|W, at two generic points; no input audit JSON."""
    N=int(label[:-2]);component=label[-2];tau=setup['tau0'];CC=tau.parent();bits=CC.precision();RF=RealField(bits)
    kernel_a={26:'w13',38:'w19',58:'w58'}[N]
    def value(z):
        w=(z-tau)/(z-tau.conjugate());s=CC(0)
        for c in reversed(setup['b']):s=s*w+c
        return (1-w)**2*s
    checks=[]
    for name,W in AL.items():
        expected=1 if name==(kernel_a if component=='a' else 'w2') else -1
        scores=[]
        for w0 in (CC(0),CC(RF(1)/10,RF(7)/100)):
            z=(tau-tau.conjugate()*w0)/(1-w0);wz=_mobius(W,z)
            gamma,wr=reducer.reduce_upper(wz);zr=(tau-tau.conjugate()*wr)/(1-wr)
            eigen=W.det()*value(zr)/((W[1,0]*z+W[1,1])**2*(gamma[1,0]*wz+gamma[1,1])**2*value(z))
            err=abs(eigen-expected);score=bits if err==0 else int(floor(-err.log2()))
            assert score>=60,(label,name,'AL check requires more terms',score)
            scores.append(score)
        checks.append(dict(operator=name,expected_eigenvalue=expected,agreement_bits=scores))
    return checks

def _g2_context(label,setup,data,out,stamp):
    N=int(label[:-2]);component=label[-2];nper,bits=int(setup['Nper']),int(setup['prec_bits'])
    path=out/('%s_periods_%s_%s.sobj'%(label,nper,bits))
    cached=load(str(path)) if path.exists() else None
    if cached is None or cached.get('source_fingerprint')!=stamp:
        audit=globals()['audit_x%s_series_modularity'%N](setup)
        periods=g2_side_periods(setup,audit)
        ordered=sorted(periods,key=lambda r:r['stats']['tail_abs_estimate'])
        values=[r['value'] for r in ordered]
        matched=g2_match_periods(values,label,bits,min_match_bits=60)
        native=label if component=='a' else label[:-1]+'2';degree=1 if component=='a' else (7 if N==26 else 5)
        E,J=EllipticCurve(label),EllipticCurve(native)
        assert E==E.optimal_curve() and J.isogeny_degree(E)==degree and matched['index']==degree
        native_match=matched if component=='a' else g2_match_periods(values,native,bits,min_match_bits=60)
        assert native_match['index']==1
        assert abs(matched['mu']/native_match['mu']-degree)<RealField(bits)(2)**(-60)
        cached=dict(source_fingerprint=stamp,modularity_audit=audit,period_match=matched,
                    native_period_match=native_match,native_quotient=native,isogeny_degree=degree)
        _g2_save(cached,path)
    RF,CC=RealField(bits),ComplexField(bits)
    AL={k:globals()['x%s_quaternion_matrix'%N](v,RF) for k,v in data['CM_AL_COORDS'].items()}
    reducer=DirichletReducer(setup['side_data'],prec_bits=bits)
    cm_center=_heegner_fixed_point(globals()['x%s_quaternion_matrix'%N](data['CM_BASE_COORDS'],RF),CC)
    assert abs(cm_center-setup['tau0'])<RF(2)**(-min(100,bits//2))
    return cached,AL,reducer

def _g2_numeric(label,setup,data,out,stamp):
    N=int(label[:-2]);nper,bits=int(setup['Nper']),int(setup['prec_bits'])
    context,AL,reducer=_g2_context(label,setup,data,out,stamp)
    path=out/('%s_numeric_%s_%s.sobj'%(label,nper,bits));payload=load(str(path)) if path.exists() else None
    if payload is None or payload.get('source_fingerprint')!=stamp:
        payload=dict(rows=[],source_fingerprint=stamp,Nper=nper,prec_bits=bits,E_label=label,
            period_match=context['period_match'],modularity_audit=context['modularity_audit'],
            basepoint=setup['tau0'],AL_coords=data['CM_AL_COORDS'])
    done={(r['D'],r['orbit_index'],r['operator']) for r in payload['rows']}
    RF,CC=RealField(bits),ComplexField(bits)
    for record in data['points']:
        tauD=_heegner_fixed_point(globals()['x%s_quaternion_matrix'%N](record['mu_coords'],RF),CC)
        for operator in sorted(AL):
            if (record['D'],record['orbit_index'],operator) in done:continue
            ans=heegner_image_from_series(operator,tauD,b=setup['b'],tau0=setup['tau0'],AL=AL,
                mu=context['period_match']['mu'],E_label=label,Nper=nper,prec_bits=bits,reducer=reducer,
                verbose=False,identity_min_bits=None,max_radius=.985)
            row=dict(record,operator=operator,tauD=tauD,zE=ans['zE'],zE_reduced=ans['zE_reduced'],
                x=ans['P'][0],y=ans['P'][1],stats=ans['stats'],identity_bits=ans['identity_bits'],is_infinity=False)
            payload['rows'].append(row);_g2_save(payload,path)
            _g2_progress(out,label,'numerical_images',Nper=nper,done=len(payload['rows']),total=3*len(data['points']))
    return payload,path,context,AL,reducer

def _g2_precondition(label,setup,payload,other,path,AL,reducer,kernel_checks):
    if label[-2]!='b':return
    bits=int(setup['prec_bits']);RF=RealField(bits)
    lookup={(r['D'],r['orbit_index'],r['operator']):r for r in other['rows']}
    radius=lambda z:abs(reducer.reduce_upper(z)[1])
    for r in payload['rows']:
        if r.get('kernel_preconditioning',{}).get('version')==1:continue
        z=r['tauD'];W=AL[r['operator']];vz=_mobius(AL['w2'],z)
        r0=max(radius(z),radius(_mobius(W,z)));r1=max(radius(vz),radius(_mobius(W,vz)))
        info=dict(version=1,kernel='w2',selected='identity',original_tau=z,evaluation_tau=z,
                  original_radius=r0,alternate_radius=r1)
        if r['operator']!='w2' and r1+RF('0.001')<r0 and r['identity_bits']<70:
            original={k:r[k] for k in ('x','y','zE','zE_reduced','stats','identity_bits')}
            other_row=lookup[r['D'],r['orbit_index'],r['operator']]
            other_original=other_row.get('kernel_preconditioning',{}).get('original',other_row)
            known=min(_g2_agreement(original[k],other_original[k]) for k in ('x','y'))
            ans=heegner_image_from_series(r['operator'],vz,b=setup['b'],tau0=setup['tau0'],AL=AL,
                mu=payload['period_match']['mu'],E_label=label,Nper=setup['Nper'],prec_bits=bits,
                reducer=reducer,verbose=False,identity_min_bits=None,max_radius=.985)
            scores=[_g2_agreement(original[k],ans['P'][j]) for j,k in enumerate(('x','y'))]
            assert min(scores)>=max(60,min(known,bits-40)-20),(label,r['D'],scores,known)
            info.update(selected='w2',evaluation_tau=vz,original=original,
                original_two_truncation_bits=known,original_xy_agreement_bits=scores)
            r.update(x=ans['P'][0],y=ans['P'][1],zE=ans['zE'],zE_reduced=ans['zE_reduced'],
                     stats=ans['stats'],identity_bits=ans['identity_bits'])
        r['kernel_preconditioning']=info;payload['kernel_validation']=kernel_checks;_g2_save(payload,path)

def _g2_match_exact(exact,items,bits,required):
    H=exact['ambient_class_field'];embedding=H.embeddings(ComplexField(bits))[0]
    pairs=[tuple(embedding(p['ambient_point'][j]) for j in (0,1)) for p in exact['points']]
    matches=[]
    for item in items:
        scores=[min(_g2_agreement(x,item['x']),_g2_agreement(y,item['y'])) for x,y in pairs]
        best=max(range(len(scores)),key=lambda j:scores[j]);assert scores[best]>=required
        matches.append(dict(orbit_index=item['orbit_index'],class_index=item['class_index'],
                            point_index=best+1,xy_match_bits=scores[best]))
    exact.update(complex_embedding=embedding,numeric_matches=matches)
    return min(m['xy_match_bits'] for m in matches)

def _g2_recognize(label,low,high,metadata,out,stamp):
    bits=int(metadata['prec_bits']);cap=bits-100;assert cap>=100
    lows={(r['D'],r['orbit_index'],r['operator']):r for r in low['rows']};groups={}
    for r in high['rows']:
        old=lows[r['D'],r['orbit_index'],r['operator']]
        r['agreement_bits']=min(_g2_agreement(r['x'],old['x']),_g2_agreement(r['y'],old['y']))
        r['is_infinity']=min(r['identity_bits'],old['identity_bits'])>=70
        r['identity_method']='period_lattice_two_truncations'
        r['usable_bits']=min(r['identity_bits'],old['identity_bits']) if r['is_infinity'] else max(0,min(r['agreement_bits']-15,cap))
        groups.setdefault((r['D'],r['orientation_index'],r['operator']),[]).append(r)
    path=out/('x%s_cm_images_%s_result.sobj'%(label[:-2],label))
    old=load(str(path)) if path.exists() else {}
    if old.get('metadata',{}).get('source_fingerprint')!=stamp:old={}
    key=lambda r:(r['D'],r['orientation_index'],r['operator'])
    rows={key(r):r for r in old.get('rows',[])};exacts={key(r):r for r in old.get('exact_rows',[])}
    fields={};R=PolynomialRing(ZZ,'x');E=EllipticCurve(label)
    for group,items in sorted(groups.items(),key=lambda kv:(abs(kv[0][0]),kv[0][1],kv[0][2])):
        if group in rows and rows[group]['status'] in ('recognized','identity'):continue
        D,orientation,operator=group;items.sort(key=lambda r:r['class_index']);h=items[0]['h']
        if not all(r['is_infinity'] for r in items) and D not in fields:fields[D]=_cm_default_hilbert_class_field(D,h)
        H=fields.get(D)
        selected=max(items,key=lambda r:(abs(r['x'].imag())>RealField(bits)(2)**(-50)*max(1,abs(r['x'])),r['usable_bits']))
        rec=_cm_recognize_x(selected,items,D,h,R,cap,tuple(range(80,cap+1,20)),10**200,50,1,15,2*h,H,None,True,True)
        rec['selected_class_index']=selected['class_index']
        row=dict(E=label,D=D,h=h,orientation_index=orientation,operator=operator,operator_label=operator[1:],
            status=rec['status'],polynomial=rec.get('polynomial'),recognition=rec,
            annotation=_cm_two_torsion_annotation(E,R,rec.get('polynomial')))
        if row['status']=='recognized':
            try:
                exact=_cm_recover_row_points(row,H);exact['orientation_index']=orientation
                row['xy_match_bits']=_g2_match_exact(exact,items,bits,max(50,int(rec['algdep_use_bits'])+30))
                row['point_field_degrees']=sorted({int(p['point_field']['degree']) for p in exact['points']});exacts[group]=exact
            except Exception as error:row['status']='exact_point_error';rec['message']=str(error)
        rows[group]=row
        result=dict(rows=list(rows.values()),exact_rows=list(exacts.values()),numeric_rows=high['rows'],metadata=metadata)
        _g2_save(result,path);_g2_progress(out,label,'recognition',done=len(rows),total=len(groups),D=int(D),status=row['status'])
        print(label,D,orientation,operator,row['status'],flush=True)
    result=dict(rows=list(rows.values()),exact_rows=list(exacts.values()),numeric_rows=high['rows'],metadata=metadata)
    _g2_save(result,path);return result

def run_genus_two(labels=tuple(G2_PRESETS),out_dir='results',pass_specs=None,require_complete=True):
    """Build all requested images from source. Resume uses only this run's files."""
    if isinstance(labels,str):labels=(labels,)
    specs=dict(G2_PRESETS);specs.update(pass_specs or {});out=Path(out_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    answer={}
    for label in labels:
        assert label in G2_PRESETS;N=int(label[:-2]);lo,hi,bits=map(int,specs[label]);assert 0<lo<hi and bits>=200
        data=_g2_data(N);stamp=_g2_fingerprint(label,(lo,hi,bits));geometry=_g2_geometry_and_kernel(data)
        low_setup,high_setup=_g2_setup_pair(label,lo,hi,bits,out,stamp)
        low,low_path,_,low_AL,low_reducer=_g2_numeric(label,low_setup,data,out,stamp)
        high,high_path,context,high_AL,high_reducer=_g2_numeric(label,high_setup,data,out,stamp)
        eigenvalues=_g2_al_eigenvalues(label,high_setup,high_AL,high_reducer)
        kernel_checks=dict(exact=geometry,eigenvalues=eigenvalues)
        _g2_precondition(label,low_setup,low,high,low_path,low_AL,low_reducer,kernel_checks)
        _g2_precondition(label,high_setup,high,low,high_path,high_AL,high_reducer,kernel_checks)
        metadata=dict(E_label=label,Nper_lo=lo,Nper_hi=hi,prec_bits=bits,base_discriminant=-11,abs_D_max=199,
            native_quotient=context['native_quotient'],isogeny_degree=int(context['isogeny_degree']),
            period_lattice_index=int(context['period_match']['index']),source_fingerprint=stamp,
            numeric_low_cache=low_path.name,numeric_high_cache=high_path.name,
            evaluation_method='endpoint_reduction_with_verified_kernel_AL_preconditioning' if label[-2]=='b' else 'endpoint_reduction',
            coefficient_strategy='guarded_high_order_series_with_shorter_truncation',
            series_working_bits=int(high_setup['working_bits']),
            analytic_identification_is_interval_certified=False)
        result=_g2_recognize(label,low,high,metadata,out,stamp);answer[label]=result
        unresolved=[(r['D'],r['operator'],r['status']) for r in result['rows'] if r['status'] not in ('recognized','identity')]
        _g2_progress(out,label,'complete' if not unresolved else 'unresolved',divisors=len(result['numeric_rows']),unresolved=unresolved)
        if require_complete and unresolved:raise RuntimeError('%s unresolved entries; checkpoints saved: %s'%(label,unresolved))
    return answer

def genus_two_sanity_check(label='26a1',Nper=20,prec_bits=80):
    """Tiny coefficient + one CM-lift test only; no periods, images or algdep."""
    assert label in G2_PRESETS and 7<=Nper<=30
    N=int(label[:-2]);data=_g2_data(N,check_coverage=False)
    setup=globals()['build_x%s_cm11_setup'%N](label[:-1] if N==26 else label[-2],Nper,prec_bits)
    assert len(setup['b'])==Nper and setup['b'][0]==1
    record=data['points'][0];RF,CC=RealField(prec_bits),ComplexField(prec_bits)
    matrix_cm=globals()['x%s_quaternion_matrix'%N](record['mu_coords'],RF)
    tau=_heegner_fixed_point(matrix_cm,CC);assert tau.imag()>0
    assert abs(_mobius(matrix_cm,tau)-tau)<RF(2)**(-prec_bits//2)
    reducer=DirichletReducer(setup['side_data'],prec_bits=prec_bits);_,w=reducer.reduce_upper(tau)
    assert abs(w)<1
    return dict(E=label,terms=Nper,prec_bits=prec_bits,CM_discriminant=int(record['D']),
                reduced_radius=float(abs(w)),full_batch_run=False,images_or_recognition_run=False)
