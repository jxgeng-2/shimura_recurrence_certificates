"""Two-truncation CM image computation and recognition for all twelve forms."""
from sage.all import *
from pathlib import Path
import hashlib,json,os,time
SIX_SOURCE_DIR=G2_DIR=REPRO_ROOT

SIX_PASS_SPECS={14:(1000,1200,1000),15:(1000,1200,1000),21:(2400,2600,1000),33:(4000,4400,1000),34:(2400,2800,1000),46:(4000,4400,1000)}
SIX_BASE_DISCRIMINANTS={14:-11,15:-7,21:-43,33:-67,34:-11,46:-163}
SIX_NATIVE_QUOTIENTS={14:('14a2',2),15:('15a1',1),21:('21a2',2),33:('33a1',1),34:('34a3',3),46:('46a2',2)}

def _six_save(value,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    pending=path.with_name(path.stem+'.pending.sobj')
    save(value,str(pending));os.replace(str(pending),str(path))

def _six_key(row):return (int(row['D']),int(row['orbit_index']),row['operator'])
def _six_group_key(row):return (int(row['D']),int(row['orientation_index']),row['operator'])
def _six_agree(a,b):
    bits=min(a.parent().precision(),b.parent().precision())
    error=abs(a-b)/max(1,abs(a),abs(b))
    return int(bits) if error==0 else max(0,int(floor(-error.log2())))
def _six_mob(g,z):return (g[0,0]*z+g[0,1])/(g[1,0]*z+g[1,1])

def six_load_inputs(level):
    """Load exact geometry and recurrence definitions only; no series generation."""
    N=int(level);disc=abs(SIX_BASE_DISCRIMINANTS[N]);stem='x%s_cm%s'%(N,disc)
    ns=dict(globals());ns.update(cm_inputs(N))
    assert int(ns['CM_DB'])==N and int(ns['CM_ABS_D_MIN'])==3 and int(ns['CM_ABS_D_MAX'])==199
    records=[dict(r) for r in ns['cm_points_all_orientations']]
    AL={k:list(v) for k,v in ns['CM_AL_COORDS'].items()}
    assert len({(int(r['D']),int(r['orbit_index'])) for r in records})==len(records)
    assert all(3<=abs(int(r['D']))<200 and fundamental_discriminant(r['D'])==r['D'] for r in records)
    assert {int(r['D']) for r in records}==set(map(int,ns['CM_SELECTED_DS']))
    geometry=dict(N=N,algebra=(ns['CM_ALGEBRA_I2'],ns['CM_ALGEBRA_J2']),
        order_basis=[list(r) for r in ns['CM_ORDER_BASIS']],AL_coords=AL,records=records)
    signature=hashlib.sha256(repr(geometry).encode()).hexdigest()
    return dict(geometry,signature=signature,namespace=ns,stem=stem,
        side_path=side_data(N),
        builder=ns['build_'+stem+'_setup'],matrix=ns['x%s_quaternion_matrix'%N],
        modularity=ns['audit_x%s_series_modularity'%N])

def _six_setup(inputs,nper,bits,output):
    path=output/('x%s_setup_%s_%s.sobj'%(inputs['N'],nper,bits))
    source=(REPRO_ROOT/'Fundamental Domain Side Data.json').read_bytes()+(REPRO_ROOT/'Compute CM Expansion Series.sage').read_bytes()
    signature=hashlib.sha256(source).hexdigest()
    if path.exists():
        setup=load(str(path))
        assert setup.get('source_signature')==signature, 'Setup source changed; use a fresh output directory.'
    else:
        print('Building guarded series',inputs['N'],nper,bits,flush=True)
        setup=inputs['builder'](inputs['side_path'],nper,bits)
        setup['source_signature']=signature;_six_save(setup,path)
    assert int(setup['Nper'])==nper and int(setup['prec_bits'])==bits
    return setup

def _six_periods(inputs,setup,reducer,output):
    N,nper,bits=inputs['N'],int(setup['Nper']),int(setup['prec_bits'])
    path=output/('x%s_periods_%s_%s.sobj'%(N,nper,bits))
    if path.exists():
        saved=load(str(path));assert saved['source_signature']==setup['source_signature'];return saved['rows']
    tau=setup['tau0'];rows=[]
    if N==21:
        # Preserve the tested X21 continued-path normalization.
        for i,g in enumerate(setup['S_matrices']):
            value,_,_,stats=integrate_between_points_from_series(tau,_six_mob(g,tau),setup['b'],tau,
                Nper=nper,prec_bits=bits,reducer=reducer,initial_subdivisions=32,max_radius=.91,verbose=False)
            rows.append(dict(side_index=i,value=value,stats=stats))
    else:
        H=antiderivative_evaluator(setup['b'],tau,nper,bits);bar=tau.conjugate()
        for i,g in enumerate(setup['S_matrices']):
            endpoint=(_six_mob(g.inverse(),tau)-tau)/(_six_mob(g.inverse(),tau)-bar)
            w=endpoint/(1+(1-abs(endpoint)**2).sqrt());z=(tau-bar*w)/(1-w)
            gz=_six_mob(g,z);gw=(gz-tau)/(gz-bar)
            stats=dict(mode='balanced_period_cocycle',max_u=max(abs(w),abs(gw)),
                tail_abs_estimate=H.tail_estimate(w)+H.tail_estimate(gw),tail_is_rigorous=False,segments=1,bad_segments=0)
            rows.append(dict(side_index=i,value=H(gw)-H(w),stats=stats))
    _six_save(dict(rows=rows,source_signature=setup['source_signature']),path)
    return rows

def six_numeric_run(level,nper,bits=1000,output_dir=None):
    """All classes/orientations and all three operators; checkpoint after each image."""
    N,nper,bits=int(level),int(nper),int(bits);inputs=six_load_inputs(N)
    output=Path(output_dir).resolve() if output_dir is not None else SIX_SOURCE_DIR/'results'
    setup=_six_setup(inputs,nper,bits,output);CC,RF=ComplexField(bits),RealField(bits)
    reducer=DirichletReducer(inputs['side_path'],prec_bits=bits)
    audit=inputs['modularity'](setup)
    assert min(r['residual_bits'] for r in audit)>=60, 'Insufficient balanced modularity accuracy.'
    period_rows=_six_periods(inputs,setup,reducer,output)
    ordered=period_rows if N==21 else sorted(period_rows,key=lambda r:r['stats']['tail_abs_estimate'])
    matched=six_match_periods([r['value'] for r in ordered],E_label=str(N)+'a1',prec_bits=bits,min_match_bits=60)
    assert int(matched['index'])==SIX_NATIVE_QUOTIENTS[N][1], 'Unexpected optimal-target period index.'
    AL={k:inputs['matrix'](v,RF) for k,v in inputs['AL_coords'].items()}
    path=output/('x%s_numeric_%s_%s.sobj'%(N,nper,bits));payload=load(str(path)) if path.exists() else {}
    if payload:
        assert payload['geometry_signature']==inputs['signature'] and payload['period_match']['mu']==matched['mu']
    rows=payload.get('rows',[]);done={_six_key(r) for r in rows}
    exact={(int(r['D']),int(r['orbit_index'])):r for r in inputs['records']}
    for r in rows:assert all(r[k]==v for k,v in exact[(int(r['D']),int(r['orbit_index']))].items())
    payload.update(rows=rows,Nper=nper,prec_bits=bits,E_label=str(N)+'a1',period_match=matched,
        modularity_audit=audit,basepoint=setup['tau0'],AL_coords=inputs['AL_coords'],geometry_signature=inputs['signature'])
    for record in inputs['records']:
        tauD=_cm_fixed_point_upper_half_plane(inputs['matrix'](record['mu_coords'],RF),CC)
        for operator in sorted(AL):
            if (int(record['D']),int(record['orbit_index']),operator) in done:continue
            ans=heegner_image_from_series(operator,tauD,b=setup['b'],tau0=setup['tau0'],AL=AL,
                mu=matched['mu'],E_label=str(N)+'a1',Nper=nper,prec_bits=bits,reducer=reducer,
                verbose=False,identity_min_bits=None,max_radius=.985)
            row=dict(record);row.update(operator=operator,tauD=tauD,zE=ans['zE'],zE_reduced=ans['zE_reduced'],
                x=ans['P'][0],y=ans['P'][1],stats=ans['stats'],identity_bits=ans['identity_bits'],is_infinity=False)
            rows.append(row);_six_save(payload,path)
            print(N,nper,len(rows),'/',3*len(inputs['records']),record['D'],record['orbit_index'],operator,flush=True)
    assert len(rows)==3*len(inputs['records'])
    _six_save(payload,path);return payload

def _six_match_exact(exact,items,bits,required):
    H=exact['ambient_class_field'];embedding=H.embeddings(ComplexField(bits))[0]
    pairs=[(embedding(p['ambient_point'][0]),embedding(p['ambient_point'][1])) for p in exact['points']]
    matches=[]
    for row in items:
        agrees=[min(_six_agree(x,row['x']),_six_agree(y,row['y'])) for x,y in pairs]
        best=max(range(len(pairs)),key=lambda k:agrees[k])
        if agrees[best]<required:raise ValueError('Both coordinates fail exact-point matching: %s < %s'%(agrees[best],required))
        matches.append(dict(orbit_index=row['orbit_index'],class_index=row['class_index'],point_index=best+1,xy_match_bits=agrees[best]))
    exact['complex_embedding']=embedding;exact['numeric_matches']=matches
    return min(m['xy_match_bits'] for m in matches)

def six_recognize(level,low_nper,high_nper,bits=1000,output_dir=None):
    N,bits=int(level),int(bits);inputs=six_load_inputs(N);label=str(N)+'a1'
    output=Path(output_dir).resolve() if output_dir is not None else SIX_SOURCE_DIR/'results'
    low_name='x%s_numeric_%s_%s.sobj'%(N,low_nper,bits);high_name='x%s_numeric_%s_%s.sobj'%(N,high_nper,bits)
    low,high=load(str(output/low_name)),load(str(output/high_name));low_by={_six_key(r):r for r in low['rows']}
    expected={(int(r['D']),int(r['orbit_index']),op) for r in inputs['records'] for op in inputs['AL_coords']}
    assert set(low_by)=={_six_key(r) for r in high['rows']}==expected, 'Both complete numerical passes are required.'
    assert low['geometry_signature']==high['geometry_signature']==inputs['signature']
    groups={}
    for original in high['rows']:
        r=dict(original);old=low_by[_six_key(r)]
        agrees=min(_six_agree(r['x'],old['x']),_six_agree(r['y'],old['y']))
        r['agreement_bits']=agrees;r['is_infinity']=min(r['identity_bits'],old['identity_bits'])>=70
        r['identity_method']=r.get('identity_method','period_lattice_two_truncations')
        r['usable_bits']=max(0,min(agrees-15,900,r.get('accuracy_cap_bits',900),old.get('accuracy_cap_bits',900))) if not r['is_infinity'] else min(r['identity_bits'],old['identity_bits'])
        groups.setdefault(_six_group_key(r),[]).append(r)
    method='w11_preconditioning_and_affine_w33_translation' if N==33 else 'endpoint_reduction'
    if N==33:assert 'affine_AL_translation_audit' in high and 'affine_AL_translation_audit' in low
    native,degree=SIX_NATIVE_QUOTIENTS[N]
    metadata=dict(E_label=label,Nper_lo=int(low_nper),Nper_hi=int(high_nper),prec_bits=bits,
        base_discriminant=SIX_BASE_DISCRIMINANTS[N],discriminant_abs_bound=200,abs_D_max=199,
        native_quotient=native,native_label=native,isogeny_degree=int(degree),period_lattice_index=int(high['period_match']['index']),
        evaluation_method=method,numeric_low_cache=low_name,numeric_high_cache=high_name,geometry_signature=inputs['signature'])
    result_path=output/('x%s_cm_images_%s_result.sobj'%(N,label));prior=load(str(result_path)) if result_path.exists() else {}
    valid=prior.get('metadata')==metadata
    prior_rows={_six_group_key(r):r for r in prior.get('rows',[]) if r['status'] in ('recognized','identity')} if valid else {}
    prior_exact={_six_group_key(r):r for r in prior.get('exact_rows',[])} if valid else {}
    rows=[];exact_rows=[];fields={};R=PolynomialRing(ZZ,'x');numeric_rows=[r for items in groups.values() for r in items]
    for key,items in sorted(groups.items(),key=lambda kv:(abs(kv[0][0]),kv[0][1],kv[0][2])):
        D,orientation,operator=key;items.sort(key=lambda r:r['class_index'])
        if key in prior_rows:
            rows.append(prior_rows[key])
            if key in prior_exact:exact_rows.append(prior_exact[key])
            continue
        h=items[0]['h']
        if D not in fields:fields[D]=_cm_default_hilbert_class_field(D,h)
        H=fields[D];selected=max(items,key=lambda r:(abs(r['x'].imag())>RealField(bits)(2)**(-50)*max(1,abs(r['x'])),r['usable_bits']))
        rec=_cm_recognize_x(selected,items,D,h,R,900,
            (80,100,120,160,200,220,240,260,280,300,320,340,360,380,400,450,500,550,600,650,700,750,800,850,900),
            10**100,50,1,15,2*h,H,None,True,True)
        rec['selected_class_index']=selected['class_index']
        row=dict(E=label,D=D,h=h,orientation_index=orientation,operator=operator,operator_label=operator[1:],
            status=rec['status'],polynomial=rec.get('polynomial'),recognition=rec,
            annotation=_cm_two_torsion_annotation(EllipticCurve(label),R,rec.get('polynomial')))
        if row['status']=='recognized':
            try:
                reusable=next((e for e in exact_rows if e['D']==D and e['x_polynomial']==row['polynomial']),None)
                exact=dict(reusable) if reusable is not None else _cm_recover_row_points(row,H)
                exact.update(operator=operator,orientation_index=orientation)
                row['xy_match_bits']=_six_match_exact(exact,items,bits,max(50,int(rec['algdep_use_bits'])+30))
                row['point_field_degrees']=sorted({int(p['point_field']['degree']) for p in exact['points']})
                exact_rows.append(exact)
            except Exception as err:row['status']='exact_point_error';rec['message']=str(err)
        rows.append(row);print('Recognition',N,D,orientation,operator,row['status'],row['polynomial'],flush=True)
        _six_save(dict(rows=rows,exact_rows=exact_rows,numeric_rows=numeric_rows,metadata=metadata),result_path)
    result=dict(rows=rows,exact_rows=exact_rows,numeric_rows=numeric_rows,metadata=metadata)
    _six_save(result,result_path);return result

def run_six_curves(levels=(14,15,21,33,34,46),pass_specs=None,output_dir=None):
    """Run requested levels/labels from source, resuming only self-generated checkpoints."""
    if isinstance(levels,(int,Integer,str)):levels=(levels,)
    output=Path(output_dir).resolve() if output_dir is not None else SIX_SOURCE_DIR/'results'
    output.mkdir(parents=True,exist_ok=True);results={}
    for level in levels:
        N=int(str(level).split('a')[0]) if isinstance(level,str) else int(level)
        if isinstance(level,str) and level not in (str(N),str(N)+'a1'):
            raise ValueError('Use the optimal target label %sa1.'%N)
        lo,hi,bits=map(int,(pass_specs or {}).get(N,SIX_PASS_SPECS[N]))
        assert lo<hi
        for nper in (lo,hi):six_numeric_run(N,nper,bits,output)
        if N==33:
            six33_prepare_pair(lo,hi,bits,output)
        results[N]=six_recognize(N,lo,hi,bits,output)
        require_complete_result(str(N)+'a1',results[N])
    return results


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
    ns=cm_inputs(N)
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
    paths=[REPRO_ROOT/name for name in ('Compute and Validate CM Images.sage','Match Periods and Normalizations.sage','Compute CM Expansion Series.sage','Fundamental Domain Side Data.json','CM Quaternion and Embedding Data.sage','Integrate Divisors and Recognize Images.sage')]
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
        _g2_progress(out,label,'complete' if result['rows'] and not unresolved else 'unresolved',divisors=len(result['numeric_rows']),unresolved=unresolved)
        if require_complete: require_complete_result(label,result)
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


def require_complete_result(label,result):
    """Keep partial checkpoints, but never report an unresolved result as complete."""
    rows=result.get('rows',[])
    unresolved=[(r.get('D'),r.get('orientation_index'),r.get('operator'),r.get('status'))
                for r in rows if r.get('status') not in ('recognized','identity')]
    if not rows or unresolved:
        raise RuntimeError('%s incomplete; checkpoints preserved; unresolved rows: %s'
                           %(label,unresolved if rows else 'empty result'))

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
    assert len(rows)==len(result['rows']) and len(exact)==len(result['exact_rows'])
    assert set(exact)=={key for key,row in rows.items() if row['status']=='recognized'}
    for key,entry in exact.items():
        ids=[m['orbit_index'] for m in entry['numeric_matches']]
        wanted_ids={r['orbit_index'] for r in result['numeric_rows'] if _twelve_key(r)==key}
        assert len(ids)==len(wanted_ids) and set(ids)==wanted_ids
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
